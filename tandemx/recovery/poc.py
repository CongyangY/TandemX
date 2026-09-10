"""Auditable two-stage recovery PoC, never an in-place assembly patcher."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import logging
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from tandemx.io.sequences import read_sequence_records
from tandemx.recovery.evidence import (
    anchor_uniqueness, extract_intervals, historical_loci, make_flanks,
    read_paf, select_spanning_candidate, spanning_interval,
)
from tandemx.recovery.read_store import collect_recruitment

LOG = logging.getLogger(__name__)
CANDIDATE_FIELDS = (
    "locus_id", "family_id", "chromosome", "start", "end",
    "original_assembly_repeat_bp", "read_derived_abundance_bp", "family_abundance_deficit_bp",
    "left_flank_uniqueness", "right_flank_uniqueness", "recruited_read_count",
    "flank_anchored_read_count", "dual_flank_read_count", "maximum_read_span_bp",
    "recovered_bp", "recovery_status", "confidence", "failure_reason", "warning",
    "candidate_span_bp", "maximum_recruited_read_length_bp", "candidate_read_id",
    "old_between_anchor_start", "old_between_anchor_end",
)
FROZEN_RULES = {
    "flank_bp": 2000, "flank_offsets_bp": [0, 2000, 5000, 10000],
    "anchor_query_coverage_min": .9, "anchor_identity_min": .98,
    "competing_hit_identity_min": .95, "repeat_read_aligned_bp_min": 500,
    "repeat_read_identity_min": .85, "minimum_dual_anchor_reads": 3,
    "span_length_tolerance": "max(100bp,1pct_median)",
    "candidate_method": "median_length_observed_dual_anchor_read_span_unpolished",
    "mapping_timeout_seconds": 1800, "maximum_recruited_ids": 100000,
    "maximum_retained_alignment_rows": 2000000, "stop_on_no_supported_candidates": True,
}


def validate_rules(config: dict) -> None:
    rules = config.get("rules", FROZEN_RULES)
    resource_fields = {"maximum_recruited_ids", "maximum_retained_alignment_rows"}
    if ({k:v for k,v in rules.items() if k not in resource_fields}
            != {k:v for k,v in FROZEN_RULES.items() if k not in resource_fields}):
        raise ValueError("Recovery config rules differ from the frozen implemented contract")
    for key, maximum in (("maximum_recruited_ids", 1000000), ("maximum_retained_alignment_rows", 10000000)):
        if type(rules.get(key)) is not int or not 1 <= rules[key] <= maximum:
            raise ValueError(f"Invalid bounded resource rule: {key}")
    size = config.get("recruitment_template_min_bp", 30000)
    if type(size) is not int or not 500 <= size <= 30000:
        raise ValueError("Recruitment template must be 500--30000 bp")


def source_identity() -> dict:
    root = Path(__file__).parent
    return {path.name: sha256(path) for path in sorted(root.glob("*.py"))}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")


def write_table(path: Path, rows: list[dict], fields: tuple | list | None = None) -> None:
    columns = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row.get(key) for key in columns})


def run_mapping(command: list[str], output: Path, timeout: int = 1800) -> None:
    """Require an exact command/output receipt to reuse an external stage."""
    receipt = output.with_suffix(".receipt.json")
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if (saved.get("command") == command and saved.get("exit_code") == 0
                and output.is_file() and saved.get("sha256") == sha256(output)):
            LOG.info("Reusing hash-verified mapping %s", output.name)
            return
        raise ValueError(f"Existing mapping receipt is not reusable: {receipt}")
    partial = output.with_suffix(".partial.paf")
    if partial.exists() or output.exists():
        raise ValueError(f"Unreceipted mapping exists; retain it and use a fresh directory: {output}")
    start = time.monotonic()
    LOG.info("Starting mapping %s", output.name)
    with partial.open("w") as stdout, output.with_suffix(".log").open("w") as stderr:
        try:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
        except subprocess.TimeoutExpired:
            write_json(receipt, dict(command=command, exit_code=None, fate="timeout",
                                    elapsed_seconds=time.monotonic()-start))
            raise
    payload = dict(command=command, exit_code=result.returncode, elapsed_seconds=time.monotonic()-start)
    if result.returncode != 0:
        write_json(receipt, dict(payload, fate="external_mapping_failure"))
        raise ValueError(f"Mapping failed ({result.returncode}); see {output.with_suffix('.log')}")
    partial.rename(output)
    write_json(receipt, dict(payload, sha256=sha256(output), bytes=output.stat().st_size))


def load_catalog(path: Path, families: set[str]) -> dict[str, str]:
    selected = {}
    for record in read_sequence_records(path):
        fields = dict(part.split("=", 1) for part in record.description.split(";") if "=" in part)
        family = fields.get("family_id", record.id)
        if family in families:
            if family in selected:
                raise ValueError(f"Duplicate family {family}")
            selected[family] = record.sequence
    if set(selected) != families:
        raise ValueError("Selected family missing from frozen catalogue")
    return selected


def write_targets(path: Path, catalog: dict[str, str], anchors: list[dict],
                  sequences: dict[str, str], minimum_bp: int) -> None:
    with path.open("w") as handle:
        for family, monomer in sorted(catalog.items()):
            handle.write(f">repeat_{family}\n{monomer*((minimum_bp+len(monomer)-1)//len(monomer))}\n")
        for anchor in anchors:
            if anchor["eligible"]:
                handle.write(f'>{anchor["anchor_id"]}\n{sequences[anchor["anchor_id"]]}\n')


def reuse_preparation(config_path: Path, outdir: Path, previous: Path) -> dict:
    """Import a locked old-flank audit with template/resource-only changes."""
    config = json.loads(config_path.read_text())
    validate_rules(config)
    state = json.loads((previous / "prepared.json").read_text())
    lock = json.loads((previous / "prepare_lock.json").read_text())
    if outdir.exists():
        raise ValueError("Reused preparation requires a new output directory")
    if (sha256(previous / "prepared.json") != lock["prepared_sha256"]
            or sha256(previous / "input_manifest.json") != lock["input_manifest_sha256"]
            or sha256(previous / "flank_audit.tsv") != state["anchor_audit_sha256"]
            or sha256(previous / "recruitment_targets.fa") != state["targets_sha256"]):
        raise ValueError("Previous locked preparation changed")
    exceptions = {"experiment_id", "recruitment_template_min_bp", "continuation_reason", "prior_attempt", "rules"}
    if {k:v for k,v in config.items() if k not in exceptions} != {k:v for k,v in state["config"].items() if k not in exceptions}:
        raise ValueError("Only recruitment template length, resource caps and continuation metadata may change")
    validate_rules(state["config"])
    catalog = load_catalog(Path(config["inputs"]["catalog"]), set(config["families"]))
    sequences = {r.id:r.sequence for r in read_sequence_records(previous / "flanks.fa")}
    outdir.mkdir(parents=True)
    for name in ("flanks.fa", "flanks_to_old.paf", "flanks_to_old.log", "flanks_to_old.receipt.json", "flank_audit.tsv"):
        shutil.copy2(previous / name, outdir / name)
    manifest = json.loads((previous / "input_manifest.json").read_text())
    manifest.update(config_sha256=sha256(config_path), imported_preparation=str(previous.resolve()),
                    imported_preparation_lock_sha256=sha256(previous / "prepare_lock.json"))
    write_json(outdir / "input_manifest.json", manifest)
    write_json(outdir / "run_config.yaml", config)
    write_targets(outdir / "recruitment_targets.fa", catalog, state["anchors"], sequences,
                  config.get("recruitment_template_min_bp", 30000))
    state.update(config=config, targets_sha256=sha256(outdir / "recruitment_targets.fa"))
    write_json(outdir / "prepared.json", state)
    write_json(outdir / "prepare_lock.json", dict(prepared_sha256=sha256(outdir / "prepared.json"),
               input_manifest_sha256=sha256(outdir / "input_manifest.json"), source=source_identity(),
               flank_audit_imported_without_rerun=True))
    return state


def prepare(config_path: Path, outdir: Path) -> dict:
    config = json.loads(config_path.read_text())
    validate_rules(config)
    config.setdefault("rules", FROZEN_RULES)
    if config.get("schema_version") != 1 or not config.get("families"):
        raise ValueError("Invalid recovery enrollment")
    if any("new" in key.lower() for key in config.get("inputs", {})):
        raise ValueError("Newer-assembly inputs are prohibited during recovery")
    required = {"old_assembly", "old_arrays", "catalog", "copy_number", "reads", "read_qc"}
    if set(config["inputs"]) != required:
        raise ValueError(f"Recovery requires exactly {sorted(required)}")
    if outdir.exists():
        raise ValueError("Recovery preparation requires a fresh directory")
    paths = {key: Path(value) for key, value in config["inputs"].items()}
    if any(not path.is_file() for path in paths.values()):
        raise ValueError("An enrolled recovery input is missing")
    minimap = Path(config["minimap2"])
    if not minimap.is_file():
        raise ValueError("Pinned minimap2 executable is missing")
    qc = json.loads(paths["read_qc"].read_text())
    if (qc.get("complete") is not True or qc.get("fastq_records_valid") is not True
            or qc.get("exact_duplicate_read_ids") != 0):
        raise ValueError("Full-read QC has not passed")
    inputs = {}
    for key, path in paths.items():
        LOG.info("Verifying %s", key)
        digest = sha256(path)
        expected = config.get("expected_sha256", {}).get(key)
        if (expected and digest != expected) or (key == "reads" and digest != qc["input_sha256"]):
            raise ValueError(f"Frozen input hash mismatch: {key}")
        inputs[key] = dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=digest)
    outdir.mkdir(parents=True)
    write_json(outdir / "run_config.yaml", config)  # JSON is a YAML subset.
    write_json(outdir / "input_manifest.json", dict(inputs=inputs, config_sha256=sha256(config_path),
               minimap2_sha256=sha256(minimap), minimap2_version=subprocess.check_output([str(minimap), "--version"], text=True).strip(),
               recovery_input_separation="newer_sequence_and_locus_coordinates_not_used",
               prior_knowledge="family_level_prior_results_known_not_fully_unobserved_validation"))
    catalog = load_catalog(paths["catalog"], set(config["families"]))
    loci = historical_loci(paths["old_arrays"], set(catalog))
    flanks = make_flanks(loci)
    sequences = extract_intervals(paths["old_assembly"], flanks)
    bank = outdir / "flanks.fa"
    with bank.open("w") as handle:
        for row in flanks:
            sequence = sequences.get(row["anchor_id"], "")
            if len(sequence) == 2000 and not set(sequence)-set("ACGT"):
                handle.write(f'>{row["anchor_id"]}\n{sequence}\n')
    mapping = outdir / "flanks_to_old.paf"
    if bank.stat().st_size:
        run_mapping([str(minimap), "-x", "asm5", "-c", "-k", "15", "-w", "5", "-N", "1000", "-p", "0.5",
                     "--secondary=yes", "-t", str(config.get("threads", 4)), str(paths["old_assembly"]), str(bank)], mapping)
    else:
        mapping.write_text("")
    hits = defaultdict(list)
    for hit in read_paf(mapping):
        hits[hit.query].append(hit)
    anchors = [anchor_uniqueness(row, sequences.get(row["anchor_id"], ""), hits[row["anchor_id"]]) for row in flanks]
    write_table(outdir / "flank_audit.tsv", anchors, ["anchor_id", "locus_id", "family_id", "chromosome", "start", "end", "side", "offset", "sequence_length", "eligible", "qualifying_hits", "uniqueness", "reason"])
    targets = outdir / "recruitment_targets.fa"
    write_targets(targets, catalog, anchors, sequences, config.get("recruitment_template_min_bp", 30000))
    state = dict(config=config, loci=loci, anchors=anchors, catalog_lengths={key: len(value) for key,value in catalog.items()},
                 targets_sha256=sha256(targets), anchor_audit_sha256=sha256(outdir / "flank_audit.tsv"))
    write_json(outdir / "prepared.json", state)
    write_json(outdir / "prepare_lock.json", dict(prepared_sha256=sha256(outdir / "prepared.json"),
               input_manifest_sha256=sha256(outdir / "input_manifest.json"), source=source_identity()))
    LOG.info("Prepared %s loci, %s eligible anchors", len(loci), sum(row["eligible"] for row in anchors))
    return state


def recruit(outdir: Path, reuse_mapping_from: Path | None = None) -> None:
    state = json.loads((outdir / "prepared.json").read_text())
    config = state["config"]
    validate_rules(config)
    lock = json.loads((outdir / "prepare_lock.json").read_text())
    if (sha256(outdir / "prepared.json") != lock["prepared_sha256"]
            or sha256(outdir / "input_manifest.json") != lock["input_manifest_sha256"]
            or sha256(outdir / "flank_audit.tsv") != state["anchor_audit_sha256"]):
        raise ValueError("Prepared recovery evidence changed")
    targets = outdir / "recruitment_targets.fa"
    if sha256(targets) != state["targets_sha256"]:
        raise ValueError("Recruitment target hash changed")
    manifest = json.loads((outdir / "input_manifest.json").read_text())
    if sha256(Path(config["minimap2"])) != manifest["minimap2_sha256"]:
        raise ValueError("Minimap2 executable changed")
    read_path = Path(config["inputs"]["reads"])
    for key, entry in manifest["inputs"].items():
        path = Path(entry["path"])
        LOG.info("Rechecking enrolled %s", key)
        if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            raise ValueError(f"Recovery input changed since preparation: {key}")
    write_json(outdir / "recruit_source.json", dict(source=source_identity(),
               prepared_sha256=sha256(outdir / "prepared.json")))
    paf = outdir / "reads_to_targets.paf"
    command = [config["minimap2"], "-x", "map-hifi", "-c", "-N", "50", "-p", "0.5", "--secondary=yes",
               "-K", "50M", "-t", str(config.get("threads", 4)), str(targets), str(read_path)]
    if reuse_mapping_from is None:
        run_mapping(command, paf)
    else:
        previous = reuse_mapping_from
        prior_paf = previous / "reads_to_targets.paf"
        prior_receipt = json.loads((previous / "reads_to_targets.receipt.json").read_text())
        prior_manifest = json.loads((previous / "input_manifest.json").read_text())
        equivalent = list(command)
        equivalent[-2] = str(previous / "recruitment_targets.fa")
        if (prior_receipt.get("exit_code") != 0 or prior_receipt.get("command") != equivalent
                or sha256(prior_paf) != prior_receipt.get("sha256")
                or sha256(previous / "recruitment_targets.fa") != sha256(targets)
                or prior_manifest["inputs"] != manifest["inputs"]
                or prior_manifest["minimap2_sha256"] != manifest["minimap2_sha256"]):
            raise ValueError("Previous mapping is not exact-input/command equivalent")
        if paf.exists():
            raise ValueError("Imported mapping destination already exists")
        for suffix in (".paf", ".receipt.json", ".log"):
            shutil.copy2(previous / f"reads_to_targets{suffix}", outdir / f"reads_to_targets{suffix}")
        write_json(outdir / "mapping_import.json", dict(source=str(previous), output_sha256=sha256(paf),
                   receipt_sha256=sha256(previous / "reads_to_targets.receipt.json"),
                   mapping_rerun=False, source_command=prior_receipt["command"],
                   continuation_scope="streamed_collection_with_larger_resource_caps_only"))
    evaluate(outdir, state, paf)


def evaluate(outdir: Path, state: dict, paf: Path) -> None:
    config, loci = state["config"], state["loci"]
    anchors = {row["anchor_id"]: row for row in state["anchors"] if row["eligible"]}
    limits = config.get("rules", FROZEN_RULES)
    repeats, hit_map, read_lengths, retained_rows = collect_recruitment(
        paf, anchors, outdir / "recruited_reads.tsv", set(config["families"]),
        limits["maximum_recruited_ids"], limits["maximum_retained_alignment_rows"])
    repeat_reads = defaultdict(set, repeats)
    anchor_hits = defaultdict(dict, hit_map)
    with Path(config["inputs"]["copy_number"]).open() as handle:
        abundance = {r["family_id"]: float(r["estimated_bp"]) for r in csv.DictReader(handle, delimiter="\t") if r["family_id"] in config["families"]}
    if set(abundance) != set(config["families"]):
        raise ValueError("A selected family is absent from the frozen abundance table")
    family_old = {family: sum(l["end"]-l["start"] for l in loci if l["family_id"] == family) for family in config["families"]}
    results, selected, validation = [], {}, []
    for locus in loci:
        family, locus_id = locus["family_id"], locus["locus_id"]
        available = [a for a in anchors.values() if a["locus_id"] == locus_id]
        left = sorted((a for a in available if a["side"] == "left"), key=lambda a: a["offset"])
        right = sorted((a for a in available if a["side"] == "right"), key=lambda a: a["offset"])
        spans = []
        # Anchor pair is selected from old assembly alone, never from read outcomes.
        if left and right:
            for read_id, mappings in anchor_hits[locus_id].items():
                left_hits = [h for h in mappings if h.target == left[0]["anchor_id"]]
                right_hits = [h for h in mappings if h.target == right[0]["anchor_id"]]
                if read_id not in repeat_reads[family] or len(left_hits) != 1 or len(right_hits) != 1:
                    continue
                interval = spanning_interval(left_hits[0], right_hits[0])
                if interval:
                    spans.append(dict(read_id=read_id, start=interval[0], end=interval[1], strand=interval[2]))
            chosen, status = select_spanning_candidate(spans)
        else:
            chosen, status = None, "unresolved_no_unique_anchor"
        if chosen:
            selected[locus_id] = chosen
        supporting = repeat_reads[family] | set(anchor_hits[locus_id])
        reason = {"partially_resolved": "unpolished_observed_read_span_not_complete_assembly",
                  "unresolved_no_unique_anchor": "no_two_sided_old_assembly_unique_anchor_pair",
                  "insufficient_read_support": "fewer_than_three_unambiguous_repeat_supporting_dual_anchor_reads",
                  "unresolved_conflicting_paths": "dual_anchor_read_span_lengths_disagree"}[status]
        row = dict(locus_id=locus_id, family_id=family, chromosome=locus["chromosome"], start=locus["start"], end=locus["end"],
                   original_assembly_repeat_bp=locus["end"]-locus["start"], read_derived_abundance_bp=abundance[family],
                   family_abundance_deficit_bp=max(abundance[family]-family_old[family], 0),
                   left_flank_uniqueness="unique_in_old_assembly" if left else "unresolved",
                   right_flank_uniqueness="unique_in_old_assembly" if right else "unresolved",
                   recruited_read_count=len(repeat_reads[family]), flank_anchored_read_count=len(anchor_hits[locus_id]),
                   dual_flank_read_count=len(spans), maximum_read_span_bp=max((r["end"]-r["start"] for r in spans), default=0),
                   recovered_bp=None, candidate_span_bp=chosen["end"]-chosen["start"] if chosen else None,
                   maximum_recruited_read_length_bp=max((read_lengths[r] for r in supporting), default=0),
                   candidate_read_id=chosen["read_id"] if chosen else None,
                   old_between_anchor_start=left[0]["end"] if left else None,
                   old_between_anchor_end=right[0]["start"] if right else None, recovery_status=status,
                   confidence="candidate" if chosen else "unresolved", failure_reason=reason,
                   warning="family_abundance_not_locus_specific;assembly_relative_anchor_uniqueness;no_newer_sequence_used;unpolished_span_includes_intervening_flank_sequence")
        results.append(row)
        validation.append(dict(locus_id=locus_id, family_id=family, newer_proxy_evaluated=False,
                               candidate_sequence_id=locus_id if chosen else None, independent_read_support=max(0,len(spans)-1),
                               recovery_status=status, validation_status="candidate_remapping_required" if chosen else "no_candidate_sequence"))
    for family in config["families"]:
        if not any(l["family_id"] == family for l in loci):
            results.append(dict(locus_id=f"{family}_no_locus", family_id=family, chromosome=None, start=None, end=None,
                                original_assembly_repeat_bp=0, read_derived_abundance_bp=abundance[family],
                                family_abundance_deficit_bp=abundance[family], left_flank_uniqueness="unavailable",
                                right_flank_uniqueness="unavailable", recruited_read_count=len(repeat_reads[family]),
                                flank_anchored_read_count=0, dual_flank_read_count=0,
                                maximum_read_span_bp=0,
                                maximum_recruited_read_length_bp=max((read_lengths[r] for r in repeat_reads[family]), default=0),
                                recovered_bp=None, recovery_status="unresolved_no_assembly_locus", confidence="unresolved",
                                failure_reason="frozen_old_assembly_localization_has_no_interval",
                                warning="not_proof_family_absent;family_abundance_not_locus_specific;no_newer_sequence_used"))
            validation.append(dict(locus_id=f"{family}_no_locus", family_id=family,
                                   newer_proxy_evaluated=False, candidate_sequence_id=None,
                                   independent_read_support=0,
                                   recovery_status="unresolved_no_assembly_locus",
                                   validation_status="no_candidate_sequence"))
    write_table(outdir / "recovery_candidates.tsv", results, CANDIDATE_FIELDS)
    write_table(outdir / "recovery_validation.tsv", validation, ["locus_id", "family_id", "newer_proxy_evaluated", "candidate_sequence_id", "independent_read_support", "recovery_status", "validation_status"])
    with (outdir / "recovery_loci.bed").open("w") as handle:
        handle.write("#chrom\tstart\tend\tlocus_id\tfamily_id\tstatus\n")
        for row in results:
            if row["chromosome"] is not None:
                handle.write("\t".join(str(row[k]) for k in ("chromosome", "start", "end", "locus_id", "family_id", "recovery_status"))+"\n")
    written = set()
    with (outdir / "recovered_sequences.fasta").open("w") as handle:
        if selected:
            for record in read_sequence_records(Path(config["inputs"]["reads"])):
                for locus_id, span in selected.items():
                    if record.id == span["read_id"]:
                        sequence = record.sequence[span["start"]:span["end"]]
                        if span["strand"] == "-":
                            sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                        handle.write(f">{locus_id};source_read={record.id};status=partially_resolved;unpolished=true\n{sequence}\n")
                        written.add(locus_id)
    if written != set(selected):
        raise ValueError("Selected candidate source read absent from original input; no candidate lock written")
    counts = defaultdict(int)
    for row in results:
        counts[row["recovery_status"]] += 1
    receipt = dict(status="candidate_generation_complete", selected_families=config["families"],
                   candidate_sequences=len(selected), recovery_status_counts=dict(counts),
                   expansion_decision="requires_remapping_and_posthoc_proxy_validation" if selected else "stop_recovery_expansion_negative_poc",
                   newer_assembly_used=False, complete_genomic_resolvability_test=False,
                   source=source_identity(), candidate_repeat_bp_status="requires_separate_repeat_localization",
                   retained_alignment_rows=retained_rows, unique_recruited_reads=len(read_lengths),
                   resource_caps={k:limits[k] for k in ("maximum_recruited_ids", "maximum_retained_alignment_rows")},
                   scope="fixed_historical_loci_2kb_flanks_at_0_2_5_10kb_and_full_HiFi_read_mapping",
                   outputs={name: sha256(outdir/name) for name in ("recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta")})
    write_json(outdir / "candidate_lock.json", receipt)
    table = "".join("<tr>"+"".join(f"<td>{html.escape(str(row.get(k, 'NA')))}</td>" for k in ("family_id", "locus_id", "recovery_status", "recruited_read_count", "flank_anchored_read_count", "dual_flank_read_count", "failure_reason"))+"</tr>" for row in results)
    (outdir / "recovery_report.html").write_text(f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Targeted recovery evidence</title><style>body{{font:16px system-ui;max-width:1200px;margin:3rem auto;padding:1rem}}table{{border-collapse:collapse}}td,th{{padding:.7rem;border-bottom:1px solid #ccc;text-align:left}}.note{{max-width:80ch}}</style><h1>Targeted repeat recovery</h1><p class="note">{len(selected)} candidate sequences from {len(config["families"])} selected families. No newer-assembly sequence was used. These bounded results do not establish that a family is genomically unresolvable. Family-wide read abundance is not a locus-specific missing-length estimate.</p><table><tr><th>Family</th><th>Locus</th><th>Status</th><th>Repeat reads</th><th>Flank reads</th><th>Dual-anchor reads</th><th>Reason</th></tr>{table}</table><p><a href="recovery_candidates.tsv">Candidate evidence</a> · <a href="flank_audit.tsv">All flank tests</a> · <a href="recruited_reads.tsv">Read evidence</a> · <a href="candidate_lock.json">Candidate lock and decision</a></p></html>')
    LOG.info("Recovery outcomes: %s", dict(counts))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("prepare", "recruit"))
    parser.add_argument("--config", type=Path)
    parser.add_argument("--reuse-preparation", type=Path, help="Import a locked flank audit; only template length and resource caps may differ")
    parser.add_argument("--reuse-mapping-from", type=Path, help="Import an exact-input, completed mapping without rerunning it")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.phase == "prepare":
        if not args.config:
            parser.error("prepare requires --config")
        if args.reuse_preparation:
            reuse_preparation(args.config, args.outdir, args.reuse_preparation)
        else:
            prepare(args.config, args.outdir)
    else:
        recruit(args.outdir, args.reuse_mapping_from)


if __name__ == "__main__":
    main()
