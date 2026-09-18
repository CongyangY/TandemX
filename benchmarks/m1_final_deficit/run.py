"""Replay frozen fixed-read, edited-assembly M1 decision panel."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
import tempfile

from benchmarks.abundance.simulate import GenomeSpec, build_genome, sample_reads
from benchmarks.m1_shared_signature.full_read_research import FullReadConfig, classify_read
from benchmarks.m1_shared_signature.run_full_read_dev import ordinary_segments
from tandemx.io.sequences import read_sequence_records_many
from tandemx.quantify.mvp import QuantifyConfig, quantify_toy_copy_number

PROTOCOL = Path(__file__).with_name("protocol.json")
COMMITTED_PROTOCOL_SHA256 = "d415717c946a0b5235bc4d7d35c5b4113f5c4ca9f2d5553c3910c06c3586e246"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_derived_estimates(records: list, units: dict[str, str], depth: float,
                           source_genome: str, truth: list[dict]) -> tuple[dict, dict]:
    assigned = {"ordinary_mapping": Counter(), "full_read": Counter()}
    anchors = {
        row["family_id"]: (source_genome[row["start"]-32:row["start"]],
                           source_genome[row["end"]:row["end"]+32])
        for row in truth
    }
    if any(len(left) != 32 or len(right) != 32 for left, right in anchors.values()):
        raise ValueError("Short source flank")
    spanners: dict[str, list[int]] = {name: [] for name in units}
    config = FullReadConfig()
    seen_read_ids: set[str] = set()
    for record in records:
        if record.id in seen_read_ids:
            raise ValueError(f"Duplicate read identifier: {record.id}")
        seen_read_ids.add(record.id)
        sequence = record.sequence.upper()
        for row in ordinary_segments(sequence, units):
            if row["status"] == "assigned":
                assigned["ordinary_mapping"][row["identity"]] += row["end"] - row["start"]
        for row in classify_read(sequence, units, config):
            if row["status"] == "assigned":
                assigned["full_read"][row["identity"]] += row["end"] - row["start"]
        reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for family, (left, right) in anchors.items():
            for oriented in (sequence, reverse):
                first = oriented.find(left)
                last = oriented.find(right)
                if first >= 0 and last >= first + len(left):
                    if oriented.find(left, first+1) >= 0 or oriented.find(right, last+1) >= 0:
                        continue  # nonunique anchor within the read
                    spanners[family].append(last - (first + len(left)))
                    break  # one molecule contributes at most once
    estimates = {
        method: {family: assigned[method][family] / depth for family in units}
        for method in assigned
    }
    estimates["robust_read_span"] = {
        family: statistics.median(values) if len(values) >= 3 else None
        for family, values in spanners.items()
    }
    return estimates, {family: dict(distinct_read_count=len(values),
                                    observed_spans_bp=values) for family, values in spanners.items()}


def run(outdir: Path) -> None:
    protocol_bytes = PROTOCOL.read_bytes()
    if sha(protocol_bytes) != COMMITTED_PROTOCOL_SHA256:
        raise ValueError("Frozen M1 final protocol changed")
    if outdir.exists():
        raise FileExistsError(outdir)
    p = json.loads(protocol_bytes)
    outdir.mkdir(parents=True)
    rows: list[dict] = []
    panels: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="m1_final_") as temp:
        tempbase = Path(temp)
        for rate in p["source_unit_substitution_rates"]:
            for seed in p["source_seeds"]:
                spec = GenomeSpec(seed=seed, periods=tuple(p["periods_bp"]),
                                  copies=tuple(p["source_copies"]), flank_bp=p["flank_bp"],
                                  unit_substitution_rate=rate)
                genome, units, truth = build_genome(spec)
                source_bp = {row["family_id"]: row["repeat_bp"] for row in truth}
                for coverage in p["nominal_coverages"]:
                    panel_id = f"s{seed}_v{int(rate*100):02d}_c{coverage}"
                    base = tempbase / panel_id
                    read_dir = base / "sampling"
                    sample = sample_reads(genome, truth, read_dir,
                                          seed=p["read_seed_offset"]+seed,
                                          coverage=coverage, read_length=p["read_length_bp"],
                                          substitution_rate=p["sequencing_substitution_rate"])
                    records = list(read_sequence_records_many([read_dir / "reads.fa"]))
                    if sum(len(x.sequence) for x in records) != sample["read_count"]*p["read_length_bp"]:
                        raise AssertionError("Sampling receipt mismatch")
                    depth = sample["actual_base_coverage"]
                    read_est, spans = read_derived_estimates(records, units, depth, genome, truth)
                    catalogue = base / "catalogue.fa"
                    catalogue.write_text("".join(f">{key}\n{unit}\n" for key, unit in units.items()))
                    production = quantify_toy_copy_number(QuantifyConfig(
                        reads=read_dir / "reads.fa", monomers=catalogue,
                        genome_size=len(genome), outdir=base / "quantify",
                        k=p["quantify_k"], haploid_depth=None, kmer_backend="python"))
                    if any(abs(row.haploid_depth-depth) > 1e-12 for row in production):
                        raise AssertionError("Methods did not share observed depth")
                    read_est["production_quantify"] = {x.family_id: x.estimated_bp for x in production}
                    warnings = {x.family_id: x.warning for x in production}
                    panels.append(dict(panel_id=panel_id, source_seed=seed, unit_substitution_rate=rate,
                                       nominal_coverage=coverage, actual_coverage=depth,
                                       source_genome_sha256=sha(genome.encode()),
                                       reads_sha256=sample["files"]["reads.fa"],
                                       catalogue_sha256=sha(catalogue.read_bytes()),
                                       read_count=sample["read_count"],
                                       read_side_estimates_bp=read_est, span_support=spans,
                                       production_warnings=warnings))
                    for fraction in p["assembly_copy_fractions"]:
                        assembly, assembly_units, edited = build_genome(spec, fraction)
                        if assembly_units != units:
                            raise AssertionError("Catalogue changed across edit")
                        edited_bp = {row["family_id"]: row["repeat_bp"] for row in edited}
                        for family in units:
                            truth_missing = source_bp[family] - edited_bp[family]
                            for method, by_family in read_est.items():
                                estimate = by_family[family]
                                predicted = None if estimate is None else estimate - edited_bp[family]
                                rows.append(dict(panel_id=panel_id, family=family, method=method,
                                                 assembly_fraction=fraction,
                                                 assembly_sha256=sha(assembly.encode()),
                                                 source_bp=source_bp[family],
                                                 edited_assembly_bp=edited_bp[family],
                                                 injected_missing_bp=truth_missing,
                                                 read_estimate_bp=estimate,
                                                 predicted_signed_deficit_bp=predicted,
                                                 signed_error_bp=None if predicted is None else predicted-truth_missing,
                                                 absolute_error_bp=None if predicted is None else abs(predicted-truth_missing),
                                                 relative_error_to_source=None if predicted is None else abs(predicted-truth_missing)/source_bp[family]))
    results = dict(protocol_sha256=sha(protocol_bytes), protocol_commit="8386133",
                   rows=rows, panels=panels)
    (outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True)+"\n")
    summary = summarize(rows, p)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")


def summarize(rows: list[dict], p: dict) -> dict:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(row["method"], []).append(row)
    methods = {}
    for name, group in grouped.items():
        observed = [r for r in group if r["absolute_error_bp"] is not None]
        eligible = [r for r in group if r["panel_id"].endswith(("_c20", "_c30"))]
        by_family = {family: [r for r in observed if r["family"] == family]
                     for family in ("f1", "f2", "f3")}
        intact = [r for r in eligible if r["assembly_fraction"] == 1]
        a_pass = (len(eligible) == len([r for r in eligible if r["absolute_error_bp"] is not None])
                  and all(r["relative_error_to_source"] <= .10 for r in intact)
                  and all(r["relative_error_to_source"] <= .20 for r in eligible))
        methods[name] = dict(cases=len(group), decision_coverage=len(observed)/len(group),
                             mae_bp=statistics.mean(r["absolute_error_bp"] for r in observed) if observed else None,
                             max_relative_error_to_source=max((r["relative_error_to_source"] for r in observed), default=None),
                             family_mae_bp={f: statistics.mean(r["absolute_error_bp"] for r in family_rows)
                                            if family_rows else None for f, family_rows in by_family.items()},
                             a_gate_pass=a_pass,
                             eligible_intact_baseline_errors=[r["signed_error_bp"] for r in intact])
    return dict(methods=methods, a_gate_pass=methods["production_quantify"]["a_gate_pass"],
                b_gate_assessable_from_fixed_read_edits=False,
                explanation="Within each panel, signed error is the constant read estimate minus source bp for every assembly edit; monotonicity is algebraic. Repeated edits must not be counted as independent evidence.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    run(parser.parse_args().outdir)
