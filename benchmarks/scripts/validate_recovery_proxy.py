#!/usr/bin/env python3
"""Post-lock proxy assessment; never generates or changes recovery candidates."""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from tandemx.discover.mvp import FastaRecord
from tandemx.io.sequences import read_sequence_records
from tandemx.locate.mvp import build_family_kmer_index, locate_record_arrays
from tandemx.quantify.mvp import read_monomer_fasta
from tandemx.recovery.evidence import extract_intervals, read_paf
from tandemx.recovery.poc import run_mapping, sha256, write_json, write_table


def validate(recovery: Path, newer: Path, outdir: Path) -> None:
    lock_path = recovery / "candidate_lock.json"
    lock = json.loads(lock_path.read_text())
    if lock.get("status") != "candidate_generation_complete":
        raise ValueError("Completed candidate lock required before opening newer assembly")
    required = {"recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed",
                "recruited_reads.tsv", "recovered_sequences.fasta"}
    if set(lock.get("outputs", {})) != required:
        raise ValueError("Candidate lock lacks the canonical output set")
    for name, expected in lock["outputs"].items():
        if sha256(recovery / name) != expected:
            raise ValueError(f"Locked recovery output changed: {name}")
    if outdir.exists():
        raise ValueError("Posthoc validation requires a new output directory")
    prepare_lock = json.loads((recovery / "prepare_lock.json").read_text())
    if (sha256(recovery / "prepared.json") != prepare_lock["prepared_sha256"]
            or sha256(recovery / "input_manifest.json") != prepare_lock["input_manifest_sha256"]):
        raise ValueError("Prepared recovery input lock changed")
    state = json.loads((recovery / "prepared.json").read_text())
    config = state["config"]
    if sha256(recovery / "recruitment_targets.fa") != state["targets_sha256"]:
        raise ValueError("Locked recruitment targets changed")
    manifest = json.loads((recovery / "input_manifest.json").read_text())
    for key in ("old_assembly", "catalog"):
        entry = manifest["inputs"][key]
        if str(Path(config["inputs"][key])) != entry["path"] or sha256(Path(entry["path"])) != entry["sha256"]:
            raise ValueError(f"Prepared source input changed: {key}")
    if sha256(Path(config["minimap2"])) != manifest["minimap2_sha256"]:
        raise ValueError("Prepared mapper changed")
    with (recovery / "recovery_candidates.tsv").open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_read_id"] != "NA"]
    if not rows:
        raise ValueError("No locked candidate requires proxy evaluation")
    outdir.mkdir(parents=True)
    write_json(outdir / "input_manifest.json", {
        "candidate_lock_sha256": sha256(lock_path), "newer_assembly": str(newer),
        "newer_sha256": sha256(newer), "recovery": str(recovery),
        "source_sha256": sha256(Path(__file__)),
        "scope": "post_lock_reference_proxy_not_absolute_copy_truth",
        "frozen_localization": {"k": 21, "min_identity": .9, "identity_model": "iid_base",
                                "catalogue": "complete_frozen_catalogue"},
    })
    flank_sequences = {r.id: r.sequence for r in read_sequence_records(recovery / "recruitment_targets.fa")}
    pairs = {}
    with (outdir / "selected_flanks.fa").open("w") as handle:
        for row in rows:
            pair = []
            for side in ("left", "right"):
                anchors = sorted((a for a in state["anchors"] if a["eligible"]
                                  and a["locus_id"] == row["locus_id"] and a["side"] == side),
                                 key=lambda a: a["offset"])
                anchor = anchors[0]
                pair.append(anchor)
                handle.write(f'>{anchor["anchor_id"]}\n{flank_sequences[anchor["anchor_id"]]}\n')
            pairs[row["locus_id"]] = pair
    paf = outdir / "flanks_to_newer.paf"
    run_mapping([config["minimap2"], "-x", "asm5", "-c", "-k", "15", "-w", "5",
                 "-N", "1000", "-p", "0.5", "--secondary=yes", "-t", "4",
                 str(newer), str(outdir / "selected_flanks.fa")], paf)
    by_anchor = defaultdict(list)
    for hit in read_paf(paf):
        if hit.query_coverage >= .9 and hit.identity >= .95:
            by_anchor[hit.query].append(hit)
    old_intervals, proxy_intervals, proxy_orientation = {}, {}, {}
    assessments = []
    for row in rows:
        locus = row["locus_id"]
        left, right = pairs[locus]
        old_intervals[locus + "_old"] = (row["chromosome"], int(row["old_between_anchor_start"]), int(row["old_between_anchor_end"]))
        lhs, rhs = by_anchor[left["anchor_id"]], by_anchor[right["anchor_id"]]
        record = {"locus_id": locus, "family_id": row["family_id"],
                  "newer_left_hits": len(lhs), "newer_right_hits": len(rhs),
                  "proxy_status": "unresolved_proxy_anchor_placement"}
        if len(lhs) == len(rhs) == 1:
            a, b = lhs[0], rhs[0]
            if a.identity >= .98 and b.identity >= .98 and a.target == b.target and a.strand == b.strand:
                start, end = (a.target_end, b.target_start) if a.strand == "+" else (b.target_end, a.target_start)
                if end > start:
                    proxy_intervals[locus + "_newer"] = (a.target, start, end)
                    proxy_orientation[locus + "_newer"] = a.strand
                    record.update(proxy_status="unique_flank_pair_reference_proxy", newer_chromosome=a.target,
                                  newer_start=start, newer_end=end, newer_span_bp=end-start)
        assessments.append(record)
    def interval_rows(intervals):
        return [dict(anchor_id=name, chromosome=chrom, start=start, end=end)
                for name, (chrom, start, end) in intervals.items()]
    sequences = extract_intervals(Path(config["inputs"]["old_assembly"]), interval_rows(old_intervals))
    sequences.update(extract_intervals(newer, interval_rows(proxy_intervals)))
    for key, strand in proxy_orientation.items():
        if strand == "-":
            sequences[key] = sequences[key].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for record in read_sequence_records(recovery / "recovered_sequences.fasta"):
        sequences[record.id.split(";")[0] + "_candidate"] = record.sequence
    with (outdir / "comparison_spans.fa").open("w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n{sequence}\n")
    monomers = list(read_monomer_fasta(Path(config["inputs"]["catalog"])))
    index, indexed, shared = build_family_kmer_index(monomers, 21)
    array_rows, repeat_bp = [], {}
    family_by_locus = {r["locus_id"]: r["family_id"] for r in rows}
    for name, sequence in sequences.items():
        family = family_by_locus[name.rsplit("_", 1)[0]]
        arrays = locate_record_arrays(FastaRecord(name, name, sequence), monomers, index, indexed, shared, 21, .9, "iid_base")
        matching = [a for a in arrays if a.family_id == family]
        repeat_bp[name] = sum(a.end-a.start for a in matching)
        array_rows.extend(dict(sequence_id=name, start=a.start, end=a.end, family_id=a.family_id,
                               confidence=a.confidence, warning=a.warning) for a in matching)
    for result in assessments:
        locus = result["locus_id"]
        result.update(old_span_bp=len(sequences[locus+"_old"]),
                      candidate_span_bp=len(sequences[locus+"_candidate"]),
                      old_localized_repeat_bp=repeat_bp[locus+"_old"],
                      candidate_localized_repeat_bp=repeat_bp[locus+"_candidate"],
                      candidate_repeat_gain_bp=repeat_bp[locus+"_candidate"]-repeat_bp[locus+"_old"])
        if locus+"_newer" in sequences:
            result["newer_localized_repeat_bp"] = repeat_bp[locus+"_newer"]
            result["old_repeat_distance_to_proxy_bp"] = abs(repeat_bp[locus+"_old"]-repeat_bp[locus+"_newer"])
            result["candidate_repeat_distance_to_proxy_bp"] = abs(repeat_bp[locus+"_candidate"]-repeat_bp[locus+"_newer"])
        result["warning"] = "same_frozen_localizer_not_independent_repeat_truth;candidate_span_includes_intervening_nonrepeat;reference_proxy_not_absolute_truth"
    fields = list(dict.fromkeys(k for row in assessments for k in row))
    write_table(outdir / "proxy_validation.tsv", assessments, fields)
    write_table(outdir / "localized_repeat_intervals.tsv", array_rows,
                ["sequence_id", "start", "end", "family_id", "confidence", "warning"])
    write_json(outdir / "proxy_completion.json", {"status": "posthoc_proxy_complete",
               "candidate_lock_sha256": sha256(lock_path),
               "outputs": {name: sha256(outdir/name) for name in
                           ("proxy_validation.tsv", "localized_repeat_intervals.tsv", "comparison_spans.fa")}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--newer", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    validate(args.recovery, args.newer, args.outdir)
