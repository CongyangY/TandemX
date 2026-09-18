"""Score frozen Macadamia injected-bp edits against original record spans."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

from benchmarks.m2_routes.native_read_pilot import (
    MAX_FLANK_ERROR_FRACTION, fasta_records, fastq_records,
    reverse_complement, trim_between_flanks,
)
from benchmarks.scripts.build_native_read_collapse import sha256_file
from benchmarks.scripts.verify_native_bp_collapse import verify


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/controlled_collapse/macadamia_assembly_candidate_v1_20260918"
MAPPING = ROOT / "benchmarks/controlled_collapse/macadamia_native_mapping_v1_20260918"
EDIT = ROOT / "benchmarks/controlled_collapse/macadamia_bp_provisional_v1"


def projected_boundaries(line: str, start: int, end: int) -> list[int]:
    fields = line.rstrip("\n").split("\t")
    cigars = [tag[5:] for tag in fields[12:] if tag.startswith("cg:Z:")]
    if len(cigars) != 1 or "tp:A:P" not in fields[12:]:
        raise ValueError("Expected one primary CIGAR")
    qlen, qstart, qend = map(int, fields[1:4])
    qpos = qstart if fields[4] == "+" else qlen - qend
    rpos = int(fields[7])
    projected = {}
    for text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigars[0]):
        length = int(text)
        if op in "M=X":
            for boundary in (start, end):
                if rpos <= boundary <= rpos + length:
                    projected[boundary] = qpos + boundary - rpos
            qpos += length
            rpos += length
        elif op in "DN":
            rpos += length
        elif op == "I":
            qpos += length
        else:
            raise ValueError("Unsupported clipping in context PAF CIGAR")
    if set(projected) != {start, end} or rpos != int(fields[8]):
        raise ValueError("Array boundary cannot be projected from context PAF")
    return [projected[start], projected[end]]


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Preserve prior score output")
    protocol_path = EDIT / "score_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if (protocol["status"] != "frozen_before_scoring" or
            protocol["case_denominator"] != 9 or
            protocol["natural_flank_bp_each_side"] != 1024 or
            protocol["flank_edit_distance_fraction_max"] != MAX_FLANK_ERROR_FRACTION or
            protocol["minimum_eligible_distinct_records"] != 3 or
            protocol["length_discordance_fraction"] != 0.05):
        raise ValueError("Scorer does not match frozen protocol")
    manifest_path = EDIT / "source_manifest.json"
    context_path = SOURCE / "selected_context.fa"
    reads_path = EDIT / "selected_original_reads.fastq.gz"
    paf_path = MAPPING / "score/context_alignments.paf.gz"
    receipt_path = EDIT / "generated/receipt.json"
    for path, key in ((context_path, "source_context_fasta_sha256"),
                      (reads_path, "original_selected_fastq_sha256"),
                      (paf_path, "context_paf_gzip_sha256"),
                      (receipt_path, "edit_receipt_sha256")):
        if sha256_file(path) != protocol[key]:
            raise ValueError(f"Frozen input changed: {path}")
    verified = verify(manifest_path, context_path, reads_path, EDIT / "generated")
    manifest = json.loads(manifest_path.read_text())
    if (manifest["status"] != "development_provisional_record_support" or
            manifest["support_tier"] != "full_reference_records_zmw_unverified"):
        raise ValueError("Provisional source tier changed")
    sources = fasta_records(context_path)
    if len(sources) != 1:
        raise ValueError("Expected one context")
    context_id, source = next(iter(sources.items()))
    start, end = protocol["array_start_0_in_context"], protocol["array_end_0_in_context"]
    if end - start != manifest["array"]["length_bp"]:
        raise ValueError("Array length changed")
    flank = protocol["natural_flank_bp_each_side"]
    left, right = source[start - flank:start], source[end:end + flank]
    reads = fastq_records(reads_path)
    if len(reads) != 7:
        raise ValueError("Selected source record denominator changed")
    projected = {}
    strands = {}
    with gzip.open(paf_path, "rt") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            read_id = fields[0]
            if read_id not in reads or "tp:A:P" not in fields[12:]:
                continue
            if fields[5] != context_id or read_id in projected:
                raise ValueError("Selected read has wrong or multiple primary context alignment")
            projected[read_id] = projected_boundaries(line, start, end)
            strands[read_id] = fields[4]
    if set(projected) != set(reads):
        raise ValueError("Selected read lacks primary context CIGAR")
    trims = []
    for read_id in sorted(reads):
        read = reads[read_id]
        oriented = read if strands[read_id] == "+" else reverse_complement(read)
        hit = trim_between_flanks(oriented, left, right)
        hit.pop("array_sequence", None)
        hit.update(read_id=read_id, strand=strands[read_id], read_length_bp=len(read),
                   original_read_sequence_sha256=hashlib.sha256(read.encode()).hexdigest(),
                   primary_PAF_projected_interval_0based=projected[read_id])
        if hit["status"] == "ELIGIBLE":
            observed = hit["trim_interval_0based"]
            if max(abs(a - b) for a, b in zip(observed, projected[read_id])) > \
                    protocol["maximum_flank_vs_cigar_boundary_difference_bp"]:
                hit["status"] = "ABSTAIN"
                hit["reason"] = "flank_hits_disagree_with_context_PAF_CIGAR"
            else:
                hit["observed_array_span_bp"] = observed[1] - observed[0]
        trims.append(hit)
    eligible = [row["observed_array_span_bp"] for row in trims if row["status"] == "ELIGIBLE"]
    median = statistics.median(eligible) if len(eligible) >= 3 else None
    receipt = json.loads(receipt_path.read_text())
    cases = []
    for case in receipt["cases"]:
        fasta = EDIT / "generated" / f"{case['case_id']}.fa"
        entries = fasta_records(fasta)
        if set(entries) != {case["case_id"]}:
            raise ValueError("Edited FASTA record ID changed")
        hit = trim_between_flanks(entries[case["case_id"]], left, right)
        assembly_bp = None
        if hit["status"] == "ELIGIBLE":
            assembly_bp = len(hit["array_sequence"])
            if (assembly_bp != case["edited_array_bp"] or
                    hit["trim_interval_0based"] != case["edited_array_context_interval"]):
                raise ValueError(f"Direct edited FASTA span conflicts with ledger: {case['case_id']}")
        if assembly_bp is None:
            state = "ABSTAIN"
            reason = "edited_assembly_natural_flanks_not_unique"
            delta = threshold = None
        elif median is None:
            state = "ABSTAIN"
            reason = "fewer_than_three_eligible_original_records"
            delta = threshold = None
        else:
            delta = assembly_bp - median
            threshold = max(2, math.ceil(max(assembly_bp, median) * 0.05))
            state = "DISCORDANT" if abs(delta) > threshold else "SUPPORTED"
            reason = "fixed_five_percent_span_threshold"
        cases.append({"case_id": case["case_id"], "injected_deleted_bp": case["injected_deleted_bp"],
                      "truth_positive_injected_bp": case["injected_deleted_bp"] > 0,
                      "direct_assembly_span_bp": assembly_bp, "median_original_record_span_bp": median,
                      "signed_delta_assembly_minus_read_bp": delta, "threshold_bp": threshold,
                      "state": state, "reason": reason})
    output.mkdir(parents=True, exist_ok=False)
    (output / "read_trims.json").write_text(json.dumps(trims, indent=2, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as stream:
        for row in cases:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    result = {
        "case_count": len(cases), "selected_original_record_count": len(reads),
        "eligible_record_count": len(eligible), "eligible_spans_bp": eligible,
        "median_original_record_span_bp": median,
        "tp": sum(row["truth_positive_injected_bp"] and row["state"] == "DISCORDANT" for row in cases),
        "fn": sum(row["truth_positive_injected_bp"] and row["state"] == "SUPPORTED" for row in cases),
        "tn": sum(not row["truth_positive_injected_bp"] and row["state"] == "SUPPORTED" for row in cases),
        "fp": sum(not row["truth_positive_injected_bp"] and row["state"] == "DISCORDANT" for row in cases),
        "abstain_cases": sum(row["state"] == "ABSTAIN" for row in cases),
        "molecule_identity": "unverified_ZMW", "physical_copy_truth": "unavailable",
        "biological_accuracy": "not_evaluated", "source_manifest_sha256": sha256_file(manifest_path),
        "score_protocol_sha256": sha256_file(protocol_path),
        "verified_edit_receipt_sha256": verified["receipt_sha256"],
        "read_trims_sha256": sha256_file(output / "read_trims.json"),
        "per_case_sha256": sha256_file(output / "per_case.jsonl"),
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
