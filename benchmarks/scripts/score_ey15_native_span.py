"""Score the frozen source-guided Ey15 native read-span development diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import re
import statistics
from pathlib import Path

from benchmarks.m2_routes.native_read_pilot import (
    MAX_FLANK_ERROR_FRACTION, fasta_records, fastq_records, reverse_complement,
    trim_between_flanks,
)
from benchmarks.scripts.build_native_read_collapse import sha256_file
from benchmarks.scripts.verify_native_bp_collapse import verify


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/controlled_collapse/ey15_native_pair_audit_20260917"
EDIT = ROOT / "benchmarks/controlled_collapse/ey15_native_bp_development_v1"
PROTOCOL = EDIT / "score_protocol.json"


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Output already exists; preserve previous result")
    protocol = json.loads(PROTOCOL.read_text())
    if (protocol["status"] != "frozen_before_scoring"
            or protocol["case_denominator"] != 9
            or protocol["natural_flank_bp_each_side"] != 1024
            or protocol["flank_edit_distance_fraction_max"] != 0.05
            or MAX_FLANK_ERROR_FRACTION != protocol["flank_edit_distance_fraction_max"]
            or protocol["minimum_eligible_distinct_native_reads"] != 3):
        raise ValueError("Score protocol does not match the frozen implementation")
    manifest_path = SOURCE / "source_eligibility_manifest.json"
    context_path = SOURCE / "ey15_9994_context.fa"
    reads_path = SOURCE / "ey15_native_spanners.fastq.gz"
    verified = verify(manifest_path, context_path, reads_path, EDIT)
    manifest = json.loads(manifest_path.read_text())
    if (manifest["benchmark_eligibility"]["frozen_M2_route"] !=
            "ineligible_420bp_family_period_exceeds_300bp_route_cap"
            or manifest["split"] != "development"):
        raise ValueError("Unexpected Ey15/M2 eligibility")
    source_record = fasta_records(context_path)
    if len(source_record) != 1:
        raise ValueError("Expected one source context")
    source = next(iter(source_record.values()))
    start = manifest["array"]["context_array_start0"]
    end = manifest["array"]["context_array_end0"]
    flank_bp = protocol["natural_flank_bp_each_side"]
    left, right = source[start-flank_bp:start], source[end:end+flank_bp]
    reads = fastq_records(reads_path)
    paf_path = SOURCE / "ey15_native_spanners_context.paf"
    if sha256_file(paf_path) != manifest["artifact_sha256"][paf_path.name]:
        raise ValueError("Native context PAF hash mismatch")
    spanning_context_alignments = {}
    for line in paf_path.read_text().splitlines():
        fields = line.split("\t")
        if (fields[5] != next(iter(source_record)) or int(fields[7]) > start-flank_bp
                or int(fields[8]) < end+flank_bp):
            continue
        if "tp:A:P" not in fields[12:]:
            raise ValueError("Spanning-context PAF alignment is not primary")
        if fields[0] in spanning_context_alignments:
            raise ValueError("Multiple spanning-context alignments for one read")
        cigar_tags = [tag[5:] for tag in fields[12:] if tag.startswith("cg:Z:")]
        if len(cigar_tags) != 1:
            raise ValueError("Missing full-context CIGAR")
        qlen, qstart, qend = map(int, fields[1:4])
        qpos = qstart if fields[4] == "+" else qlen-qend
        rpos = int(fields[7])
        projected = {}
        for count_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_tags[0]):
            n = int(count_text)
            if op in "M=X":
                for boundary in (start, end):
                    if rpos <= boundary <= rpos+n:
                        projected[boundary] = qpos + boundary-rpos
                qpos += n
                rpos += n
            elif op in "DN":
                rpos += n
            elif op == "I":
                qpos += n
            else:
                raise ValueError("Unsupported clipping in full-context PAF CIGAR")
        if set(projected) != {start, end} or rpos != int(fields[8]):
            raise ValueError("Array boundary cannot be projected from full-context CIGAR")
        spanning_context_alignments[fields[0]] = [projected[start], projected[end]]
    expected_ids = {row["read_id"] for row in manifest["per_read"]}
    if len(expected_ids) != 7 or set(reads) != expected_ids or set(spanning_context_alignments) != expected_ids:
        raise ValueError("Native read count or IDs changed")
    trim_rows = []
    for row in sorted(manifest["per_read"], key=lambda item: item["read_id"]):
        read = reads[row["read_id"]]
        if (len(read) != row["read_length_bp"]
                or min(row["left_natural_flank_bp"], row["right_natural_flank_bp"]) < 1000
                or row["whole_assembly_chromosome"] != manifest["array"]["chromosome"]
                or row["whole_assembly_mapq"] < 20):
            raise ValueError("Source read eligibility changed")
        oriented = read if row["strand"] == "+" else reverse_complement(read)
        hit = trim_between_flanks(oriented, left, right)
        hit.pop("array_sequence", None)
        hit.update(read_id=row["read_id"], strand=row["strand"],
                   original_read_sha256=hashlib.sha256(read.encode()).hexdigest(),
                   read_length_bp=len(read),
                   spanning_context_PAF_projected_array_interval_0based=spanning_context_alignments[row["read_id"]])
        if hit["status"] == "ELIGIBLE":
            a, b = hit["trim_interval_0based"]
            if max(abs(a-spanning_context_alignments[row["read_id"]][0]),
                   abs(b-spanning_context_alignments[row["read_id"]][1])) > 10:
                hit["status"] = "ABSTAIN"
                hit["reason"] = "flank_hits_disagree_with_full_context_PAF_CIGAR"
            else:
                hit["observed_array_span_bp"] = b-a
        trim_rows.append(hit)
    eligible_spans = [row["observed_array_span_bp"] for row in trim_rows if row["status"] == "ELIGIBLE"]
    median = statistics.median(eligible_spans) if len(eligible_spans) >= 3 else None
    receipt = json.loads((EDIT / "receipt.json").read_text())
    if len(receipt["cases"]) != 9 or receipt["original_native_fastq_sha256"] != sha256_file(reads_path):
        raise ValueError("Edited case receipt changed")
    scores = []
    for case in receipt["cases"]:
        edited_path = EDIT / f"{case['case_id']}.fa"
        edited_records = fasta_records(edited_path)
        if set(edited_records) != {case["case_id"]}:
            raise ValueError("Unexpected edited FASTA record ID")
        assembly_hit = trim_between_flanks(edited_records[case["case_id"]], left, right)
        if assembly_hit["status"] != "ELIGIBLE":
            raise ValueError(f"Source-guided edited assembly flank trim failed: {case['case_id']}")
        assembly_bp = len(assembly_hit["array_sequence"])
        if (assembly_bp != case["edited_array_bp"] or
                assembly_hit["trim_interval_0based"] != case["edited_array_context_interval"]):
            raise ValueError(f"Edited assembly span disagrees with exact edit ledger: {case['case_id']}")
        truth_positive = case["injected_deleted_bp"] > 0
        if median is None:
            score = {"state": "ABSTAIN", "reason": "fewer_than_three_eligible_native_reads",
                     "signed_bp_delta_assembly_minus_read": None, "threshold_bp": None}
        else:
            delta = assembly_bp - median
            threshold = max(2, math.ceil(max(assembly_bp, median) * 0.05))
            score = {"state": "DISCORDANT" if abs(delta) > threshold else "SUPPORTED",
                     "reason": "fixed_five_percent_span_threshold",
                     "signed_bp_delta_assembly_minus_read": delta, "threshold_bp": threshold}
        scores.append({"case_id": case["case_id"], "injected_deleted_bp": case["injected_deleted_bp"],
                       "truth_positive_injected_bp": truth_positive,
                       "assembly_array_bp": assembly_bp,
                       "assembly_flank_hit_interval_0based": assembly_hit["trim_interval_0based"],
                       "median_native_read_span_bp": median,
                       "length_baseline": score,
                       "frozen_m2_state": "INELIGIBLE_PERIOD_GT_300BP",
                       "biological_audit_state": "NOT_EVALUATED_EXACT_DNA_PAIRING_UNVERIFIED"})
    output.mkdir(parents=True, exist_ok=False)
    (output / "read_trims.json").write_text(json.dumps(trim_rows, indent=2, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as handle:
        for score in scores:
            handle.write(json.dumps(score, sort_keys=True) + "\n")
    tp = sum(row["truth_positive_injected_bp"] and row["length_baseline"]["state"] == "DISCORDANT"
             for row in scores)
    fn = sum(row["truth_positive_injected_bp"] and row["length_baseline"]["state"] == "SUPPORTED"
             for row in scores)
    tn = sum(not row["truth_positive_injected_bp"] and row["length_baseline"]["state"] == "SUPPORTED"
             for row in scores)
    fp = sum(not row["truth_positive_injected_bp"] and row["length_baseline"]["state"] == "DISCORDANT"
             for row in scores)
    abstain = sum(row["length_baseline"]["state"] == "ABSTAIN" for row in scores)
    summary = {"schema_version": 1, "case_count": len(scores), "distinct_source_read_count": len(reads),
               "eligible_read_count": len(eligible_spans), "eligible_read_spans_bp": eligible_spans,
               "median_native_read_span_bp": median, "tp": tp, "fn": fn, "tn": tn, "fp": fp,
               "abstain_cases": abstain, "frozen_m2_evaluated_cases": 0,
               "physical_copy_truth": "unavailable", "biological_accuracy": "not_evaluated",
               "source_manifest_sha256": sha256_file(manifest_path),
               "native_context_PAF_sha256": sha256_file(paf_path),
               "original_native_fastq_sha256": sha256_file(reads_path),
               "edit_receipt_sha256": sha256_file(EDIT / "receipt.json"),
               "score_protocol_sha256": sha256_file(PROTOCOL),
               "scorer_sha256": sha256_file(Path(__file__)),
               "flank_matcher_source_sha256": sha256_file(ROOT / "benchmarks/m2_routes/native_read_pilot.py"),
               "edlib_version": importlib.metadata.version("edlib"),
               "read_trims_sha256": sha256_file(output / "read_trims.json"),
               "per_case_sha256": sha256_file(output / "per_case.jsonl"),
               "verification_receipt_sha256": verified["receipt_sha256"]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
