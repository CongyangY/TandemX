"""Score the frozen C3 equal-length development challenge without M2 tuning."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any

from benchmarks.m2_routes.native_read_pilot import (
    B1,
    EDIT,
    FLANK_BP,
    ROOT,
    SOURCE,
    fasta_records,
    fastq_records,
    path_summary,
    primary_paf_array_interval,
    reverse_complement,
    sha256,
    trim_between_flanks,
)


CHALLENGE = ROOT / "benchmarks/controlled_collapse/native_equal_length_development_v1"
PILOT = ROOT / "benchmarks/m2_routes/native_read_evidence_20260917"
PROTOTYPE = ROOT / "benchmarks/m2_routes/alignment/prototype.py"
GENERATOR = ROOT / "benchmarks/scripts/build_native_equal_length_challenge.py"


def _assert_hash(path: Path, expected: str) -> None:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(f"SHA-256 mismatch: {path}: {observed} != {expected}")


def _read_fixed_inputs() -> tuple[dict[str, Any], list[dict[str, Any]],
                                  dict[str, str], dict[str, Any]]:
    protocol_path = CHALLENGE / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    receipt_path = CHALLENGE / "receipt.json"
    _assert_hash(receipt_path, protocol["generated_receipt_sha256"])
    _assert_hash(GENERATOR, protocol["generator_sha256"])
    _assert_hash(PROTOTYPE, protocol["frozen_m2_prototype_sha256"])
    _assert_hash(ROOT / "benchmarks/m2_routes/native_read_pilot.py",
                 protocol["native_trim_pilot_sha256"])
    _assert_hash(EDIT / "source_config.json", protocol["source_config_sha256"])
    _assert_hash(SOURCE / "colcen_native_spanners.fastq.gz", protocol["native_read_fastq_sha256"])
    if protocol["registered_before_scoring"] is not True or protocol["case_denominator"] != 6:
        raise ValueError("unexpected frozen challenge protocol")
    if protocol["pairing_status"] != "lineage_context_supported_donor_unverified":
        raise ValueError("biological pairing unexpectedly changed")
    receipt = json.loads(receipt_path.read_text())
    if len(receipt["cases"]) != 6 or receipt["truth_scope"] != \
            "injected_equal_length_structure_only":
        raise ValueError("unexpected challenge case manifest")
    _assert_hash(SOURCE / "colcen_native_contexts.fa", receipt["source_context_sha256"])
    _assert_hash(SOURCE / "colcen_native_spanners.fastq.gz",
                 receipt["original_native_fastq_sha256"])
    pilot_receipt = json.loads((PILOT / "receipt.json").read_text())
    _assert_hash(PILOT / "read_trims.json", pilot_receipt["read_trims_sha256"])
    _assert_hash(PILOT / "per_case.jsonl", pilot_receipt["per_case_sha256"])
    _assert_hash(B1, pilot_receipt["frozen_b1_inputs_sha256"])
    if pilot_receipt["eligible_read_count"] != 3 or \
            pilot_receipt["pairing_status"] != protocol["pairing_status"]:
        raise ValueError("prior native trim pilot eligibility changed")
    with B1.open() as handle:
        monomers = json.loads(handle.readline())["candidate_monomers"]
    if set(monomers) != {"C1", "C2", "C3"}:
        raise ValueError("supplied monomer catalogue changed")
    return protocol, receipt["cases"], monomers, pilot_receipt


def _reconstruct_native_arrays(original: str, pilot_receipt: dict[str, Any]) -> \
        tuple[dict[str, str], list[dict[str, Any]]]:
    start, end = 5000, 8560
    left = original[start-FLANK_BP:start]
    right = original[end:end+FLANK_BP]
    raw_reads = fastq_records(SOURCE / "colcen_native_spanners.fastq.gz")
    trims = json.loads((PILOT / "read_trims.json").read_text())
    source_receipt = json.loads((SOURCE / "source_and_alignment_receipt.json").read_text())
    _assert_hash(SOURCE / "source_and_alignment_receipt.json",
                 pilot_receipt["source_receipt_sha256"])
    paf_path = SOURCE / "colcen_native_spanners.paf"
    _assert_hash(paf_path, pilot_receipt["native_PAF_sha256"])
    alignments = {row["read_id"]: row for row in source_receipt["native_molecules"]
                  if row["context"] == "C3"}
    if len(trims) != 3 or {r["read_id"] for r in trims} != \
            set(pilot_receipt["read_ids"]) or set(alignments) != set(pilot_receipt["read_ids"]):
        raise ValueError("native C3 read ID set changed")
    arrays: dict[str, str] = {}
    for archived in trims:
        read_id = archived["read_id"]
        read = raw_reads[read_id]
        if hashlib.sha256(read.encode()).hexdigest() != archived["original_read_sha256"]:
            raise ValueError("original read sequence hash changed")
        if archived["status"] != "ELIGIBLE" or archived["strand"] != alignments[read_id]["strand"]:
            raise ValueError("archived native read ineligible or strand changed")
        oriented = read if archived["strand"] == "+" else reverse_complement(read)
        trim = trim_between_flanks(oriented, left, right)
        projected = primary_paf_array_interval(paf_path, read_id, "C3", start, end)
        if trim["status"] != "ELIGIBLE" or trim["trim_interval_0based"] != \
                archived["trim_interval_0based"] or projected != \
                archived["primary_PAF_projected_interval_0based"] or \
                projected != trim["trim_interval_0based"]:
            raise ValueError("native read natural-flank/PAF trim changed")
        arrays[read_id] = str(trim["array_sequence"])
    return arrays, trims


def classify_paths(assembly: dict[str, Any], reads: dict[str, dict[str, Any]]) -> \
        dict[str, str]:
    """Technical path decision only; never call M2 biological audit as verified."""
    if assembly["state"] != "RESOLVED":
        return {"state": "ABSTAIN", "reason": "assembly_" + str(assembly["reason"])}
    if len(reads) < 3 or any(row["state"] != "RESOLVED" for row in reads.values()):
        return {"state": "ABSTAIN", "reason": "insufficient_resolved_native_read_paths"}
    read_paths = {tuple(row["labels"]) for row in reads.values()}
    if len(read_paths) != 1:
        return {"state": "ABSTAIN", "reason": "mixed_native_read_paths"}
    return {"state": "SUPPORTED" if tuple(assembly["labels"]) in read_paths
            else "DISCORDANT", "reason": "technical_label_orientation_path_comparison_only"}


def classify_length(assembly_bp: int, read_lengths: list[int]) -> dict[str, Any]:
    if len(read_lengths) < 3:
        return {"state": "ABSTAIN", "reason": "fewer_than_three_eligible_reads"}
    median_bp = statistics.median(read_lengths)
    threshold = max(2, math.ceil(max(assembly_bp, median_bp) * 0.05))
    delta = assembly_bp - median_bp
    return {"state": "DISCORDANT" if abs(delta) > threshold else "SUPPORTED",
            "reason": "frozen_five_percent_span_threshold", "median_read_span_bp": median_bp,
            "signed_bp_delta_assembly_minus_read": delta, "threshold_bp": threshold}


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if len(rows) != 6:
        raise ValueError("expected six full-denominator challenge cases")
    positive = [row for row in rows if row["truth_status"] == "engineered_structural_positive"]
    negative = [row for row in rows if row["truth_status"] == "intact_negative"]
    if len(positive) != 5 or len(negative) != 1:
        raise ValueError("unexpected positive/negative case denominator")
    out: dict[str, Any] = {"case_count": 6, "positive_count": 5, "negative_count": 1,
                           "path_identifiable_positive_count": sum(
                               row["path_identifiability"] == "IDENTIFIABLE" for row in positive),
                           "biological_accuracy": "not_evaluated"}
    for route in ("m2_technical", "length_baseline"):
        out[route] = {
            "positive_detected": sum(row[route]["state"] == "DISCORDANT" for row in positive),
            "positive_missed_supported": sum(row[route]["state"] == "SUPPORTED" for row in positive),
            "positive_abstained": sum(row[route]["state"] == "ABSTAIN" for row in positive),
            "negative_supported": sum(row[route]["state"] == "SUPPORTED" for row in negative),
            "negative_false_positive": sum(row[route]["state"] == "DISCORDANT" for row in negative),
            "negative_abstained": sum(row[route]["state"] == "ABSTAIN" for row in negative),
            "identifiable_positive_detected": sum(row["path_identifiability"] ==
                                                  "IDENTIFIABLE" and row[route]["state"] ==
                                                  "DISCORDANT" for row in positive),
        }
    return out


def run(output: Path) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError("equal-length score output exists; preserve prior evidence")
    protocol, cases, monomers, pilot_receipt = _read_fixed_inputs()
    contexts = fasta_records(SOURCE / "colcen_native_contexts.fa")
    original = next(sequence for name, sequence in contexts.items() if name.startswith("C3|"))
    original_array = original[5000:8560]
    if hashlib.sha256(original_array.encode()).hexdigest() != cases[0]["source_array_sha256"]:
        raise ValueError("source C3 array hash mismatch")
    native_arrays, trims = _reconstruct_native_arrays(original, pilot_receipt)
    read_paths = {read_id: path_summary(sequence, monomers)
                  for read_id, sequence in sorted(native_arrays.items())}
    source_path = path_summary(original_array, monomers)
    rows: list[dict[str, Any]] = []
    for case in cases:
        case_id = case["case_id"]
        case_path = CHALLENGE / (case_id + ".fa")
        _assert_hash(case_path, case["edited_fasta_sha256"])
        edited = next(iter(fasta_records(case_path).values()))
        if edited[:5000] != original[:5000] or edited[8560:] != original[8560:]:
            raise ValueError("edited natural flanks changed")
        if case["edited_array_interval"] != [5000, 8560] or len(edited[5000:8560]) != 3560:
            raise ValueError("equal-length case interval changed")
        _assert_array = hashlib.sha256(edited[5000:8560].encode()).hexdigest()
        if _assert_array != case["edited_array_sha256"] or \
                case["source_array_sha256"] != hashlib.sha256(original_array.encode()).hexdigest():
            raise ValueError("controlled array hash mismatch")
        assembly_trim = trim_between_flanks(edited, original[5000-FLANK_BP:5000],
                                            original[8560:8560+FLANK_BP])
        if assembly_trim["status"] != "ELIGIBLE" or \
                assembly_trim["trim_interval_0based"] != [5000, 8560]:
            raise ValueError("edited assembly natural-flank boundaries ambiguous")
        assembly_path = path_summary(str(assembly_trim["array_sequence"]), monomers)
        if source_path["state"] != "RESOLVED" or assembly_path["state"] != "RESOLVED":
            identifiability = "UNRESOLVED_PATH"
        elif source_path["labels"] == assembly_path["labels"]:
            identifiability = "NOT_IDENTIFIABLE_BY_OPERATIONAL_LABEL_ORIENTATION_PATH"
        else:
            identifiability = "IDENTIFIABLE"
        rows.append({"case_id": case_id, "operation": case["operation"],
                     "truth_status": case["truth_status"], "injected_bp_delta": case["injected_bp_delta"],
                     "path_identifiability": identifiability, "assembly_path": assembly_path,
                     "native_read_paths": read_paths,
                     "m2_technical": classify_paths(assembly_path, read_paths),
                     "length_baseline": classify_length(3560,
                                                        [len(sequence) for sequence in native_arrays.values()]),
                     "biological_audit_state": "NOT_EVALUATED_PAIRING_UNVERIFIED"})
    summary = summarize(rows)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "per_case.jsonl").open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    log = ["validated frozen protocol, generated receipt, six case FASTAs and array hashes",
           "validated source FASTQ, context, PAF, prior trim receipt and B1 monomer catalogue",
           "reconstructed all three original C3 read trims and matched archived plus PAF coordinates",
           "ran unchanged M2 decompose defaults and frozen 5% span threshold",
           "biological audit not evaluated; same-study donor/haplotype pairing unverified"]
    (output / "run_log.txt").write_text("\n".join(log) + "\n")
    receipt = {"protocol_sha256": sha256(CHALLENGE / "protocol.json"),
               "challenge_receipt_sha256": sha256(CHALLENGE / "receipt.json"),
               "native_fastq_sha256": sha256(SOURCE / "colcen_native_spanners.fastq.gz"),
               "native_context_sha256": sha256(SOURCE / "colcen_native_contexts.fa"),
               "native_PAF_sha256": sha256(SOURCE / "colcen_native_spanners.paf"),
               "prior_trim_receipt_sha256": sha256(PILOT / "receipt.json"),
               "prior_trim_table_sha256": sha256(PILOT / "read_trims.json"),
               "M2_prototype_sha256": sha256(PROTOTYPE),
               "score_script_sha256": sha256(Path(__file__)),
               "read_ids": sorted(native_arrays), "read_trim_intervals_0based": {
                   r["read_id"]: r["trim_interval_0based"] for r in trims},
               "case_count": 6, "per_case_sha256": sha256(output / "per_case.jsonl"),
               "summary_sha256": sha256(output / "summary.json"),
               "run_log_sha256": sha256(output / "run_log.txt"),
               "biological_accuracy": "not_evaluated", "split": protocol["split"]}
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
