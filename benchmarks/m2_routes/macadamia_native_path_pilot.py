"""Frozen M2 path feasibility on the provisional Macadamia exact-bp edits."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmarks.m2_routes.native_read_pilot import (
    fasta_records, fastq_records, path_summary, reverse_complement,
)
from benchmarks.scripts.build_native_read_collapse import sha256_file


ROOT = Path(__file__).resolve().parents[2]
EDIT = ROOT / "benchmarks/controlled_collapse/macadamia_bp_provisional_v1"
SOURCE = ROOT / "benchmarks/controlled_collapse/macadamia_assembly_candidate_v1_20260918"
DISCOVERY_MONOMERS = Path(
    "/Volumes/T7/Codex/TandemX/results/"
    "macadamia_jansenii_donor_matched_collapse_v1_20260909/"
    "run/discover/monomers.fa"
)


def technical_decision(assembly: dict, reads: dict[str, dict], minimum_reads: int = 3) -> dict:
    if len(reads) < minimum_reads:
        return {"state": "ABSTAIN", "reason": "fewer_than_three_eligible_records"}
    if assembly["state"] != "RESOLVED":
        return {"state": "ABSTAIN", "reason": "assembly_" + str(assembly["reason"])}
    if any(value["state"] != "RESOLVED" for value in reads.values()):
        return {"state": "ABSTAIN", "reason": "one_or_more_original_read_paths_unresolved"}
    paths = {tuple(value["labels"]) for value in reads.values()}
    if len(paths) != 1:
        return {"state": "ABSTAIN", "reason": "mixed_original_read_paths"}
    return {"state": "SUPPORTED" if next(iter(paths)) == tuple(assembly["labels"])
            else "DISCORDANT", "reason": "technical_label_orientation_path_comparison_only"}


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Preserve prior M2 pilot output")
    protocol_path = EDIT / "m2_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if (protocol["status"] != "frozen_before_M2_scoring" or
            protocol["case_denominator"] != 9 or
            protocol["selected_original_record_denominator"] != 7 or
            protocol["template_length_bp"] != 144):
        raise ValueError("M2 pilot does not match frozen protocol")
    paths = {
        "source_context_fasta_sha256": SOURCE / "selected_context.fa",
        "original_selected_fastq_sha256": EDIT / "selected_original_reads.fastq.gz",
        "span_score_protocol_sha256": EDIT / "score_protocol.json",
        "edit_receipt_sha256": EDIT / "generated/receipt.json",
        "existing_discovery_monomers_fasta_sha256": DISCOVERY_MONOMERS,
        "single_family_monomer_fasta_sha256": EDIT / "monomer_template.fa",
        "frozen_m2_prototype_sha256": ROOT / "benchmarks/m2_routes/alignment/prototype.py",
    }
    for key, path in paths.items():
        if sha256_file(path) != protocol[key]:
            raise ValueError(f"Frozen M2 input changed: {path}")
    template = fasta_records(EDIT / "monomer_template.fa")
    if len(template) != 1 or len(next(iter(template.values()))) != 144:
        raise ValueError("Expected one 144-bp monomer template")
    template_sequence = next(iter(template.values()))
    source_monomers = fasta_records(DISCOVERY_MONOMERS)
    match = [seq for header, seq in source_monomers.items()
             if "family_id=TXF000219;monomer_id=TXM000219;" in header]
    if match != [template_sequence]:
        raise ValueError("Supplied monomer does not match existing discovery")
    monomers = {"TXM000219": template_sequence}
    read_trims_path = EDIT / "span_score/read_trims.json"
    span_summary_path = EDIT / "span_score/summary.json"
    trim_rows = json.loads(read_trims_path.read_text())
    span_summary = json.loads(span_summary_path.read_text())
    reads = fastq_records(EDIT / "selected_original_reads.fastq.gz")
    if len(trim_rows) != 7 or set(reads) != {row["read_id"] for row in trim_rows}:
        raise ValueError("Frozen read denominator changed")
    if span_summary["read_trims_sha256"] != sha256_file(read_trims_path) or \
            span_summary["case_count"] != 9 or span_summary["eligible_record_count"] != 7:
        raise ValueError("Frozen span scorer evidence changed")
    read_paths = {}
    for row in trim_rows:
        rid = row["read_id"]
        sequence = reads[rid]
        if row["status"] != "ELIGIBLE" or \
                hashlib.sha256(sequence.encode()).hexdigest() != row["original_read_sequence_sha256"]:
            raise ValueError("Original read eligibility or hash changed")
        oriented = sequence if row["strand"] == "+" else reverse_complement(sequence)
        start, end = row["trim_interval_0based"]
        if end - start != row["observed_array_span_bp"]:
            raise ValueError("Original read span changed")
        read_paths[rid] = path_summary(oriented[start:end], monomers)
    receipt = json.loads((EDIT / "generated/receipt.json").read_text())
    score_rows = [json.loads(line) for line in (EDIT / "span_score/per_case.jsonl").read_text().splitlines()]
    if len(receipt["cases"]) != 9 or [row["case_id"] for row in score_rows] != \
            [case["case_id"] for case in receipt["cases"]] or \
            span_summary["per_case_sha256"] != sha256_file(EDIT / "span_score/per_case.jsonl"):
        raise ValueError("Frozen edited-case denominator changed")
    rows = []
    for case, score in zip(receipt["cases"], score_rows):
        fasta = EDIT / "generated" / f"{case['case_id']}.fa"
        if sha256_file(fasta) != case["edited_fasta_sha256"]:
            raise ValueError("Edited FASTA changed")
        sequences = fasta_records(fasta)
        if set(sequences) != {case["case_id"]}:
            raise ValueError("Edited FASTA case ID changed")
        interval = case["edited_array_context_interval"]
        if interval[1] - interval[0] != score["direct_assembly_span_bp"]:
            raise ValueError("M2 assembly trim conflicts with direct span scorer")
        assembly = path_summary(sequences[case["case_id"]][interval[0]:interval[1]], monomers)
        rows.append({"case_id": case["case_id"], "injected_deleted_bp": case["injected_deleted_bp"],
                     "assembly_path": assembly, "m2_technical": technical_decision(assembly, read_paths),
                     "biological_audit_state": "NOT_EVALUATED_PAIRING_AND_ZMW_UNVERIFIED"})
    output.mkdir(parents=True, exist_ok=False)
    (output / "read_paths.json").write_text(json.dumps(read_paths, indent=2, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    result = {
        "case_count": len(rows), "eligible_original_record_count": len(read_paths),
        "original_read_resolved_count": sum(x["state"] == "RESOLVED" for x in read_paths.values()),
        "assembly_resolved_count": sum(x["assembly_path"]["state"] == "RESOLVED" for x in rows),
        "technical_discordant_count": sum(x["m2_technical"]["state"] == "DISCORDANT" for x in rows),
        "technical_supported_count": sum(x["m2_technical"]["state"] == "SUPPORTED" for x in rows),
        "technical_abstain_count": sum(x["m2_technical"]["state"] == "ABSTAIN" for x in rows),
        "biological_accuracy": "not_evaluated", "molecule_identity": "unverified_ZMW",
        "m2_protocol_sha256": sha256_file(protocol_path),
        "read_paths_sha256": sha256_file(output / "read_paths.json"),
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
