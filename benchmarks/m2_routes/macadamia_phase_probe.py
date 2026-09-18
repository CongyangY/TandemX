"""Post-result phase probe for Macadamia; never a benchmark prediction."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import edlib

from benchmarks.m2_routes.alignment.prototype import decompose
from benchmarks.m2_routes.native_read_pilot import fasta_records, fastq_records, reverse_complement
from benchmarks.scripts.build_native_read_collapse import sha256_file


ROOT = Path(__file__).resolve().parents[2]
EDIT = ROOT / "benchmarks/controlled_collapse/macadamia_bp_provisional_v1"


def select_phase(reference_array: str, motif: str) -> tuple[str, dict]:
    if len(motif) != 144 or len(reference_array) < 144:
        raise ValueError("Phase probe requires 144-bp motif and source array")
    candidates = []
    for strand, base in (("+", motif), ("-", reverse_complement(motif))):
        for shift in range(144):
            rotated = base[shift:] + base[:shift]
            distance = edlib.align(rotated, reference_array[:144], mode="NW",
                                   task="distance")["editDistance"]
            candidates.append((distance, strand, shift, rotated))
    best = min(candidates)
    ties = sum(row[0] == best[0] for row in candidates)
    return best[3], {"first_144bp_edit_distance": best[0], "strand": best[1],
                     "rotation_bp": best[2], "best_phase_tie_count": ties}


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Preserve exploratory phase probe")
    context = ROOT / "benchmarks/controlled_collapse/macadamia_assembly_candidate_v1_20260918/selected_context.fa"
    monomer = EDIT / "monomer_template.fa"
    reads_path = EDIT / "selected_original_reads.fastq.gz"
    protocol = json.loads((EDIT / "m2_protocol.json").read_text())
    if sha256_file(context) != protocol["source_context_fasta_sha256"] or \
            sha256_file(monomer) != protocol["single_family_monomer_fasta_sha256"] or \
            sha256_file(reads_path) != protocol["original_selected_fastq_sha256"]:
        raise ValueError("Frozen source changed")
    source = next(iter(fasta_records(context).values()))
    motif = next(iter(fasta_records(monomer).values()))
    phase, selection = select_phase(source[3000:6155], motif)
    reads = fastq_records(reads_path)
    trims = json.loads((EDIT / "span_score/read_trims.json").read_text())
    cases = json.loads((EDIT / "generated/receipt.json").read_text())["cases"]
    rows = []
    for row in trims:
        sequence = reads[row["read_id"]]
        oriented = sequence if row["strand"] == "+" else reverse_complement(sequence)
        start, end = row["trim_interval_0based"]
        result = decompose(oriented[start:end], {"TXM000219": phase})
        rows.append({"kind": "original_record", "id": row["read_id"],
                     "state": result.state, "reason": result.reason,
                     "copy_count": len(result.copies),
                     "labels": [copy.label + copy.orientation for copy in result.copies]})
    for case in cases:
        fasta = EDIT / "generated" / (case["case_id"] + ".fa")
        if sha256_file(fasta) != case["edited_fasta_sha256"]:
            raise ValueError("Edited FASTA changed")
        sequence = fasta_records(fasta)[case["case_id"]]
        start, end = case["edited_array_context_interval"]
        result = decompose(sequence[start:end], {"TXM000219": phase})
        rows.append({"kind": "edited_assembly", "id": case["case_id"],
                     "state": result.state, "reason": result.reason,
                     "copy_count": len(result.copies),
                     "labels": [copy.label + copy.orientation for copy in result.copies]})
    output.mkdir(parents=True, exist_ok=False)
    with (output / "per_input.jsonl").open("w") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    summary = {"status": "post_result_exploratory_phase_probe_not_frozen_M2_score",
               "selected_phase_from": "intact_assembly_first_144bp",
               "phase": selection,
               "phase_template_sha256": hashlib.sha256(phase.encode()).hexdigest(),
               "original_record_resolved": sum(row["state"] == "RESOLVED" for row in rows
                                               if row["kind"] == "original_record"),
               "original_record_count": len(trims),
               "edited_assembly_resolved": sum(row["state"] == "RESOLVED" for row in rows
                                               if row["kind"] == "edited_assembly"),
               "edited_assembly_count": len(cases),
               "per_input_sha256": sha256_file(output / "per_input.jsonl"),
               "cannot_promote_frozen_M2_accuracy": True}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
