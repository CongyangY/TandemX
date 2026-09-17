"""Independently reconstruct the Ey15 native-context exact-bp edits."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.build_native_read_collapse import read_named_fasta, sha256_file


EXPECTED_LABELS = {
    "terminal_000", "terminal_025", "terminal_050", "terminal_075", "terminal_100",
    "internal_025", "internal_050", "internal_075", "boundary_left_025",
}


def verify(manifest_path: Path, context_path: Path, reads_path: Path, generated: Path) -> dict:
    manifest = json.loads(manifest_path.read_text())
    receipt = json.loads((generated / "receipt.json").read_text())
    if receipt["source_manifest_sha256"] != sha256_file(manifest_path):
        raise ValueError("Source manifest changed")
    if (receipt["source_context_sha256"] != sha256_file(context_path)
            or receipt["original_native_fastq_sha256"] != sha256_file(reads_path)
            or manifest["artifact_sha256"][context_path.name] != sha256_file(context_path)
            or manifest["artifact_sha256"][reads_path.name] != sha256_file(reads_path)):
        raise ValueError("Native source/read hash mismatch")
    array = manifest["array"]
    name = (f"E1|{array['chromosome']}:{array['context_start0']}-{array['context_end0']}|"
            f"array:{array['start0']}-{array['end0']}|family:{array['family_id']}")
    source = read_named_fasta(context_path, name)
    start, end = array["context_array_start0"], array["context_array_end0"]
    if end - start != array["length_bp"] or len(source) != array["context_end0"] - array["context_start0"]:
        raise ValueError("Source array coordinate mismatch")
    observed = set()
    for case in receipt["cases"]:
        case_id = case["case_id"]
        if case_id in observed:
            raise ValueError("Duplicate case")
        observed.add(case_id)
        left, right = case["source_deleted_context_interval"]
        removed = right - left
        if not (start <= left <= right <= end):
            raise ValueError(f"Deletion outside array: {case_id}")
        expected = source[:left] + source[right:]
        if (case["source_array_context_interval"] != [start, end]
                or case["source_array_genome_interval"] != [array["start0"], array["end0"]]
                or case["source_deleted_genome_interval"] !=
                    [array["context_start0"] + left, array["context_start0"] + right]
                or case["edited_array_context_interval"] != [start, end - removed]
                or case["edited_array_bp"] != end - start - removed
                or case["injected_deleted_bp"] != removed
                or case["edited_breakpoint0"] != left
                or case["source_array_sha256"] != hashlib.sha256(source[start:end].encode()).hexdigest()
                or case["edited_array_sha256"] != hashlib.sha256(expected[start:end-removed].encode()).hexdigest()):
            raise ValueError(f"Edit ledger mismatch: {case_id}")
        path = generated / f"{case_id}.fa"
        if sha256_file(path) != case["edited_fasta_sha256"]:
            raise ValueError(f"Edited FASTA hash mismatch: {case_id}")
        if path.read_text(encoding="ascii") != f">{case_id}\n{expected}\n":
            raise ValueError(f"Edited FASTA sequence mismatch: {case_id}")
    if observed != {"E1_" + label for label in EXPECTED_LABELS}:
        raise ValueError("Case denominator mismatch")
    present = {p.name for p in generated.iterdir()}
    if present - {"verification.json", "score_protocol.json"} != {"receipt.json"} | {name + ".fa" for name in observed}:
        raise ValueError("Unexpected or missing generated file")
    return {"status": "verified_injected_bp_only", "case_count": len(observed),
            "original_native_fastq_sha256": sha256_file(reads_path),
            "receipt_sha256": sha256_file(generated / "receipt.json")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--generated", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(verify(args.manifest, args.context, args.reads, args.generated),
                                   indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
