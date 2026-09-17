"""Independently reconstruct every archived native-context deletion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.scripts.build_native_read_collapse import read_named_fasta, sha256_file


def verify(config_path: Path, generated: Path) -> dict:
    config = json.loads(config_path.read_text())
    receipt = json.loads((generated / "receipt.json").read_text())
    if (receipt["config_sha256"] != sha256_file(config_path)
            or receipt["raw_read_bundle_sha256"] != config["raw_read_bundle_sha256"]
            or sha256_file(Path(config["raw_read_bundle"])) != config["raw_read_bundle_sha256"]):
        raise ValueError("Frozen source/read receipt mismatch")
    source_rows = {row["array_id"]: row for row in config["contexts"]}
    actual_ids = set()
    for case in receipt["cases"]:
        case_id = case["case_id"]
        if case_id in actual_ids:
            raise ValueError("Duplicate edited case")
        actual_ids.add(case_id)
        row = source_rows[case["array_id"]]
        source_path = Path(row["reference_fasta"])
        if sha256_file(source_path) != row["reference_fasta_sha256"]:
            raise ValueError("Source context hash mismatch")
        original = read_named_fasta(source_path, row["reference_record_id"])
        left, right = case["source_deleted_interval"]
        expected = original[:left] + original[right:]
        if (right - left != case["injected_deleted_bp"]
                or case["source_array_interval"] != [row["array_start0"], row["array_end0"]]
                or case["edited_array_bp"] != row["array_end0"] - row["array_start0"] - (right-left)):
            raise ValueError(f"Edit ledger mismatch: {case_id}")
        fasta_path = generated / f"{case_id}.fa"
        if sha256_file(fasta_path) != case["edited_fasta_sha256"]:
            raise ValueError(f"Edited FASTA hash mismatch: {case_id}")
        lines = fasta_path.read_text(encoding="ascii").splitlines()
        if len(lines) != 2 or lines[0] != f">{case_id}" or lines[1] != expected:
            raise ValueError(f"Edited sequence mismatch: {case_id}")
    expected_ids = {f"{array_id}_{suffix}" for array_id in source_rows for suffix in (
        "terminal_000", "terminal_025", "terminal_050", "terminal_075", "terminal_100",
        "internal_025", "internal_050", "internal_075", "boundary_left_025")}
    if actual_ids != expected_ids:
        raise ValueError("Case denominator mismatch")
    return {"status": "verified_injected_edit_only", "case_count": len(actual_ids),
            "context_count": len(source_rows), "raw_read_bundle_sha256": config["raw_read_bundle_sha256"],
            "receipt_sha256": sha256_file(generated / "receipt.json")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--generated", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.config, args.generated)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
