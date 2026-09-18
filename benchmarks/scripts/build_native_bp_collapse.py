"""Edit a nominated native assembly interval by exact base-pair amounts.

The interval need not have verified monomer boundaries. Original reads stay
unchanged, and the ledger only states the injected assembly difference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.build_native_read_collapse import read_named_fasta, sha256_file


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def edit_specs(array_bp: int) -> list[tuple[str, int, int]]:
    if array_bp < 100:
        raise ValueError("Native array must be at least 100 bp")
    specs = []
    for percent in (0, 25, 50, 75, 100):
        removed = array_bp * percent // 100
        specs.append((f"terminal_{percent:03d}", array_bp - removed, array_bp))
    for percent in (25, 50, 75):
        removed = array_bp * percent // 100
        first = (array_bp - removed) // 2
        specs.append((f"internal_{percent:03d}", first, first + removed))
    specs.append(("boundary_left_025", 0, array_bp // 4))
    return specs


def generate(manifest_path: Path, context_path: Path, reads_path: Path, outdir: Path) -> dict:
    manifest = json.loads(manifest_path.read_text())
    status = manifest.get("status")
    if manifest.get("split") != "development" or status not in {
            "development_source_enrolled_native_interval_supported",
            "development_provisional_record_support"}:
        raise ValueError("Expected supported development manifest")
    if (status == "development_provisional_record_support" and
            manifest.get("support_tier") != "full_reference_records_zmw_unverified"):
        raise ValueError("Provisional support tier must identify unverified ZMWs")
    array = manifest["array"]
    if sha256_file(context_path) != manifest["artifact_sha256"][context_path.name]:
        raise ValueError("Context FASTA hash mismatch")
    if sha256_file(reads_path) != manifest["artifact_sha256"][reads_path.name]:
        raise ValueError("Original native FASTQ hash mismatch")
    name = manifest.get("context_record_id") or (
        f"E1|{array['chromosome']}:{array['context_start0']}-{array['context_end0']}|"
        f"array:{array['start0']}-{array['end0']}|family:{array['family_id']}")
    case_prefix = manifest.get("case_prefix", "E1")
    if not case_prefix or not case_prefix.replace("_", "").isalnum():
        raise ValueError("Invalid case prefix")
    source = read_named_fasta(context_path, name)
    start, end = array["context_array_start0"], array["context_array_end0"]
    if (start < 1000 or end > len(source) - 1000 or end <= start
            or end - start != array["length_bp"]
            or len(source) != array["context_end0"] - array["context_start0"]):
        raise ValueError("Array/context coordinate mismatch")
    source_array = source[start:end]
    outdir.mkdir(parents=True, exist_ok=False)
    cases = []
    for label, local_start, local_end in edit_specs(len(source_array)):
        delete_start, delete_end = start + local_start, start + local_end
        edited = source[:delete_start] + source[delete_end:]
        removed = delete_end - delete_start
        case_id = f"{case_prefix}_{label}"
        path = outdir / f"{case_id}.fa"
        path.write_text(f">{case_id}\n{edited}\n", encoding="ascii")
        if len(edited) != len(source) - removed:
            raise AssertionError("Edited bp conservation failed")
        cases.append({"case_id": case_id,
                      "edit_type": label.rsplit("_", 1)[0],
                      "requested_deletion_percent": int(label.rsplit("_", 1)[1]),
                      "source_context_genome_interval": [array["context_start0"], array["context_end0"]],
                      "source_array_genome_interval": [array["start0"], array["end0"]],
                      "source_array_context_interval": [start, end],
                      "source_deleted_context_interval": [delete_start, delete_end],
                      "source_deleted_genome_interval": [array["context_start0"] + delete_start,
                                                        array["context_start0"] + delete_end],
                      "edited_array_context_interval": [start, end - removed],
                      "edited_breakpoint0": delete_start,
                      "source_array_bp": len(source_array),
                      "edited_array_bp": len(source_array) - removed,
                      "injected_deleted_bp": removed,
                      "source_array_sha256": sha256_sequence(source_array),
                      "edited_array_sha256": sha256_sequence(edited[start:end - removed]),
                      "edited_fasta_sha256": sha256_file(path),
                      "monomer_copy_truth": "unknown_noninteger_period_and_unverified_boundaries",
                      "read_pairing_status": manifest["source_pairing"],
                      "truth_scope": "injected_bp_delta_only"})
    if sha256_file(reads_path) != manifest["artifact_sha256"][reads_path.name]:
        raise ValueError("Original native FASTQ changed during generation")
    receipt = {"schema_version": 1, "split": "development", "material_id": manifest["sample"],
               "source_manifest_sha256": sha256_file(manifest_path),
               "source_context_sha256": sha256_file(context_path),
               "original_native_fastq_sha256": sha256_file(reads_path),
               "generator_sha256": sha256_file(Path(__file__)),
               "raw_reads_modified": False, "truth_scope": "injected_bp_delta_only",
               "independent_donor_count": 1, "cases": cases}
    if status == "development_provisional_record_support":
        receipt["source_support_tier"] = manifest["support_tier"]
        receipt["independent_molecule_count"] = None
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    generate(args.manifest, args.context, args.reads, args.outdir)


if __name__ == "__main__":
    main()
