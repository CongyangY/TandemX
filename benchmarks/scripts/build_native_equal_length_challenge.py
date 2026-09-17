"""Freeze equal-length structural edits on a native assembly context.

The source reads remain unchanged. These cases test whether a path method
detects structure beyond the trivial read--assembly span length difference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.build_native_read_collapse import read_named_fasta, sha256_file


def digest_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def transform(tiles: list[str], operation: str) -> list[str]:
    output = tiles.copy()
    if operation == "intact":
        pass
    elif operation == "invert_tile_05":
        output[5] = reverse_complement(output[5])
    elif operation == "invert_block_05_09":
        output[5:9] = [reverse_complement(tile) for tile in reversed(output[5:9])]
    elif operation == "swap_adjacent_05_06":
        output[5], output[6] = output[6], output[5]
    elif operation == "swap_distant_05_15":
        output[5], output[15] = output[15], output[5]
    elif operation == "replace_tile_05_with_06":
        output[5] = output[6]
    else:
        raise ValueError("Unknown structural operation")
    return output


OPERATIONS = ("intact", "invert_tile_05", "invert_block_05_09",
              "swap_adjacent_05_06", "swap_distant_05_15",
              "replace_tile_05_with_06")


def generate(config_path: Path, outdir: Path) -> dict:
    config = json.loads(config_path.read_text())
    if config.get("split") != "development" or config.get("read_pairing_status") != \
            "lineage_context_supported_donor_unverified":
        raise ValueError("Native structural input is development only")
    reads = Path(config["raw_read_bundle"])
    if sha256_file(reads) != config["raw_read_bundle_sha256"]:
        raise ValueError("Original read-bundle hash mismatch")
    context = next((row for row in config["contexts"] if row["array_id"] == "C3"), None)
    if context is None:
        raise ValueError("C3 source context absent")
    source_path = Path(context["reference_fasta"])
    if sha256_file(source_path) != context["reference_fasta_sha256"]:
        raise ValueError("Source context hash mismatch")
    source = read_named_fasta(source_path, context["reference_record_id"])
    start, end, unit = context["array_start0"], context["array_end0"], context["operational_unit_bp"]
    array = source[start:end]
    if (unit != 178 or len(array) != 20 * unit or digest_sequence(array) != context["array_sequence_sha256"]):
        raise ValueError("Frozen C3 20-tile source changed")
    tiles = [array[index:index + unit] for index in range(0, len(array), unit)]
    outdir.mkdir(parents=True, exist_ok=False)
    rows = []
    for operation in OPERATIONS:
        edited_tiles = transform(tiles, operation)
        edited_array = "".join(edited_tiles)
        if len(edited_array) != len(array) or (operation != "intact" and edited_array == array):
            raise ValueError("Structural edit is silent or changes length")
        edited = source[:start] + edited_array + source[end:]
        path = outdir / f"C3_{operation}.fa"
        path.write_text(f">C3_{operation}\n{edited}\n", encoding="ascii")
        rows.append({"case_id": f"C3_{operation}", "operation": operation,
                     "source_array_interval": [start, end],
                     "edited_array_interval": [start, end],
                     "injected_bp_delta": 0,
                     "source_array_sha256": digest_sequence(array),
                     "edited_array_sha256": digest_sequence(edited_array),
                     "edited_fasta_sha256": sha256_file(path),
                     "truth_status": "intact_negative" if operation == "intact" else "engineered_structural_positive",
                     "truth_scope": "injected_equal_length_structure_only"})
    if sha256_file(reads) != config["raw_read_bundle_sha256"]:
        raise ValueError("Original read bundle changed during generation")
    receipt = {"schema_version": 1, "split": "development",
               "source_config_sha256": sha256_file(config_path),
               "generator_sha256": sha256_file(Path(__file__)),
               "source_context_sha256": sha256_file(source_path),
               "original_native_fastq_sha256": config["raw_read_bundle_sha256"],
               "source_assembly_lineage_count": 1,
               "read_pairing_status": config["read_pairing_status"],
               "truth_scope": "injected_equal_length_structure_only", "cases": rows}
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    generate(args.config, args.outdir)


if __name__ == "__main__":
    main()
