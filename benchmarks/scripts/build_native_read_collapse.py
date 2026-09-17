"""Generate exact local assembly edits while retaining an unchanged native-read bundle.

The injected delta is known. The reference assembly and reads are separate
observations; this command does not declare biological copy-number truth.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_one_fasta(path: Path) -> tuple[str, str]:
    if path.stat().st_size > 2_000_000:
        raise ValueError("Native context exceeds the bounded 2-MB input cap")
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Expected exactly one FASTA record")
    name = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).upper()
    if not name or not sequence or set(sequence) - set("ACGT"):
        raise ValueError("Native context must have one named ACGT sequence")
    return name, sequence


def fastq_read_ids(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    identifiers: set[str] = set()
    with opener(path, "rt", encoding="ascii") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            identifier = header[1:].strip().split()[0] if header.startswith("@") else ""
            if (not identifier or identifier in identifiers or not sequence
                    or set(sequence.upper()) - set("ACGT") or not separator.startswith("+")
                    or len(quality) != len(sequence)):
                raise ValueError("Invalid, duplicate, or truncated native FASTQ record")
            identifiers.add(identifier)
    if not identifiers:
        raise ValueError("Native read bundle is empty")
    return identifiers


def edit_specs(unit_count: int) -> list[tuple[str, int, int]]:
    if unit_count < 4 or unit_count % 4:
        raise ValueError("Array requires a positive multiple of four operational units")
    specs = []
    for percent in (0, 25, 50, 75, 100):
        deleted = unit_count * percent // 100
        specs.append((f"terminal_{percent:03d}", unit_count - deleted, unit_count))
    for percent in (25, 50, 75):
        deleted = unit_count * percent // 100
        start = (unit_count - deleted) // 2
        specs.append((f"internal_{percent:03d}", start, start + deleted))
    specs.append(("boundary_left_025", 0, unit_count // 4))
    return specs


def generate(config_path: Path, outdir: Path) -> dict:
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1 or config.get("split") != "development":
        raise ValueError("Expected a declared development source config")
    rows = config.get("contexts")
    if not isinstance(rows, list) or not rows:
        raise ValueError("No contexts declared")
    read_bundle = Path(config["raw_read_bundle"])
    if sha256_file(read_bundle) != config["raw_read_bundle_sha256"]:
        raise ValueError("Native raw-read bundle hash mismatch")
    read_ids = fastq_read_ids(read_bundle)
    pairing = config["read_pairing_status"]
    if pairing not in {"lineage_context_supported_donor_unverified", "verified_same_donor_haplotype"}:
        raise ValueError("Unrecognized read-pairing status")
    if pairing == "verified_same_donor_haplotype" and not config.get("independent_pairing_evidence"):
        raise ValueError("Verified donor/haplotype pairing requires evidence")
    checked = []
    ids = set()
    for row in rows:
        array_id = row["array_id"]
        if array_id in ids:
            raise ValueError("Duplicate array ID")
        ids.add(array_id)
        if not set(row["supporting_read_ids"]) <= read_ids or not row["supporting_read_ids"]:
            raise ValueError(f"Declared native read support is absent: {array_id}")
        source = Path(row["reference_fasta"])
        if sha256_file(source) != row["reference_fasta_sha256"]:
            raise ValueError(f"Native source FASTA hash mismatch: {array_id}")
        name, sequence = read_one_fasta(source)
        genomic_interval = row["source_genome_interval"]
        if (not isinstance(genomic_interval, list) or len(genomic_interval) != 2
                or not all(isinstance(value, int) for value in genomic_interval)
                or genomic_interval[0] < 0
                or genomic_interval[1] - genomic_interval[0] != len(sequence)):
            raise ValueError(f"Native source genome interval mismatch: {array_id}")
        left, right, unit = row["array_start0"], row["array_end0"], row["operational_unit_bp"]
        if (name != row["reference_record_id"] or not all(isinstance(x, int) for x in (left, right, unit))
                or unit < 4 or left < 500 or right > len(sequence) - 500 or right <= left
                or (right - left) % unit):
            raise ValueError(f"Invalid array/flank coordinates: {array_id}")
        original_array = sequence[left:right]
        if sha256_sequence(original_array) != row["array_sequence_sha256"]:
            raise ValueError(f"Array sequence mismatch: {array_id}")
        checked.append((row, sequence, original_array, unit, left, right))
    outdir.mkdir(parents=True, exist_ok=False)
    cases = []
    for row, source, original, unit, left, right in checked:
        unit_count = len(original) // unit
        for label, first, after in edit_specs(unit_count):
            deleted_start = left + first * unit
            deleted_end = left + after * unit
            edited = source[:deleted_start] + source[deleted_end:]
            removed = deleted_end - deleted_start
            case_id = f"{row['array_id']}_{label}"
            fasta = outdir / f"{case_id}.fa"
            fasta.write_text(f">{case_id}\n{edited}\n", encoding="ascii")
            if len(edited) != len(source) - removed or edited[:deleted_start] != source[:deleted_start]:
                raise AssertionError("Edit conservation failed")
            cases.append(dict(case_id=case_id, array_id=row["array_id"],
                              edit_type=label.rsplit("_", 1)[0], deletion_percent=round(100 * removed / len(original)),
                              source_reference_fasta=row["reference_fasta"],
                              source_reference_sha256=row["reference_fasta_sha256"],
                              source_record_id=row["reference_record_id"],
                              source_genome_interval=row["source_genome_interval"],
                              source_array_interval=[left, right],
                              source_deleted_interval=[deleted_start, deleted_end],
                              edited_array_interval=[left, right - removed],
                              edited_breakpoint0=deleted_start,
                              original_array_bp=len(original), edited_array_bp=len(original) - removed,
                              injected_deleted_bp=removed,
                              operational_unit_bp=unit,
                              operational_unit_count=unit_count,
                              unit_boundary_status="operational_not_independently_verified_native_monomer_truth",
                              truth_scope="injected_edit_delta_only",
                              native_read_bundle_sha256=config["raw_read_bundle_sha256"],
                              supporting_read_ids=row["supporting_read_ids"],
                              read_pairing_status=config["read_pairing_status"],
                              edited_fasta_sha256=sha256_file(fasta)))
    if sha256_file(read_bundle) != config["raw_read_bundle_sha256"]:
        raise ValueError("Native read bundle changed during edit generation")
    receipt = dict(schema_version=1, split="development", material_id=config["material_id"],
                   config_sha256=sha256_file(config_path), generator_sha256=sha256_file(Path(__file__)),
                   raw_read_bundle=str(read_bundle),
                   raw_read_bundle_sha256=config["raw_read_bundle_sha256"],
                   raw_reads_modified=False, truth_scope="injected_edit_delta_only",
                   read_pairing_status=config["read_pairing_status"],
                   independent_donor_count=config["independent_donor_count"], cases=cases)
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    generate(args.config, args.outdir)


if __name__ == "__main__":
    main()
