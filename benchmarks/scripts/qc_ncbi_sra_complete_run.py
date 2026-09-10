"""Verify an official NCBI SRA object and its converted FASTQ records."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.download_http_ranges import digest_file as digest_md5_sha256
from tandemx.io.sequences import count_sequence_records_many


def metadata_row(metadata: Path, accession: str) -> dict[str, str]:
    with metadata.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("run_accession") == accession
        ]
    if len(rows) != 1:
        raise ValueError(f"Expected one metadata row for {accession}; found {len(rows)}")
    return rows[0]


def qc_run(
    metadata: Path,
    accession: str,
    sra: Path,
    fastqs: list[Path],
    vdb_validate: Path,
    output_json: Path,
    *,
    threads: int,
    object_type: str = "normalized_sra",
) -> dict[str, object]:
    if output_json.exists():
        raise ValueError(f"Refusing to overwrite existing output: {output_json}")
    if threads < 1 or not fastqs:
        raise ValueError("threads must be positive and at least one FASTQ is required")
    if object_type not in {"normalized_sra", "sra_lite"}:
        raise ValueError("object_type must be normalized_sra or sra_lite")
    row = metadata_row(metadata, accession)
    expected_bytes = int(row[f"{object_type}_bytes"])
    if not sra.is_file() or sra.stat().st_size != expected_bytes:
        observed = sra.stat().st_size if sra.exists() else 0
        raise ValueError(f"SRA size differs from NCBI: {observed}/{expected_bytes}")
    observed_md5, observed_sha256 = digest_md5_sha256(sra)
    if observed_md5 != row[f"{object_type}_md5"]:
        raise ValueError("SRA MD5 differs from official NCBI metadata")
    validation = subprocess.run(
        [str(vdb_validate), str(sra)],
        check=False,
        capture_output=True,
        text=True,
    )
    if validation.returncode != 0:
        raise ValueError(
            f"vdb-validate failed ({validation.returncode}): "
            f"{(validation.stderr or validation.stdout).strip()}"
        )
    missing = [str(path) for path in fastqs if not path.is_file()]
    if missing:
        raise ValueError(f"Converted FASTQ files are missing: {missing}")
    stats = count_sequence_records_many(fastqs, threads=threads)
    expected_records = int(row["sequence_record_count"])
    expected_bases = int(row["base_count"])
    if stats.record_count != expected_records or stats.total_bases != expected_bases:
        raise ValueError(
            f"Converted FASTQ statistics differ from NCBI: records "
            f"{stats.record_count}/{expected_records}; bases {stats.total_bases}/{expected_bases}"
        )
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "complete",
        "run_accession": accession,
        "source_object_type": object_type,
        "source_sra_bytes": expected_bytes,
        "source_sra_md5": observed_md5,
        "source_sra_sha256": observed_sha256,
        "vdb_validate_passed": True,
        "converted_FASTQ_files": [
            {"path": str(path), "bytes": path.stat().st_size} for path in fastqs
        ],
        "FASTQ_record_count": stats.record_count,
        "total_bases": stats.total_bases,
        "max_read_length": stats.max_read_length,
        "NCBI_record_count_matched": True,
        "NCBI_base_count_matched": True,
        "metadata_sha256": digest_file(metadata),
        "script_sha256": digest_file(Path(__file__)),
        "evidence_limit": (
            "official_SRA_Lite_MD5_plus_vdb_validate_and_converted_FASTQ_aggregate_counts;"
            "original_quality_scores_are_not_present;FASTQ_files_are_derived_not_publisher_objects"
            if object_type == "sra_lite"
            else "official_normalized_SRA_MD5_plus_vdb_validate_and_converted_FASTQ_aggregate_counts;"
            "FASTQ_files_are_derived_not_publisher_objects"
        ),
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sra", required=True, type=Path)
    parser.add_argument("--fastq", required=True, action="append", type=Path)
    parser.add_argument("--vdb-validate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--object-type",
        choices=("normalized_sra", "sra_lite"),
        default="normalized_sra",
    )
    args = parser.parse_args()
    qc_run(
        args.metadata,
        args.accession,
        args.sra,
        args.fastq,
        args.vdb_validate,
        args.output,
        threads=args.threads,
        object_type=args.object_type,
    )


if __name__ == "__main__":
    main()
