"""Verify a complete ENA run by file hashes and native FASTQ statistics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.fetch_ena_complete import files_for_run, verify_file
from tandemx.io.sequences import count_sequence_records_many


def qc_run(
    metadata: Path,
    accession: str,
    fastq_dir: Path,
    output_json: Path,
    *,
    threads: int,
) -> dict[str, object]:
    if output_json.exists():
        raise ValueError(f"Refusing to overwrite existing output: {output_json}")
    if threads < 1:
        raise ValueError("threads must be positive")
    row, files = files_for_run(metadata, accession)
    paths = [fastq_dir / entry["filename"] for entry in files]
    verified_files = []
    for path, entry in zip(paths, files):
        if not path.is_file():
            raise ValueError(f"Complete ENA file is missing: {path}")
        verified_files.append(
            {
                "filename": path.name,
                **verify_file(path, entry["expected_bytes"], entry["expected_md5"]),
            }
        )
    stats = count_sequence_records_many(paths, threads=threads)
    expected_reads = int(row["read_count"])
    expected_bases = int(row["base_count"])
    layout = row.get("library_layout", "").strip().upper()
    read_count_interpretation = "ENA_read_count_equals_FASTQ_records"
    read_count_reconciled = stats.record_count == expected_reads
    if not read_count_reconciled and layout == "PAIRED" and stats.record_count == 2 * expected_reads:
        read_count_reconciled = True
        read_count_interpretation = "ENA_read_count_equals_spots;two_FASTQ_records_per_spot"
    if not read_count_reconciled or stats.total_bases != expected_bases:
        raise ValueError(
            f"FASTQ statistics differ from ENA: FASTQ records {stats.record_count}/"
            f"ENA read_count field {expected_reads} (layout={layout or 'unknown'}); "
            f"bases {stats.total_bases}/{expected_bases}"
        )
    result: dict[str, object] = {
        "schema_version": 2,
        "status": "complete",
        "run_accession": accession,
        "file_count": len(paths),
        "FASTQ_record_count": stats.record_count,
        "ENA_read_count_field": expected_reads,
        "read_count_interpretation": read_count_interpretation,
        "total_bases": stats.total_bases,
        "max_read_length": stats.max_read_length,
        "ENA_read_or_spot_count_reconciled": True,
        "ENA_base_count_matched": True,
        "all_ENA_file_MD5_matched": True,
        "files": verified_files,
        "metadata_sha256": digest_file(metadata),
        "script_sha256": digest_file(Path(__file__)),
        "evidence_limit": "official_file_MD5_plus_FASTQ_structure_and_aggregate_counts;paired_read_names_and_global_duplicate_IDs_not_recomputed",
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    qc_run(
        args.metadata,
        args.accession,
        args.fastq_dir,
        args.output,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
