from __future__ import annotations

import hashlib
from pathlib import Path

from benchmarks.scripts.qc_ncbi_sra_complete_run import qc_run


def test_qc_ncbi_sra_matches_object_and_converted_counts(tmp_path: Path) -> None:
    sra = tmp_path / "ERR1234567"
    sra.write_bytes(b"SRA-test-object")
    fastqs = [tmp_path / "reads_1.fastq", tmp_path / "reads_2.fastq"]
    for index, path in enumerate(fastqs, start=1):
        path.write_text(f"@r/{index}\nACGT\n+\nIIII\n", encoding="utf-8")
    validator = tmp_path / "vdb-validate"
    validator.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    validator.chmod(0o755)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "run_accession\tnormalized_sra_bytes\tnormalized_sra_md5\t"
        "sequence_record_count\tbase_count\n"
        f"ERR1234567\t{sra.stat().st_size}\t"
        f"{hashlib.md5(sra.read_bytes()).hexdigest()}\t2\t8\n",
        encoding="utf-8",
    )
    result = qc_run(
        metadata,
        "ERR1234567",
        sra,
        fastqs,
        validator,
        tmp_path / "qc.json",
        threads=2,
    )
    assert result["status"] == "complete"
    assert result["vdb_validate_passed"] is True
    assert result["FASTQ_record_count"] == 2
    assert result["total_bases"] == 8


def test_qc_ncbi_sra_lite_records_quality_boundary(tmp_path: Path) -> None:
    sra = tmp_path / "SRR1234567.lite.1"
    sra.write_bytes(b"SRA-lite-test-object")
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")
    validator = tmp_path / "vdb-validate"
    validator.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    validator.chmod(0o755)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "run_accession\tsra_lite_bytes\tsra_lite_md5\tsequence_record_count\tbase_count\n"
        f"SRR1234567\t{sra.stat().st_size}\t"
        f"{hashlib.md5(sra.read_bytes()).hexdigest()}\t1\t4\n",
        encoding="utf-8",
    )
    result = qc_run(
        metadata,
        "SRR1234567",
        sra,
        [fastq],
        validator,
        tmp_path / "qc.json",
        threads=1,
        object_type="sra_lite",
    )
    assert result["source_object_type"] == "sra_lite"
    assert "original_quality_scores_are_not_present" in result["evidence_limit"]
