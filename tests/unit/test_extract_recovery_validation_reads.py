from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from benchmarks.scripts.extract_recovery_validation_reads import extract
from tandemx.recovery.poc import sha256


def _fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for identifier, sequence in records:
            handle.write(f"@{identifier} comment\n{sequence}\n+\n{'I' * len(sequence)}\n")


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    evidence = tmp_path / "recruited_reads.tsv"
    evidence.write_text(
        "read_id\ttarget\tlocus_id\tfamily_id\tevidence_type\n"
        "candidate\tleft\tL0004\tF1\tflank_anchored\n"
        "support1\tright\tL0004\tF1\tflank_anchored\n"
        "support1\tleft\tL0004\tF1\tflank_anchored\n"
        "support2\tleft\tL0004\tF1\tflank_anchored\n",
        encoding="utf-8",
    )
    source = tmp_path / "reads.fastq.gz"
    _fastq(source, [("candidate", "AAAA"), ("support1", "CCCC"), ("support2", "GGGG"), ("other", "TTTT")])
    qc = tmp_path / "qc.json"
    qc.write_text(json.dumps({"complete": True, "fastq_records_valid": True, "input_sha256": sha256(source)}), encoding="utf-8")
    return evidence, source, qc


def test_extracts_distinct_non_candidate_support_reads_with_hash_receipt(tmp_path: Path) -> None:
    evidence, source, qc = _inputs(tmp_path)
    result = extract(evidence, source, qc, tmp_path / "out", "L0004", "candidate")

    assert result["complete"] is True
    assert result["extracted_independent_read_ids"] == 2
    assert (tmp_path / "out" / "reads.fa").read_text(encoding="ascii") == ">support1\nCCCC\n>support2\nGGGG\n"
    assert (tmp_path / "out" / "IDs.tsv").read_text(encoding="utf-8") == "read_id\nsupport1\nsupport2\n"
    receipt = json.loads((tmp_path / "out" / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["source_fastq_sha256"] == sha256(source)


def test_missing_or_duplicate_source_id_keeps_partials_and_records_failure(tmp_path: Path) -> None:
    evidence, source, qc = _inputs(tmp_path)
    _fastq(source, [("candidate", "AAAA"), ("support1", "CCCC"), ("support1", "CCCC")])
    qc.write_text(json.dumps({"complete": True, "fastq_records_valid": True, "input_sha256": sha256(source)}), encoding="utf-8")
    outdir = tmp_path / "out"

    with pytest.raises(ValueError, match="missing IDs: support2; duplicate IDs: support1"):
        extract(evidence, source, qc, outdir, "L0004", "candidate")
    assert (outdir / "reads.fa.partial").is_file()
    assert not (outdir / "reads.fa").exists()
    assert json.loads((outdir / "receipt.json").read_text(encoding="utf-8"))["complete"] is False
