from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from tandemx.io.sequences import (
    SequenceFormatError,
    SequenceStats,
    count_sequence_records_many,
    parse_seqkit_stats_table,
    read_sequence_records,
    read_sequence_records_many,
)


def collect(path: Path) -> list[tuple[str, str, str | None]]:
    return [(record.id, record.sequence, record.quality) for record in read_sequence_records(path)]


def test_reads_fasta_fastq_and_gzip_variants(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fa"
    fastq = tmp_path / "reads.fastq"
    fasta_gz = tmp_path / "reads.fasta.gz"
    fastq_gz = tmp_path / "reads.fq.gz"

    fasta.write_text(">r1 description\nACGT\nACGT\n>r2\nNNAA\n", encoding="utf-8")
    fastq.write_text("@q1 description\nACGT\n+\n!!!!\n@q2\nNNAA\n+\n####\n", encoding="utf-8")
    with gzip.open(fasta_gz, "wt", encoding="utf-8") as handle:
        handle.write(fasta.read_text(encoding="utf-8"))
    with gzip.open(fastq_gz, "wt", encoding="utf-8") as handle:
        handle.write(fastq.read_text(encoding="utf-8"))

    assert collect(fasta) == [("r1", "ACGTACGT", None), ("r2", "NNAA", None)]
    assert collect(fasta_gz) == collect(fasta)
    assert collect(fastq) == [("q1", "ACGT", "!!!!"), ("q2", "NNAA", "####")]
    assert collect(fastq_gz) == collect(fastq)


def test_read_sequence_records_many_streams_files_in_order(tmp_path: Path) -> None:
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">r1\nACGT\n", encoding="utf-8")
    second.write_text(">r2\nTTAA\n", encoding="utf-8")

    observed = [
        (record.id, record.sequence)
        for record in read_sequence_records_many((first, second))
    ]

    assert observed == [("r1", "ACGT"), ("r2", "TTAA")]


def test_read_sequence_records_many_rejects_duplicate_ids_across_files(tmp_path: Path) -> None:
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">r1\nACGT\n", encoding="utf-8")
    second.write_text(">r1\nTTAA\n", encoding="utf-8")

    with pytest.raises(SequenceFormatError, match="Duplicate sequence id across input read files: r1"):
        list(read_sequence_records_many((first, second)))


@pytest.mark.parametrize(
    ("name", "text", "message"),
    [
        ("empty.fa", "", "empty or contains no records"),
        ("bad.fa", "ACGT\n", "sequence before header"),
        ("dup.fa", ">r1\nACGT\n>r1 other\nACGT\n", "Duplicate sequence id"),
        ("bad.fastq", "@r1\nACGT\n-\n!!!!\n", "Invalid FASTQ separator"),
        ("mismatch.fq", "@r1\nACGT\n+\n!!!\n", "sequence and quality lengths differ"),
    ],
)
def test_sequence_reader_reports_clear_format_errors(tmp_path: Path, name: str, text: str, message: str) -> None:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")

    with pytest.raises(SequenceFormatError, match=message):
        list(read_sequence_records(path))


def test_fasta_invalid_base_reports_source_line_after_record_level_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad_base.fa"
    path.write_text(">r1\nACGT\nACXT\n", encoding="utf-8")

    with pytest.raises(SequenceFormatError, match=r"Invalid base.*line 3"):
        list(read_sequence_records(path))


def test_parse_seqkit_stats_table_aggregates_rows() -> None:
    stats = parse_seqkit_stats_table(
        "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
        "a.fa\tFASTA\tDNA\t2\t12\t4\t6.0\t8\n"
        "b.fa\tFASTA\tDNA\t3\t21\t5\t7.0\t10\n"
    )

    assert stats == SequenceStats(record_count=5, total_bases=33, max_read_length=10)


def test_count_sequence_records_many_prefers_seqkit_when_available(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">r1\nACGT\n", encoding="utf-8")
    second.write_text(">r2\nTTAA\n", encoding="utf-8")

    monkeypatch.setattr("tandemx.io.sequences.shutil.which", lambda name: "/usr/bin/seqkit" if name == "seqkit" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        assert command[:4] == ["/usr/bin/seqkit", "stats", "-T", "-j"]
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
                "first.fa\tFASTA\tDNA\t1\t4\t4\t4.0\t4\n"
                "second.fa\tFASTA\tDNA\t1\t4\t4\t4.0\t4\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("tandemx.io.sequences.subprocess.run", fake_run)

    stats = count_sequence_records_many((first, second), threads=6)

    assert stats == SequenceStats(record_count=2, total_bases=8, max_read_length=4)


def test_count_sequence_records_many_uses_rust_before_python_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">r1\nACGT\n", encoding="utf-8")
    second.write_text(">r2\nTTAA\n", encoding="utf-8")

    monkeypatch.setattr("tandemx.io.sequences.shutil.which", lambda name: None)
    monkeypatch.setattr("tandemx.io.sequences.rust_backend_available", lambda: True)

    called = {"rust": False}

    def fake_rust(paths: tuple[Path, ...], *, threads: int) -> SequenceStats:
        called["rust"] = True
        assert paths == (first, second)
        assert threads == 6
        return SequenceStats(record_count=2, total_bases=8, max_read_length=4)

    monkeypatch.setattr("tandemx.io.sequences.rust_count_sequence_paths_stats", fake_rust)

    stats = count_sequence_records_many((first, second), threads=6)

    assert called["rust"] is True
    assert stats == SequenceStats(record_count=2, total_bases=8, max_read_length=4)
