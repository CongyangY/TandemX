from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def periodic_fastq(path: Path, count: int = 8) -> None:
    monomer = "ACGTTCAGGACTAACCGTGA"
    sequence = monomer * 20
    with path.open("wt", encoding="utf-8") as handle:
        for index in range(count):
            handle.write(f"@read_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_inspect_and_real_read_pilot_benchmark(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fq"
    periodic_fastq(reads)
    stats = tmp_path / "read_stats.tsv"

    inspect_result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/inspect_reads.py",
            "--reads",
            str(reads),
            "--output",
            str(stats),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert inspect_result.returncode == 0, inspect_result.stderr
    stats_row = read_tsv(stats)[0]
    assert stats_row["read_count"] == "8"
    assert stats_row["format"] == "fastq"

    outdir = tmp_path / "benchmark"
    benchmark_result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/run_real_read_pilot_benchmark.py",
            "--reads",
            str(reads),
            "--max-reads",
            "4,8",
            "--outdir",
            str(outdir),
            "--min-period",
            "20",
            "--max-period",
            "30",
            "--progress-every",
            "2",
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert benchmark_result.returncode == 0, benchmark_result.stderr
    rows = read_tsv(outdir / "tmpfq_benchmark_summary.tsv")
    assert [row["processed_reads"] for row in rows] == ["4", "8"]
    assert all(row["output_validated"] == "true" for row in rows)
    assert all(float(row["mb_per_second"]) > 0 for row in rows)
    assert all(row["peak_memory_mb"] == "NA" for row in rows)
    assert all(row["backend"] == "python" for row in rows)
