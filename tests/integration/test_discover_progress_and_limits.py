from __future__ import annotations

import gzip
import subprocess
import sys
import time
from pathlib import Path


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def periodic_reads(count: int, length: int = 400) -> str:
    monomer = "ACGTTCAGGACTAACCGTGA"
    sequence = (monomer * ((length // len(monomer)) + 1))[:length]
    return "".join(f">read_{index}\n{sequence}\n" for index in range(1, count + 1))


def test_max_reads_progress_and_incremental_files(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(periodic_reads(40), encoding="utf-8")
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--min-period",
        "20",
        "--max-period",
        "30",
        "--min-support-reads",
        "2",
        "--max-reads",
        "10",
        "--progress-every",
        "5",
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "run.log").is_file()
    assert (outdir / "candidate_reads.tsv").is_file()
    assert "discover | scan_reads" in result.stderr
    assert "reads/min" in result.stderr
    assert "elapsed" in result.stderr
    assert "total est" in result.stderr
    assert "remaining" in result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "processed_reads=10" in log
    assert "reads_per_second=" in log
    candidate_lines = (outdir / "candidate_reads.tsv").read_text(encoding="utf-8").splitlines()
    assert len(candidate_lines) == 11


def test_discover_precounts_reads_for_bounded_progress_without_limits(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(periodic_reads(12), encoding="utf-8")
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--min-period",
        "20",
        "--max-period",
        "30",
        "--min-support-reads",
        "2",
        "--progress-every",
        "6",
        "--count-threads",
        "1",
    )

    assert result.returncode == 0, result.stderr
    assert "discover | count_inputs" in result.stderr
    assert "/12 reads" in result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "input_summary read_files=1 total_reads=12" in log


def test_max_read_bases_stops_before_exceeding_limit(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(periodic_reads(10, length=400), encoding="utf-8")
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--min-period",
        "20",
        "--max-period",
        "30",
        "--min-support-reads",
        "2",
        "--max-read-bases",
        "1200",
        "--progress-every",
        "1",
    )

    assert result.returncode == 0, result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "processed_reads=3 processed_bases=1200" in log
    assert "limit_reached=max_read_bases" in log


def test_run_log_and_candidate_table_exist_while_discover_is_running(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(periodic_reads(3000, length=800), encoding="utf-8")
    outdir = tmp_path / "discover"
    process = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "tandemx.cli",
            "discover",
            "--reads",
            str(reads),
            "--outdir",
            str(outdir),
            "--min-period",
            "20",
            "--max-period",
            "30",
            "--min-support-reads",
            "2",
            "--progress-every",
            "100",
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    observed_while_running = False
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline and process.poll() is None:
        if (outdir / "run.log").is_file() and (outdir / "candidate_reads.tsv").is_file():
            observed_while_running = True
            process.terminate()
            break
        time.sleep(0.01)
    stdout, stderr = process.communicate(timeout=10)

    assert observed_while_running, (stdout, stderr)
    assert (outdir / "candidate_reads.tsv").read_text(encoding="utf-8").startswith("read_id\t")


def test_gzip_fastq_discover_and_kmer_too_long_warning(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq.gz"
    monomer = "ACGTTCAGGACTAACCGTGA"
    sequence = monomer * 8
    with gzip.open(reads, "wt", encoding="utf-8") as handle:
        for index in range(1, 6):
            handle.write(f"@read_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")

    valid_outdir = tmp_path / "valid"
    valid = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(valid_outdir),
        "--min-period",
        "20",
        "--max-period",
        "30",
        "--min-support-reads",
        "2",
    )
    assert valid.returncode == 0, valid.stderr

    invalid_outdir = tmp_path / "invalid"
    invalid = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(invalid_outdir),
        "--kmer-size",
        "1000",
    )
    assert invalid.returncode != 0
    assert "No reads were long enough for --kmer-size 1000" in invalid.stderr
    assert "skipped_short_kmer=5" in (invalid_outdir / "run.log").read_text(encoding="utf-8")
