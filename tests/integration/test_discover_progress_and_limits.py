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
    assert "discover scan_reads" in result.stderr
    assert "r/min" in result.stderr
    assert "elapsed" in result.stderr
    assert "est" in result.stderr
    assert "rem" in result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "processed_reads=10" in log
    assert "reads_per_second=" in log
    candidate_lines = (outdir / "candidate_reads.tsv").read_text(encoding="utf-8").splitlines()
    assert len(candidate_lines) == 11


def test_discover_avoids_unneeded_full_input_precount(tmp_path: Path) -> None:
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
    assert "discover scan_reads" in result.stderr
    assert "counting" not in result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "input_count_skipped reason=not_required_for_discovery_budget" in log


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


def test_auto_discovery_budget_uses_genome_size_and_round_robin(tmp_path: Path) -> None:
    first = tmp_path / "reads_a.fa"
    second = tmp_path / "reads_b.fa"
    first.write_text(periodic_reads(5, length=600), encoding="utf-8")
    monomer_b = "TTGCAACCTTGGAACCGTTAGGCAATCGTA"
    sequence_b = (monomer_b * 40)[:600]
    second.write_text(
        "".join(f">read_b_{index}\n{sequence_b}\n" for index in range(1, 6)),
        encoding="utf-8",
    )
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(first),
        str(second),
        "--outdir",
        str(outdir),
        "--min-period",
        "18",
        "--max-period",
        "32",
        "--min-support-reads",
        "1",
        "--genome-size",
        "1200",
        "--enable-auto-discovery-budget",
        "--target-discovery-coverage",
        "1.0",
        "--progress-every",
        "1",
    )

    assert result.returncode == 0, result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "auto_discovery_budget enabled reason=genome_size_x_target_coverage" in log
    assert "processed_reads=2 processed_bases=1200" in log
    candidates = (outdir / "candidate_reads.tsv").read_text(encoding="utf-8")
    assert "read_1" in candidates
    assert "read_b_1" in candidates
    families = (outdir / "families.tsv").read_text(encoding="utf-8")
    assert "TXF000001" in families
    assert "TXF000002" in families


def test_genome_size_alone_does_not_enable_auto_discovery_budget(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(periodic_reads(6, length=400), encoding="utf-8")
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
        "1",
        "--genome-size",
        "1200",
        "--progress-every",
        "2",
    )

    assert result.returncode == 0, result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "auto_discovery_budget enabled" not in log
    assert "processed_reads=6 processed_bases=2400" in log
