from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_tiny_synthetic_benchmark(tmp_path: Path) -> None:
    outdir = tmp_path / "benchmark"
    result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/run_synthetic_benchmark.py",
            "--config",
            "benchmarks/configs/synthetic_scale.yaml",
            "--scale",
            "tiny",
            "--outdir",
            str(outdir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    benchmark_summary = outdir / "benchmark_summary.tsv"
    accuracy_summary = outdir / "accuracy_summary.tsv"
    assert benchmark_summary.is_file()
    assert accuracy_summary.is_file()

    rows = read_tsv(benchmark_summary)
    assert {row["command"] for row in rows} == {"simulate", "discover", "quantify", "locate", "probe", "validate"}
    assert all(row["exit_status"] == "0" for row in rows)
    assert all(float(row["runtime_seconds"]) > 0 for row in rows)
    assert all(float(row["peak_memory_mb"]) > 0 for row in rows)
    validate_rows = [row for row in rows if row["command"] == "validate"]
    assert validate_rows and validate_rows[0]["output_validated"] == "true"
    assert int(validate_rows[0]["recovered_family_count"]) >= 2
    discover_rows = [row for row in rows if row["command"] == "discover"]
    assert discover_rows[0]["processed_reads"] == "1000"
    assert int(discover_rows[0]["processed_bases"]) > 0
    assert int(discover_rows[0]["candidate_reads"]) > 0
    assert float(discover_rows[0]["reads_per_second"]) > 0
    assert float(discover_rows[0]["mb_per_second"]) > 0
    assert discover_rows[0]["algorithm_mode"] == "spacing_prefilter"

    accuracy_rows = read_tsv(accuracy_summary)
    assert len(accuracy_rows) >= 2
    assert all(row["recovered_closest_length"] for row in accuracy_rows)
    assert all(int(row["length_error_bp"]) == 0 for row in accuracy_rows)
    assert all(float(row["recovered_sequence_identity"]) >= 0.99 for row in accuracy_rows)
    assert all(float(row["copy_number_relative_error"]) < 0.15 for row in accuracy_rows)
    assert all(row["matching_method"] == "sequence_identity_length_aware" for row in accuracy_rows)

    reads = outdir / "tiny" / "simulated" / "reads.fa"
    assert reads.is_file()
