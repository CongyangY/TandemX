from __future__ import annotations

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


def test_discover_1000_read_pilot_runtime(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    simulated = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "515",
        "--num-reads",
        "1000",
        "--read-length",
        "1800",
        "--background-length",
        "2500",
        "--monomer-lengths",
        "421,729",
        "--copies",
        "8,6",
        "--error-rate",
        "0.005",
    )
    assert simulated.returncode == 0, simulated.stderr

    started = time.perf_counter()
    result = run_cli(
        "discover",
        "--reads",
        str(toy / "reads.fa"),
        "--outdir",
        str(discover),
        "--min-period",
        "300",
        "--max-period",
        "760",
        "--min-support-reads",
        "5",
        "--min-repeat-span",
        "800",
        "--progress-every",
        "250",
    )
    runtime = time.perf_counter() - started

    assert result.returncode == 0, result.stderr
    assert runtime < 30
    family_lines = (discover / "families.tsv").read_text(encoding="utf-8").splitlines()
    assert len(family_lines) >= 3
    log = (discover / "run.log").read_text(encoding="utf-8")
    assert "algorithm_mode=spacing_prefilter" in log
    assert "processed_reads=1000" in log
