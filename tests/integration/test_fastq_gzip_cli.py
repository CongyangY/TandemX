from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def test_discover_and_quantify_accept_gzipped_fastq_reads(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq.gz"
    repeat = "ACGTACGA" * 8
    with gzip.open(reads, "wt", encoding="utf-8") as handle:
        for index in range(1, 4):
            handle.write(f"@r{index};strand=+\n{repeat}\n+\n{'I' * len(repeat)}\n")

    discover = tmp_path / "discover"
    discover_result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(discover),
        "--min-monomer-len",
        "8",
        "--max-monomer-len",
        "8",
        "--min-support-reads",
        "2",
        "--min-repeat-span",
        "32",
    )
    assert discover_result.returncode == 0, discover_result.stderr
    assert (discover / "monomers.fa").is_file()

    quantify = tmp_path / "quantify"
    quantify_result = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalog",
        str(discover / "monomers.fa"),
        "--genome-size",
        "512",
        "--k",
        "4",
        "--outdir",
        str(quantify),
    )
    assert quantify_result.returncode == 0, quantify_result.stderr
    assert (quantify / "copy_number.tsv").is_file()
