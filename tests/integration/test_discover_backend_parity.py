from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest
import tandemx

from tandemx.cli import main
from tandemx.discover.rust_backend import rust_backend_available


pytestmark = pytest.mark.skipif(
    not rust_backend_available(),
    reason="TandemX Rust extension is not installed",
)


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_python_and_rust_discover_backend_parity(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    simulated = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "909",
        "--num-reads",
        "300",
        "--read-length",
        "1800",
        "--monomer-lengths",
        "421,729",
        "--copies",
        "8,6",
        "--error-rate",
        "0.005",
    )
    assert simulated.returncode == 0, simulated.stderr

    candidate_counts: dict[str, int] = {}
    family_lengths: dict[str, list[int]] = {}
    for backend in ("python", "rust"):
        outdir = tmp_path / backend
        result = run_cli(
            "discover",
            "--reads",
            str(toy / "reads.fa"),
            "--outdir",
            str(outdir),
            "--min-period",
            "300",
            "--max-period",
            "760",
            "--min-repeat-span",
            "800",
            "--min-support-reads",
            "3",
            "--kmer-backend",
            backend,
        )
        assert result.returncode == 0, result.stderr
        candidate_counts[backend] = len(read_rows(outdir / "candidate_reads.tsv"))
        family_lengths[backend] = sorted(
            int(row["monomer_length_bp"])
            for row in read_rows(outdir / "families.tsv")
        )

    difference = abs(candidate_counts["python"] - candidate_counts["rust"])
    assert difference / max(candidate_counts.values()) <= 0.20
    assert len(family_lengths["python"]) == len(family_lengths["rust"])
    assert all(
        abs(python_length - rust_length) <= 5
        for python_length, rust_length in zip(
            family_lengths["python"], family_lengths["rust"]
        )
    )


def test_rust_cli_reports_unavailable_extension(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read_1\nACGTTCAGGACTAACCGTGAACGTTCAGGACTAACCGTGA\n", encoding="utf-8")
    outdir = tmp_path / "discover"
    monkeypatch.setitem(sys.modules, "tandemx._rust_core", None)
    monkeypatch.delattr(tandemx, "_rust_core", raising=False)

    with pytest.raises(SystemExit) as error:
        main(
            [
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
                "--min-repeat-span",
                "20",
                "--kmer-backend",
                "rust",
            ]
        )

    assert error.value.code == 2
    assert "Rust backend is unavailable" in capsys.readouterr().err
    assert "rust_backend_unavailable=" in (outdir / "run.log").read_text(encoding="utf-8")
