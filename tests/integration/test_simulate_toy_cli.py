from __future__ import annotations

import hashlib
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


def file_hashes(directory: Path) -> dict[str, str]:
    hashes = {}
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.name not in {"run_config.yaml", "run.log"}:
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def test_simulate_help() -> None:
    result = run_cli("simulate", "toy", "--help")
    assert result.returncode == 0
    assert "--monomer-lengths" in result.stdout
    assert "--error-rate" in result.stdout


def test_simulate_toy_cli_outputs_and_reproducibility(tmp_path: Path) -> None:
    out1 = tmp_path / "toy1"
    out2 = tmp_path / "toy2"
    args = [
        "simulate",
        "toy",
        "--seed",
        "123",
        "--num-reads",
        "12",
        "--read-length",
        "900",
        "--background-length",
        "1500",
        "--monomer-lengths",
        "566,350",
        "--copies",
        "6,4",
        "--error-rate",
        "0.02",
    ]

    result1 = run_cli(*args, "--outdir", str(out1))
    result2 = run_cli(*args, "--outdir", str(out2))

    assert result1.returncode == 0, result1.stderr
    assert result2.returncode == 0, result2.stderr
    assert "wrote toy dataset" in result1.stdout

    expected = {
        "reads.fa",
        "assembly.fa",
        "truth_monomers.fa",
        "truth_arrays.bed",
        "truth_copy_number.tsv",
        "simulation_config.yaml",
        "run_config.yaml",
        "run.log",
    }
    assert expected == {path.name for path in out1.iterdir()}
    assert file_hashes(out1) == file_hashes(out2)

    reads = (out1 / "reads.fa").read_text(encoding="utf-8")
    assert "strand=+" in reads
    assert "strand=-" in reads

    truth = (out1 / "truth_copy_number.tsv").read_text(encoding="utf-8")
    assert "TXF000001\tTXM000001\t566" in truth
    assert "TXF000002\tTXM000002\t350" in truth
    assert "simulated_under_assembly" in truth

    config = (out1 / "simulation_config.yaml").read_text(encoding="utf-8")
    assert 'command: "tandemx simulate toy"' in config
    assert "seed: 123" in config
    assert "monomer_length_bp: 566" in config


def test_simulate_toy_rejects_mismatched_lengths_and_copies(tmp_path: Path) -> None:
    result = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(tmp_path / "bad"),
        "--monomer-lengths",
        "566,350",
        "--copies",
        "6",
    )
    assert result.returncode != 0
    assert "must have the same number of values" in result.stderr
