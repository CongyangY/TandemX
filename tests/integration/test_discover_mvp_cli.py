from __future__ import annotations

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


def parse_family_lengths(path: Path) -> list[int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [int(line.split("\t")[2]) for line in lines[1:] if line.strip()]


def test_discover_mvp_recovers_toy_family_lengths(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    simulate_result = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "19",
        "--num-reads",
        "80",
        "--read-length",
        "1200",
        "--background-length",
        "2000",
        "--monomer-lengths",
        "566,350",
        "--copies",
        "9,7",
        "--error-rate",
        "0.01",
    )
    assert simulate_result.returncode == 0, simulate_result.stderr

    discover_result = run_cli(
        "discover",
        "--reads",
        str(toy / "reads.fa"),
        "--outdir",
        str(discover),
        "--min-monomer-len",
        "300",
        "--max-monomer-len",
        "700",
        "--min-support-reads",
        "3",
        "--min-repeat-span",
        "600",
    )
    assert discover_result.returncode == 0, discover_result.stderr
    assert "wrote" in discover_result.stdout

    assert (discover / "candidate_reads.tsv").is_file()
    assert (discover / "monomers.fa").is_file()
    assert (discover / "families.tsv").is_file()

    lengths = parse_family_lengths(discover / "families.tsv")
    assert any(abs(length - 566) <= 5 for length in lengths)
    assert any(abs(length - 350) <= 5 for length in lengths)

    candidate_text = (discover / "candidate_reads.tsv").read_text(encoding="utf-8")
    assert "strand=+" not in candidate_text
    assert "\t+\t" in candidate_text
    assert "\t-\t" in candidate_text

    monomers_text = (discover / "monomers.fa").read_text(encoding="utf-8")
    assert ">family_id=TXF" in monomers_text
    assert "length_bp=566" in monomers_text
    assert "length_bp=350" in monomers_text

    config = (discover / "run_config.yaml").read_text(encoding="utf-8")
    assert 'status: "discover_mvp_completed"' in config


def test_discover_runs_de_novo_with_only_reads_and_outdir(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    sequence = "ACGTACGTACGTACGTACGT" * 8
    reads.write_text(
        "".join(f">read_{index}\n{sequence}\n" for index in range(1, 6)),
        encoding="utf-8",
    )
    outdir = tmp_path / "discover"

    result = run_cli("discover", "--reads", str(reads), "--outdir", str(outdir))

    assert result.returncode == 0, result.stderr
    assert (outdir / "candidate_reads.tsv").is_file()
    assert (outdir / "monomers.fa").is_file()
    assert (outdir / "families.tsv").is_file()


def test_discover_collapse_mode_writes_optional_outputs(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    sequence = "ACGTACGTACGTACGTACGT" * 8
    reads.write_text(
        "".join(f">read_{index}\n{sequence}\n" for index in range(1, 6)),
        encoding="utf-8",
    )
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--collapse-redundant-families",
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "collapsed_families.tsv").is_file()
    assert (outdir / "collapsed_monomers.fa").is_file()
    assert (outdir / "family_collapse.tsv").is_file()
    validate = run_cli("validate", "--project", str(outdir))
    assert validate.returncode == 0, validate.stderr


def test_discover_mvp_accepts_fastq(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    sequence = "ACGTACGA" * 8
    reads.write_text(
        "".join(f"@r{index};strand=+\n{sequence}\n+\n{'I' * len(sequence)}\n" for index in range(1, 6)),
        encoding="utf-8",
    )

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(tmp_path / "discover"),
        "--min-monomer-len",
        "8",
        "--max-monomer-len",
        "8",
        "--min-support-reads",
        "2",
        "--min-repeat-span",
        "32",
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "discover" / "families.tsv").is_file()
