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


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_empty_fasta_fails_clearly(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fa", "")
    result = run_cli("discover", "--reads", str(reads), "--outdir", str(tmp_path / "discover"))

    assert result.returncode != 0
    assert "empty or contains no records" in result.stderr
    assert not (tmp_path / "discover" / "families.tsv").exists()


def test_invalid_fastq_fails_clearly(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fastq", "@r1\nACGT\n-\n!!!!\n")
    result = run_cli("discover", "--reads", str(reads), "--outdir", str(tmp_path / "discover"))

    assert result.returncode != 0
    assert "Invalid FASTQ separator" in result.stderr


def test_fastq_quality_length_mismatch_fails_clearly(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fq", "@r1\nACGT\n+\n!!!\n")
    result = run_cli("discover", "--reads", str(reads), "--outdir", str(tmp_path / "discover"))

    assert result.returncode != 0
    assert "sequence and quality lengths differ" in result.stderr


def test_low_complexity_short_period_reads_are_flagged(tmp_path: Path) -> None:
    reads = write(
        tmp_path / "reads.fa",
        ">r1\n" + "A" * 500 + "\n"
        ">r2\n" + "AT" * 250 + "\n"
        ">r3\n" + "A" * 500 + "\n",
    )
    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(tmp_path / "discover"),
        "--min-monomer-len",
        "2",
        "--max-monomer-len",
        "20",
        "--min-support-reads",
        "1",
        "--min-repeat-span",
        "100",
    )

    assert result.returncode == 0, result.stderr
    families = (tmp_path / "discover" / "families.tsv").read_text(encoding="utf-8")
    candidates = (tmp_path / "discover" / "candidate_reads.tsv").read_text(encoding="utf-8")
    assert "low_complexity_family" in families
    assert "low_complexity_candidate" in candidates


def test_missing_catalog_fails_before_quantify(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fa", ">r1\nACGTACGT\n")
    missing = tmp_path / "missing.fa"
    result = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalog",
        str(missing),
        "--genome-size",
        "100",
        "--outdir",
        str(tmp_path / "quantify"),
    )

    assert result.returncode != 0
    assert f"Input file does not exist: {missing}" in result.stderr


def test_empty_catalog_fails_quantify(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fa", ">r1\nACGTACGTACGT\n")
    catalog = write(tmp_path / "monomers.fa", "")
    result = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalog",
        str(catalog),
        "--genome-size",
        "100",
        "--outdir",
        str(tmp_path / "quantify"),
    )

    assert result.returncode != 0
    assert "empty or contains no records" in result.stderr


def test_empty_assembly_fails_locate(tmp_path: Path) -> None:
    assembly = write(tmp_path / "assembly.fa", "")
    catalog = write(tmp_path / "monomers.fa", ">family_id=TXF000001;length_bp=12\nACGTACGTACGT\n")
    result = run_cli(
        "locate",
        "--assembly",
        str(assembly),
        "--catalog",
        str(catalog),
        "--outdir",
        str(tmp_path / "locate"),
    )

    assert result.returncode != 0
    assert "empty or contains no records" in result.stderr


def test_kmer_longer_than_read_fails_quantify(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fa", ">r1\nACGT\n")
    catalog = write(tmp_path / "monomers.fa", ">family_id=TXF000001;length_bp=12\nACGTACGTACGT\n")
    result = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalog",
        str(catalog),
        "--genome-size",
        "100",
        "--k",
        "9",
        "--outdir",
        str(tmp_path / "quantify"),
    )

    assert result.returncode != 0
    assert "--k is greater than all read lengths" in result.stderr


def test_existing_outdir_overwrites_known_outputs(tmp_path: Path) -> None:
    reads = write(tmp_path / "reads.fa", ">r1\nACGTACGTACGTACGT\n")
    outdir = tmp_path / "discover"
    outdir.mkdir()
    (outdir / "families.tsv").write_text("old result\n", encoding="utf-8")

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--min-monomer-len",
        "4",
        "--max-monomer-len",
        "8",
        "--min-support-reads",
        "1",
        "--min-repeat-span",
        "8",
    )

    assert result.returncode == 0, result.stderr
    assert "old result" not in (outdir / "families.tsv").read_text(encoding="utf-8")
