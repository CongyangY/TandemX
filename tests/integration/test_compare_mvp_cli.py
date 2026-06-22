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


def parse_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def write_copy_number(path: Path) -> None:
    path.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
        "haploid_depth\testimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF_consistent\t100\t10\t10\t1\t10\t1000\thigh\t\n"
        "TXF_collapse\t100\t10\t10\t1\t10\t1000\thigh\t\n"
        "TXF_over\t100\t10\t10\t1\t10\t1000\thigh\t\n"
        "TXF_reads_only\t100\t10\t10\t1\t10\t1000\thigh\t\n",
        encoding="utf-8",
    )


def write_arrays(path: Path) -> None:
    path.write_text(
        "chr1\t0\t900\tTXF_consistent\t900\t.\thigh\t\n"
        "chr1\t1000\t1300\tTXF_collapse\t900\t.\thigh\t\n"
        "chr1\t2000\t3800\tTXF_over\t900\t.\thigh\t\n"
        "chr1\t4000\t4500\tTXF_assembly_only\t900\t.\thigh\t\n",
        encoding="utf-8",
    )


def test_compare_mvp_writes_statuses_and_validates(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    arrays = tmp_path / "arrays.bed"
    outdir = tmp_path / "compare"
    write_copy_number(copy_number)
    write_arrays(arrays)

    result = run_cli(
        "compare",
        "--copy-number",
        str(copy_number),
        "--arrays",
        str(arrays),
        "--outdir",
        str(outdir),
    )

    assert result.returncode == 0, result.stderr
    output = outdir / "assembly_vs_read_cn.tsv"
    assert output.is_file()
    rows = {row["family_id"]: row for row in parse_tsv(output)}
    assert rows["TXF_consistent"]["status"] == "consistent"
    assert rows["TXF_collapse"]["status"] == "possible_collapse"
    assert rows["TXF_over"]["status"] == "possible_overexpansion"
    assert rows["TXF_reads_only"]["status"] == "reads_only"
    assert rows["TXF_assembly_only"]["status"] == "assembly_only"
    assert rows["TXF_assembly_only"]["warning"] == "missing_read_estimate"
    assert run_cli("validate", "--project", str(outdir)).returncode == 0


def test_compare_mvp_missing_copy_number_errors_clearly(tmp_path: Path) -> None:
    arrays = tmp_path / "arrays.bed"
    write_arrays(arrays)

    result = run_cli(
        "compare",
        "--copy-number",
        str(tmp_path / "missing_copy_number.tsv"),
        "--arrays",
        str(arrays),
        "--outdir",
        str(tmp_path / "compare"),
    )

    assert result.returncode != 0
    assert "Input file does not exist" in result.stderr
    assert "missing_copy_number.tsv" in result.stderr


def test_compare_mvp_missing_arrays_errors_clearly(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    write_copy_number(copy_number)

    result = run_cli(
        "compare",
        "--copy-number",
        str(copy_number),
        "--arrays",
        str(tmp_path / "missing_arrays.bed"),
        "--outdir",
        str(tmp_path / "compare"),
    )

    assert result.returncode != 0
    assert "Input file does not exist" in result.stderr
    assert "missing_arrays.bed" in result.stderr


def test_compare_mvp_rejects_assembly_density_as_family_input(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    density = tmp_path / "repeat_density.bedgraph"
    write_copy_number(copy_number)
    density.write_text("chr1\t0\t100\t0.5\n", encoding="utf-8")

    result = run_cli(
        "compare",
        "--copy-number",
        str(copy_number),
        "--assembly-density",
        str(density),
        "--outdir",
        str(tmp_path / "compare"),
    )

    assert result.returncode != 0
    assert "--arrays is required" in result.stderr
    assert "does not include family_id" in result.stderr
