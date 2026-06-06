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


def test_locate_mvp_flags_simulated_collapse(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    quantify = tmp_path / "quantify"
    locate = tmp_path / "locate"
    num_reads = 120
    read_length = 1200
    source_length = 7744
    haploid_depth = (num_reads * read_length) / source_length

    assert run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "31",
        "--num-reads",
        str(num_reads),
        "--read-length",
        str(read_length),
        "--background-length",
        "2000",
        "--monomer-lengths",
        "566,350",
        "--copies",
        "9,7",
        "--error-rate",
        "0.005",
    ).returncode == 0
    assert run_cli(
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
    ).returncode == 0
    assert run_cli(
        "quantify",
        "--reads",
        str(toy / "reads.fa"),
        "--catalog",
        str(discover / "monomers.fa"),
        "--genome-size",
        str(source_length),
        "--haploid-depth",
        f"{haploid_depth:.6f}",
        "--k",
        "21",
        "--outdir",
        str(quantify),
    ).returncode == 0

    result = run_cli(
        "locate",
        "--assembly",
        str(toy / "assembly.fa"),
        "--catalog",
        str(discover / "monomers.fa"),
        "--copy-number",
        str(quantify / "copy_number.tsv"),
        "--window-size",
        "500",
        "--step-size",
        "250",
        "--k",
        "21",
        "--outdir",
        str(locate),
    )
    assert result.returncode == 0, result.stderr
    assert (locate / "repeat_density.bedgraph").is_file()
    assert (locate / "arrays.bed").is_file()
    assert (locate / "assembly_vs_read_cn.tsv").is_file()

    comparisons = parse_tsv(locate / "assembly_vs_read_cn.tsv")
    assert any(row["status"] == "possible_collapse" for row in comparisons)

    for line in (locate / "arrays.bed").read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        start = int(parts[1])
        end = int(parts[2])
        assert start >= 0
        assert end > start
        assert parts[5] in {"+", "-", "."}

    density_line = (locate / "repeat_density.bedgraph").read_text(encoding="utf-8").splitlines()[0]
    assert len(density_line.split("\t")) == 4
