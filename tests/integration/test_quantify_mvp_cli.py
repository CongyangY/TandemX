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


def test_quantify_mvp_estimates_toy_copy_number(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    quantify = tmp_path / "quantify"
    num_reads = 120
    read_length = 1200
    background_length = 2000
    copies = (9, 7)
    monomer_lengths = (566, 350)
    spacer_bp = 100 * len(copies)
    source_length = background_length + sum(length * copy for length, copy in zip(monomer_lengths, copies)) + spacer_bp
    haploid_depth = (num_reads * read_length) / source_length

    assert run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "23",
        "--num-reads",
        str(num_reads),
        "--read-length",
        str(read_length),
        "--background-length",
        str(background_length),
        "--monomer-lengths",
        "566,350",
        "--copies",
        "9,7",
        "--error-rate",
        "0.005",
    ).returncode == 0
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

    quantify_result = run_cli(
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
    )
    assert quantify_result.returncode == 0, quantify_result.stderr
    assert (quantify / "copy_number.tsv").is_file()

    estimates = parse_tsv(quantify / "copy_number.tsv")
    truth = {row["monomer_length_bp"]: row for row in parse_tsv(toy / "truth_copy_number.tsv")}
    assert len(estimates) >= 2
    for estimate in estimates:
        length = estimate["monomer_length"]
        if length not in truth:
            continue
        observed = float(estimate["estimated_copy_number"])
        expected = float(truth[length]["read_copies"])
        assert abs(observed - expected) / expected <= 0.35
        assert float(estimate["estimated_bp"]) > 0

    config = (quantify / "run_config.yaml").read_text(encoding="utf-8")
    assert 'status: "quantify_mvp_completed"' in config
