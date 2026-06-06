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


def test_probe_mvp_ranks_toy_candidates(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    quantify = tmp_path / "quantify"
    locate = tmp_path / "locate"
    probe = tmp_path / "probe"
    source_length = 7744
    haploid_depth = (120 * 1200) / source_length

    assert run_cli("simulate", "toy", "--outdir", str(toy), "--seed", "41", "--num-reads", "120", "--error-rate", "0.005").returncode == 0
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
        "--outdir",
        str(quantify),
    ).returncode == 0
    assert run_cli(
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
        "--outdir",
        str(locate),
    ).returncode == 0

    result = run_cli(
        "probe",
        "--catalog",
        str(discover / "monomers.fa"),
        "--assembly",
        str(toy / "assembly.fa"),
        "--copy-number",
        str(quantify / "copy_number.tsv"),
        "--arrays",
        str(locate / "arrays.bed"),
        "--min-len",
        "80",
        "--max-len",
        "300",
        "--outdir",
        str(probe),
    )

    assert result.returncode == 0, result.stderr
    assert (probe / "probes.fa").is_file()
    assert (probe / "probes.rank.tsv").is_file()
    assert (probe / "in_silico_fish.tsv").is_file()

    ranked = parse_tsv(probe / "probes.rank.tsv")
    assert ranked
    scores = [float(row["probe_score"]) for row in ranked]
    assert scores == sorted(scores, reverse=True)
    assert all(float(row["gc_content"]) < 0.8 for row in ranked)

    signals = parse_tsv(probe / "in_silico_fish.tsv")
    assert signals
    assert all(int(row["end"]) > int(row["start"]) for row in signals)

    config = (probe / "run_config.yaml").read_text(encoding="utf-8")
    assert 'status: "probe_mvp_completed"' in config
