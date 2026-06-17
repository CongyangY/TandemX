from __future__ import annotations

import random
import re
import subprocess
import sys
from pathlib import Path


SEEDS = (1, 7, 13, 42, 99)


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


def truth_monomer_lengths(path: Path) -> list[int]:
    lengths = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            continue
        match = re.search(r"length_bp=(\d+)", line)
        if match:
            lengths.append(int(match.group(1)))
    return lengths


def randomized_case(seed: int) -> tuple[tuple[int, int], tuple[int, int]]:
    rng = random.Random(seed)
    first = rng.randint(180, 700)
    second = min(900, first + rng.randint(100, 180))
    lengths = (first, second)
    copies = (rng.randint(8, 11), rng.randint(6, 10))
    return lengths, copies


def test_randomized_toy_workflows_recover_supported_outputs(tmp_path: Path) -> None:
    for seed in SEEDS:
        lengths, copies = randomized_case(seed)
        case = tmp_path / f"seed_{seed}"
        toy = case / "toy"
        discover = case / "discover"
        quantify = case / "quantify"
        locate = case / "locate"
        probe = case / "probe"
        background_length = 3000
        read_length = max(900, (max(lengths) * 2) + 150)
        num_reads = 80
        source_length = background_length + sum(length * copy for length, copy in zip(lengths, copies)) + 100 * len(copies)
        haploid_depth = (num_reads * read_length) / source_length

        simulate_result = run_cli(
            "simulate",
            "toy",
            "--outdir",
            str(toy),
            "--seed",
            str(seed),
            "--num-reads",
            str(num_reads),
            "--read-length",
            str(read_length),
            "--background-length",
            str(background_length),
            "--monomer-lengths",
            ",".join(str(length) for length in lengths),
            "--copies",
            ",".join(str(copy) for copy in copies),
            "--error-rate",
            "0.002",
        )
        assert simulate_result.returncode == 0, simulate_result.stderr

        discover_result = run_cli(
            "discover",
            "--reads",
            str(toy / "reads.fa"),
            "--outdir",
            str(discover),
            "--min-monomer-len",
            str(max(20, min(lengths) - 30)),
            "--max-monomer-len",
            str(max(lengths) + 30),
            "--min-support-reads",
            "3",
            "--min-repeat-span",
            str(max(500, min(lengths))),
        )
        assert discover_result.returncode == 0, f"seed={seed} lengths={lengths} stderr={discover_result.stderr}"

        truth_lengths = truth_monomer_lengths(toy / "truth_monomers.fa")
        assert sorted(truth_lengths) == list(lengths)
        families = parse_tsv(discover / "families.tsv")
        discovered_lengths = [int(row["monomer_length_bp"]) for row in families]
        assert len(families) >= 2, f"seed={seed} discovered={discovered_lengths}"
        for length in truth_lengths:
            assert any(abs(observed - length) <= 10 for observed in discovered_lengths), (
                f"seed={seed} truth={length} discovered={discovered_lengths}"
            )

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
        assert quantify_result.returncode == 0, f"seed={seed} stderr={quantify_result.stderr}"
        copy_rows = parse_tsv(quantify / "copy_number.tsv")
        assert copy_rows
        assert all(float(row["estimated_copy_number"]) > 0 for row in copy_rows)

        locate_result = run_cli(
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
        assert locate_result.returncode == 0, f"seed={seed} stderr={locate_result.stderr}"
        arrays = (locate / "arrays.bed").read_text(encoding="utf-8").splitlines()
        assert arrays
        for line in arrays:
            fields = line.split("\t")
            assert int(fields[1]) >= 0
            assert int(fields[2]) > int(fields[1])

        probe_result = run_cli(
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
        assert probe_result.returncode == 0, f"seed={seed} stderr={probe_result.stderr}"
        probe_rows = parse_tsv(probe / "probes.rank.tsv")
        assert any(float(row["probe_score"]) > 0 and row["warning"] == "" for row in probe_rows)

        validate_result = run_cli("validate", "--project", str(case))
        assert validate_result.returncode == 0, f"seed={seed} stderr={validate_result.stderr}"
