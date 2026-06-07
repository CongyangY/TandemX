from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


CORE_ALGORITHM_FILES = [
    Path("tandemx/discover/mvp.py"),
    Path("tandemx/quantify/mvp.py"),
    Path("tandemx/locate/mvp.py"),
    Path("tandemx/probe/mvp.py"),
]


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


def parse_truth_monomer_lengths(path: Path) -> list[int]:
    lengths = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            continue
        match = re.search(r"length_bp=(\d+)", line)
        if match:
            lengths.append(int(match.group(1)))
    return lengths


def read_fasta_sequences(path: Path) -> list[str]:
    sequences = []
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if parts:
                sequences.append("".join(parts))
            parts = []
        else:
            parts.append(line.strip())
    if parts:
        sequences.append("".join(parts))
    return sequences


def max_base_fraction(sequence: str) -> float:
    return max(sequence.count(base) for base in "ACGT") / len(sequence)


def test_core_algorithm_files_do_not_hardcode_default_toy_lengths() -> None:
    forbidden = re.compile(r"(?<!\d)(566|350)(?!\d)")
    for path in CORE_ALGORITHM_FILES:
        text = path.read_text(encoding="utf-8")
        assert not forbidden.search(text), f"{path} contains a hardcoded default toy length"


def test_non_default_repeat_lengths_end_to_end(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    discover = tmp_path / "discover"
    quantify = tmp_path / "quantify"
    locate = tmp_path / "locate"
    probe = tmp_path / "probe"

    num_reads = 180
    read_length = 1800
    background_length = 2500
    monomer_lengths = (421, 729)
    copies = (8, 6)
    source_length = background_length + sum(
        length * copy for length, copy in zip(monomer_lengths, copies)
    ) + 100 * len(copies)
    haploid_depth = (num_reads * read_length) / source_length

    simulate_result = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "77",
        "--num-reads",
        str(num_reads),
        "--read-length",
        str(read_length),
        "--background-length",
        str(background_length),
        "--monomer-lengths",
        "421,729",
        "--copies",
        "8,6",
        "--error-rate",
        "0.005",
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
        "900",
        "--min-support-reads",
        "3",
        "--min-repeat-span",
        "800",
    )
    assert discover_result.returncode == 0, discover_result.stderr

    truth_lengths = parse_truth_monomer_lengths(toy / "truth_monomers.fa")
    assert sorted(truth_lengths) == [421, 729]
    families = parse_tsv(discover / "families.tsv")
    discovered_lengths = [int(row["monomer_length_bp"]) for row in families]
    assert len(families) >= 2
    for truth_length in truth_lengths:
        assert any(abs(length - truth_length) <= 10 for length in discovered_lengths)

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
    copy_rows = parse_tsv(quantify / "copy_number.tsv")
    assert len(copy_rows) >= 2
    for row in copy_rows:
        assert int(row["diagnostic_kmer_count"]) > 0
        assert float(row["estimated_copy_number"]) > 0
        assert "confidence" in row
        assert "warning" in row

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
    assert locate_result.returncode == 0, locate_result.stderr
    assert (locate / "repeat_density.bedgraph").is_file()
    assert (locate / "arrays.bed").is_file()
    assert (locate / "assembly_vs_read_cn.tsv").is_file()

    array_lines = (locate / "arrays.bed").read_text(encoding="utf-8").splitlines()
    assert array_lines
    for line in array_lines:
        fields = line.split("\t")
        start = int(fields[1])
        end = int(fields[2])
        assert start >= 0
        assert end > start
    comparison_rows = parse_tsv(locate / "assembly_vs_read_cn.tsv")
    statuses = {row["status"] for row in comparison_rows}
    warnings = {row["warning"] for row in comparison_rows}
    assert "possible_collapse" in statuses or any(warning for warning in warnings)

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
    assert probe_result.returncode == 0, probe_result.stderr
    assert (probe / "probes.fa").is_file()
    assert (probe / "probes.rank.tsv").is_file()
    assert (probe / "in_silico_fish.tsv").is_file()

    probe_rows = parse_tsv(probe / "probes.rank.tsv")
    assert probe_rows
    scores = [float(row["probe_score"]) for row in probe_rows]
    assert scores == sorted(scores, reverse=True)
    probe_sequences = read_fasta_sequences(probe / "probes.fa")
    assert probe_sequences
    assert all(max_base_fraction(sequence) < 0.8 for sequence in probe_sequences)
