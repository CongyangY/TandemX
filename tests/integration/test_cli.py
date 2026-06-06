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


def write_file(path: Path, text: str = "placeholder\n") -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_top_level_help() -> None:
    result = run_cli("--help")
    assert result.returncode == 0
    assert "discover" in result.stdout
    assert "quantify" in result.stdout
    assert "visualize" in result.stdout
    assert "simulate" in result.stdout


def test_subcommand_help() -> None:
    for command in ("discover", "quantify", "locate", "probe", "compare", "visualize", "simulate"):
        result = run_cli(command, "--help")
        assert result.returncode == 0
        assert "usage:" in result.stdout


def test_missing_input_file_errors(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fastq"
    result = run_cli(
        "discover",
        "--reads",
        str(missing),
        "--outdir",
        str(tmp_path / "out"),
    )
    assert result.returncode != 0
    assert f"Input file does not exist: {missing}" in result.stderr


def test_directory_input_path_errors(tmp_path: Path) -> None:
    input_dir = tmp_path / "reads_dir"
    input_dir.mkdir()

    result = run_cli(
        "discover",
        "--reads",
        str(input_dir),
        "--outdir",
        str(tmp_path / "out"),
    )

    assert result.returncode != 0
    assert f"Input path is not a file: {input_dir}" in result.stderr


def test_discover_writes_run_files(tmp_path: Path) -> None:
    reads = write_file(tmp_path / "reads.fa", ">r1\nACGTACGTACGTACGT\n")
    outdir = tmp_path / "discover"

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

    assert result.returncode == 0
    assert "wrote" in result.stdout
    assert outdir.is_dir()
    assert (outdir / "run_config.yaml").is_file()
    assert (outdir / "run.log").is_file()
    assert (outdir / "candidate_reads.tsv").is_file()
    assert (outdir / "monomers.fa").is_file()
    assert (outdir / "families.tsv").is_file()
    config = (outdir / "run_config.yaml").read_text(encoding="utf-8")
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert 'command: "tandemx discover"' in config
    assert 'subcommand: "discover"' in config
    assert 'status: "discover_mvp_completed"' in config
    assert "cwd:" in config
    assert "argv:" in config
    assert "python_version:" in config
    assert "platform:" in config
    assert "func" not in config
    assert "command=tandemx discover" in log
    assert "status=discover_mvp_completed" in log
    assert f"output_directory={outdir}" in log


def test_each_subcommand_writes_run_config(tmp_path: Path) -> None:
    reads = write_file(tmp_path / "reads.fa", ">r1\nACGTACGTACGTACGT\n")
    catalogue = write_file(tmp_path / "repeat_catalogue.tsv")
    monomers = write_file(tmp_path / "monomers.fasta", ">m1\nACGT\n")
    assembly = write_file(tmp_path / "assembly.fa", ">chr1\nACGTACGT\n")
    copy_number = write_file(tmp_path / "copy_number.tsv")
    density = write_file(tmp_path / "repeat_density.tsv")
    probes = write_file(tmp_path / "probe_candidates.tsv")
    comparison = write_file(tmp_path / "copy_number_comparison.tsv")

    quantify_result = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalogue",
        str(monomers),
        "--genome-size",
        "16",
        "--k",
        "4",
        "--outdir",
        str(tmp_path / "quantify"),
    )
    quantify_outdir = tmp_path / "quantify"
    assert quantify_result.returncode == 0, quantify_result.stderr
    assert "copy-number estimates" in quantify_result.stdout
    assert (quantify_outdir / "copy_number.tsv").is_file()
    quantify_config = (quantify_outdir / "run_config.yaml").read_text(encoding="utf-8")
    assert 'command: "tandemx quantify"' in quantify_config
    assert 'subcommand: "quantify"' in quantify_config
    assert 'status: "quantify_mvp_completed"' in quantify_config
    assert "func" not in quantify_config

    locate_result = run_cli(
        "locate",
        "--assembly",
        str(assembly),
        "--catalogue",
        str(monomers),
        "--window-size",
        "4",
        "--step-size",
        "4",
        "--k",
        "4",
        "--outdir",
        str(tmp_path / "locate"),
    )
    locate_outdir = tmp_path / "locate"
    assert locate_result.returncode == 0, locate_result.stderr
    assert "arrays" in locate_result.stdout
    assert (locate_outdir / "repeat_density.bedgraph").is_file()
    assert (locate_outdir / "arrays.bed").is_file()
    assert (locate_outdir / "assembly_vs_read_cn.tsv").is_file()
    locate_config = (locate_outdir / "run_config.yaml").read_text(encoding="utf-8")
    assert 'command: "tandemx locate"' in locate_config
    assert 'subcommand: "locate"' in locate_config
    assert 'status: "locate_mvp_completed"' in locate_config
    assert "func" not in locate_config

    commands = [
        (
            "probe",
            [
                "probe",
                "--catalogue",
                str(catalogue),
                "--monomers",
                str(monomers),
                "--copy-number",
                str(copy_number),
                "--locations",
                str(density),
                "--outdir",
                str(tmp_path / "probe"),
            ],
        ),
        (
            "compare",
            [
                "compare",
                "--read-copy-number",
                str(copy_number),
                "--assembly-density",
                str(density),
                "--outdir",
                str(tmp_path / "compare"),
            ],
        ),
        (
            "visualize",
            [
                "visualize",
                "--catalogue",
                str(catalogue),
                "--copy-number",
                str(copy_number),
                "--locations",
                str(density),
                "--probes",
                str(probes),
                "--comparison",
                str(comparison),
                "--outdir",
                str(tmp_path / "visualize"),
            ],
        ),
    ]

    for command, args in commands:
        result = run_cli(*args)
        outdir = tmp_path / command
        assert result.returncode == 0, result.stderr
        assert "not implemented yet" in result.stdout
        assert (outdir / "run_config.yaml").is_file()
        assert (outdir / "run.log").is_file()
        config = (outdir / "run_config.yaml").read_text(encoding="utf-8")
        assert f'command: "tandemx {command}"' in config
        assert f'subcommand: "{command}"' in config
        assert 'status: "skeleton_not_implemented"' in config
        assert "func" not in config
