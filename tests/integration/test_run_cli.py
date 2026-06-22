from __future__ import annotations

import csv
import json
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


def simulate_toy(tmp_path: Path) -> Path:
    toy = tmp_path / "toy"
    result = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "404",
        "--num-reads",
        "40",
        "--read-length",
        "1200",
    )
    assert result.returncode == 0, result.stderr
    return toy


def read_summary(outdir: Path) -> list[dict[str, str]]:
    with (outdir / "pipeline_summary.tsv").open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_manifest(outdir: Path) -> list[dict[str, str]]:
    with (outdir / "output_manifest.tsv").open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reads_only_run_writes_summaries(tmp_path: Path) -> None:
    toy = simulate_toy(tmp_path)
    outdir = tmp_path / "reads_only"
    result = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(outdir),
        "--steps",
        "discover,quantify,validate",
        "--kmer-backend",
        "rust",
        "--profile",
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "discover" / "monomers.fa").is_file()
    assert (outdir / "quantify" / "copy_number.tsv").is_file()
    assert (outdir / "pipeline_summary.json").is_file()
    assert (outdir / "output_manifest.tsv").is_file()
    assert (outdir / "run_report.md").is_file()
    assert "output_manifest.tsv:" in result.stdout
    assert "run_report.md:" in result.stdout
    assert [row["step"] for row in read_summary(outdir)] == ["discover", "quantify", "validate"]
    assert len(json.loads((outdir / "pipeline_summary.json").read_text(encoding="utf-8"))) == 3
    assert (outdir / "profiles" / "discover.prof").is_file()
    quantify_command = read_summary(outdir)[1]["command"]
    assert "--max-reads" not in quantify_command
    manifest_paths = {Path(row["file_path"]).relative_to(outdir).as_posix() for row in read_manifest(outdir)}
    assert "discover/candidate_reads.tsv" in manifest_paths
    assert "discover/monomers.fa" in manifest_paths
    assert "discover/families.tsv" in manifest_paths
    assert "discover/family_similarity.tsv" in manifest_paths
    assert "quantify/copy_number.tsv" in manifest_paths
    report = (outdir / "run_report.md").read_text(encoding="utf-8")
    assert "## Output counts" in report
    assert "Family similarity rows:" in report
    assert "Repeat annotation summary: not generated" in report
    assert "Validation status: passed" in report


def test_full_toy_run_and_missing_assembly_skips(tmp_path: Path) -> None:
    toy = simulate_toy(tmp_path)
    full_outdir = tmp_path / "full"
    full = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--assembly",
        str(toy / "assembly.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(full_outdir),
        "--steps",
        "discover,quantify,locate,compare,probe,visualize,validate",
        "--kmer-backend",
        "rust",
    )
    assert full.returncode == 0, full.stderr
    assert [row["step"] for row in read_summary(full_outdir)] == [
        "discover",
        "quantify",
        "locate",
        "compare",
        "probe",
        "visualize",
        "validate",
    ]
    assert (full_outdir / "compare" / "assembly_vs_read_cn.tsv").is_file()
    assert (full_outdir / "visualize" / "catalogue_summary.svg").is_file()
    manifest_paths = {Path(row["file_path"]).relative_to(full_outdir).as_posix() for row in read_manifest(full_outdir)}
    assert "compare/assembly_vs_read_cn.tsv" in manifest_paths
    report = (full_outdir / "run_report.md").read_text(encoding="utf-8")
    assert "Compare rows:" in report
    assert "Compare status summary:" in report

    skip_outdir = tmp_path / "skip"
    skipped = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(skip_outdir),
        "--steps",
        "discover,locate,compare,probe,validate",
        "--kmer-backend",
        "rust",
    )
    assert skipped.returncode == 0, skipped.stderr
    rows = {row["step"]: row for row in read_summary(skip_outdir)}
    assert rows["locate"]["notes"] == "skipped_missing_assembly"
    assert rows["compare"]["notes"] == "skipped_missing_assembly"
    assert rows["probe"]["notes"] == "skipped_missing_assembly"
    skipped_manifest = read_manifest(skip_outdir)
    assert any(
        row["step"] == "locate" and row["notes"] == "skipped_missing_assembly"
        for row in skipped_manifest
    )
    assert any(
        row["step"] == "probe" and row["notes"] == "skipped_missing_assembly"
        for row in skipped_manifest
    )


def test_compare_step_reports_missing_arrays(tmp_path: Path) -> None:
    toy = simulate_toy(tmp_path)
    outdir = tmp_path / "missing_compare_inputs"
    result = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--assembly",
        str(toy / "assembly.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(outdir),
        "--steps",
        "discover,quantify,compare",
        "--kmer-backend",
        "rust",
    )

    assert result.returncode != 0
    rows = {row["step"]: row for row in read_summary(outdir)}
    assert rows["compare"]["notes"].startswith("missing_input:")
    assert "locate/arrays.bed" in rows["compare"]["notes"]
    report = (outdir / "run_report.md").read_text(encoding="utf-8")
    assert "compare: missing_input:" in report


def test_missing_genome_size_and_force_resume_behavior(tmp_path: Path) -> None:
    toy = simulate_toy(tmp_path)
    missing_outdir = tmp_path / "missing"
    missing = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--outdir",
        str(missing_outdir),
        "--steps",
        "discover,quantify",
        "--kmer-backend",
        "rust",
    )
    assert missing.returncode != 0
    assert read_summary(missing_outdir)[-1]["notes"] == "--genome-size is required when quantify is selected"

    outdir = tmp_path / "resume"
    base_args = (
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(outdir),
        "--steps",
        "discover,quantify,validate",
        "--kmer-backend",
        "rust",
    )
    first = run_cli(*base_args)
    assert first.returncode == 0, first.stderr
    blocked = run_cli(*base_args)
    assert blocked.returncode != 0
    assert read_summary(outdir)[0]["notes"] == "output_exists_use_force_or_resume"

    resumed = run_cli(*base_args, "--resume")
    assert resumed.returncode == 0, resumed.stderr
    resume_rows = {row["step"]: row for row in read_summary(outdir)}
    assert resume_rows["discover"]["notes"] == "skipped_validated_resume"
    assert resume_rows["quantify"]["notes"] == "skipped_validated_resume"

    forced = run_cli(*base_args, "--force")
    assert forced.returncode == 0, forced.stderr
    force_rows = read_summary(outdir)
    assert float(force_rows[0]["runtime_seconds"]) > 0


def test_run_passes_read_limits_to_discover_and_quantify(tmp_path: Path) -> None:
    toy = simulate_toy(tmp_path)
    outdir = tmp_path / "limited"
    result = run_cli(
        "run",
        "--reads",
        str(toy / "reads.fa"),
        "--genome-size",
        "1000000",
        "--outdir",
        str(outdir),
        "--steps",
        "discover,quantify",
        "--kmer-backend",
        "rust",
        "--max-reads",
        "20",
    )

    assert result.returncode == 0, result.stderr
    rows = read_summary(outdir)
    assert all("--max-reads 20" in row["command"] for row in rows)


def test_discover_has_no_known_repeat_input() -> None:
    result = run_cli("discover", "--help")
    assert result.returncode == 0
    assert "--known" not in result.stdout
    assert "known_repeats" not in result.stdout
