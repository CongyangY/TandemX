from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from tandemx.run_compare import compare_run_directories, load_run_artifacts


def write_discover_output(
    outdir: Path,
    reads: Path,
    *,
    max_reads: int,
    min_support_reads: int,
    min_period: int,
    max_period: int,
    top_periods: int,
    family_lengths: list[int],
    candidate_count: int,
) -> None:
    outdir.mkdir(parents=True)
    argv = [
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--max-reads",
        str(max_reads),
        "--min-support-reads",
        str(min_support_reads),
    ]
    lines = [
        'command: "tandemx discover"',
        'subcommand: "discover"',
        'version: "0.1.0"',
        'timestamp_utc: "2026-01-01T00:00:00+00:00"',
        'cwd: "/tmp"',
        "argv:",
        *(f'  - "{item}"' for item in argv),
        'python_version: "3.11.0"',
        'platform: "test"',
        'status: "discover_mvp_completed"',
        "parameters:",
        "  chunk_size: 1000",
        '  command: "discover"',
        '  kmer_backend: "python"',
        "  kmer_size: 11",
        f"  max_monomer_len: {max_period}",
        "  max_pairs_per_kmer: 100",
        "  max_read_bases: null",
        f"  max_reads: {max_reads}",
        f"  min_monomer_len: {min_period}",
        "  min_read_length: 1",
        "  min_repeat_span: 100",
        "  min_seed_occurrences: 2",
        "  min_spacing_support: 2",
        f"  min_support_reads: {min_support_reads}",
        f'  outdir: "{outdir}"',
        "  progress_every: 1000",
        f'  reads: "{reads}"',
        "  sample_rate: 1.0",
        "  seed: 1",
        f"  top_periods: {top_periods}",
    ]
    (outdir / "run_config.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")

    with (outdir / "candidate_reads.tsv").open("wt", encoding="utf-8") as handle:
        handle.write(
            "read_id\tcandidate_id\tread_start\tread_end\tstrand\tperiod_bp\trepeat_span_bp\t"
            "unit_count\tscore\tlow_complexity_flag\tconfidence\twarning\n"
        )
        for index in range(candidate_count):
            handle.write(f"read_{index}\tTXC{index:06d}\t0\t100\t.\t20\t100\t5\t1.0\tfalse\tmedium\t\n")

    with (outdir / "families.tsv").open("wt", encoding="utf-8") as handle:
        handle.write(
            "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\t"
            "support_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\n"
        )
        for index, length in enumerate(family_lengths, start=1):
            handle.write(
                f"TXF{index:06d}\tTXM{index:06d}\t{length}\tmd5-{index}\t0.5\t"
                f"{min_support_reads}\t{length * 10}\t0.9\tfalse\tmedium\t\n"
            )

    with (outdir / "monomers.fa").open("wt", encoding="utf-8") as handle:
        for index, length in enumerate(family_lengths, start=1):
            handle.write(f">family_id=TXF{index:06d};monomer_id=TXM{index:06d};length_bp={length};confidence=medium\n")
            handle.write("A" * length + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_compare_run_directories_reports_parameter_and_result_differences(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read_1\nACGTACGT\n", encoding="utf-8")
    run_a = tmp_path / "discover_benchmark" / "discover_100"
    run_b = tmp_path / "pipeline_run" / "discover"
    write_discover_output(
        run_a,
        reads,
        max_reads=100,
        min_support_reads=1,
        min_period=50,
        max_period=1000,
        top_periods=3,
        family_lengths=[155, 354, 706],
        candidate_count=5,
    )
    write_discover_output(
        run_b,
        reads,
        max_reads=100,
        min_support_reads=5,
        min_period=20,
        max_period=2000,
        top_periods=5,
        family_lengths=[354],
        candidate_count=5,
    )
    (tmp_path / "pipeline_run" / "pipeline_summary.tsv").write_text(
        "run_id\tstep\truntime_seconds\nrun-1\tdiscover\t1.0\n",
        encoding="utf-8",
    )

    rows = compare_run_directories(tmp_path / "discover_benchmark", tmp_path / "pipeline_run", tmp_path / "compare")
    rows_by_item = {row["item"]: row for row in rows}

    assert rows_by_item["reads"]["same"] == "true"
    assert rows_by_item["max_reads"]["same"] == "true"
    assert rows_by_item["min_support_reads"]["same"] == "false"
    assert rows_by_item["candidate_reads_count"]["run_a_value"] == "5"
    assert rows_by_item["family_count"]["run_a_value"] == "3"
    assert rows_by_item["family_count"]["run_b_value"] == "1"
    assert rows_by_item["monomer_lengths_bp"]["run_a_value"] == "155,354,706"
    assert rows_by_item["directly_comparable"]["run_a_value"] == "false"
    assert "min_support_reads" in rows_by_item["reason_not_directly_comparable"]["notes"]
    assert (tmp_path / "compare" / "compare_runs.tsv").is_file()
    assert (tmp_path / "compare" / "compare_runs.md").is_file()


def test_load_run_artifacts_accepts_standalone_discover_directory(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read_1\nACGTACGT\n", encoding="utf-8")
    discover = tmp_path / "discover"
    write_discover_output(
        discover,
        reads,
        max_reads=10,
        min_support_reads=1,
        min_period=20,
        max_period=100,
        top_periods=2,
        family_lengths=[20],
        candidate_count=1,
    )

    artifacts = load_run_artifacts(discover, "single")

    assert artifacts.discover_dir == discover
    assert artifacts.candidate_count == 1
    assert artifacts.monomer_lengths == [20]


def test_compare_tandemx_runs_script_writes_reports(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read_1\nACGTACGT\n", encoding="utf-8")
    run_a = tmp_path / "a" / "discover"
    run_b = tmp_path / "b" / "discover"
    write_discover_output(
        run_a,
        reads,
        max_reads=10,
        min_support_reads=1,
        min_period=20,
        max_period=100,
        top_periods=2,
        family_lengths=[20],
        candidate_count=1,
    )
    write_discover_output(
        run_b,
        reads,
        max_reads=10,
        min_support_reads=1,
        min_period=20,
        max_period=100,
        top_periods=2,
        family_lengths=[20],
        candidate_count=1,
    )
    outdir = tmp_path / "script_compare"

    result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/compare_tandemx_runs.py",
            "--run-a",
            str(tmp_path / "a"),
            "--run-b",
            str(tmp_path / "b"),
            "--outdir",
            str(outdir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert "directly_comparable=true" in result.stdout
    assert read_tsv(outdir / "compare_runs.tsv")
    assert (outdir / "compare_runs.md").read_text(encoding="utf-8").startswith("# TandemX run comparison")
