#!/usr/bin/env python3
"""Run synthetic TandemX benchmark scales from a YAML configuration."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError as exc:  # pragma: no cover - exercised only outside the dev env
    raise SystemExit("PyYAML is required to read benchmark configs. Use the tandemx-dev environment.") from exc


SUMMARY_FIELDS = [
    "benchmark_id",
    "scale",
    "seed",
    "read_count",
    "read_length",
    "total_read_bp",
    "monomer_lengths",
    "command",
    "runtime_seconds",
    "exit_status",
    "output_validated",
    "recovered_family_count",
    "processed_reads",
    "processed_bases",
    "candidate_reads",
    "candidates_per_mb",
    "reads_per_second",
    "mb_per_second",
    "peak_memory_mb",
    "algorithm_mode",
    "notes",
]

ACCURACY_FIELDS = [
    "benchmark_id",
    "expected_monomer_length",
    "recovered_closest_length",
    "length_error_bp",
    "expected_read_copy_bp",
    "estimated_read_copy_bp",
    "copy_number_relative_error",
    "locate_status",
    "notes",
]


@dataclass(frozen=True)
class CommandResult:
    name: str
    runtime_seconds: float
    exit_status: int
    notes: str


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run TandemX synthetic scale benchmarks.")
    parser.add_argument("--config", required=True, type=Path, help="Path to synthetic_scale.yaml.")
    parser.add_argument("--scale", required=True, help="Benchmark scale to run, for example tiny or small.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for benchmark results.")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    benchmarks = config.get("benchmarks", {})
    if args.scale not in benchmarks:
        available = ", ".join(sorted(benchmarks))
        raise SystemExit(f"Unknown scale '{args.scale}'. Available scales: {available}")

    scale_config = benchmarks[args.scale]
    return run_benchmark(args.scale, scale_config, args.outdir)


def load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"Config file does not exist: {path}")
    with path.open("rt", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise SystemExit(f"Config file must contain a mapping: {path}")
    return data


def run_benchmark(scale: str, config: dict[str, Any], outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    scale_dir = outdir / scale
    logs_dir = scale_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    benchmark_id = str(config["benchmark_id"])
    read_count = int(config["read_count"])
    read_length = int(config["read_length"])
    total_read_bp = int(config.get("total_read_bp", read_count * read_length))
    monomer_lengths = [int(value) for value in config["monomer_lengths"]]
    copy_numbers = [int(value) for value in config["copy_numbers"]]
    seed = int(config["seed"])
    source_length = int(config["background_length"]) + sum(
        length * copies for length, copies in zip(monomer_lengths, copy_numbers)
    ) + 100 * len(copy_numbers)
    haploid_depth = total_read_bp / source_length

    commands = build_commands(config, scale_dir, source_length, haploid_depth)
    results: list[CommandResult] = []
    exit_code = 0
    for name, command in commands:
        result = run_command(name, command, logs_dir)
        results.append(result)
        if result.exit_status != 0:
            exit_code = result.exit_status
            break

    output_validated = any(result.name == "validate" and result.exit_status == 0 for result in results)
    recovered_family_count = count_recovered_families(scale_dir / "discover" / "families.tsv")
    discover_metrics = parse_discover_metrics(scale_dir / "discover" / "run.log")
    summary_rows = [
        summary_row(
            config=config,
            scale=scale,
            command_result=result,
            total_read_bp=total_read_bp,
            output_validated=output_validated if result.name == "validate" else False,
            recovered_family_count=recovered_family_count,
            discover_metrics=discover_metrics if result.name == "discover" else {},
        )
        for result in results
    ]
    write_tsv(outdir / "benchmark_summary.tsv", SUMMARY_FIELDS, summary_rows)
    accuracy_rows = build_accuracy_rows(benchmark_id, scale_dir)
    write_tsv(outdir / "accuracy_summary.tsv", ACCURACY_FIELDS, accuracy_rows)

    missing = missing_expected_outputs(scale_dir, config.get("expected_output_files", []))
    if missing and exit_code == 0:
        exit_code = 1
        append_summary_note(outdir / "benchmark_summary.tsv", f"missing_expected_outputs={';'.join(missing)}")
    return exit_code


def build_commands(
    config: dict[str, Any],
    scale_dir: Path,
    source_length: int,
    haploid_depth: float,
) -> list[tuple[str, list[str]]]:
    simulated = scale_dir / "simulated"
    discover = scale_dir / "discover"
    quantify = scale_dir / "quantify"
    locate = scale_dir / "locate"
    probe = scale_dir / "probe"
    monomer_lengths = ",".join(str(value) for value in config["monomer_lengths"])
    copy_numbers = ",".join(str(value) for value in config["copy_numbers"])
    discover_config = config["discover"]
    locate_config = config["locate"]
    probe_config = config["probe"]
    cli = [sys.executable, "-m", "tandemx.cli"]

    return [
        (
            "simulate",
            [
                *cli,
                "simulate",
                "toy",
                "--outdir",
                str(simulated),
                "--seed",
                str(config["seed"]),
                "--num-reads",
                str(config["read_count"]),
                "--read-length",
                str(config["read_length"]),
                "--background-length",
                str(config["background_length"]),
                "--monomer-lengths",
                monomer_lengths,
                "--copies",
                copy_numbers,
                "--error-rate",
                str(config["error_rate"]),
            ],
        ),
        (
            "discover",
            [
                *cli,
                "discover",
                "--reads",
                str(simulated / "reads.fa"),
                "--outdir",
                str(discover),
                "--min-period",
                str(discover_config["min_period"]),
                "--max-period",
                str(discover_config["max_period"]),
                "--min-support-reads",
                str(discover_config["min_support_reads"]),
                "--min-repeat-span",
                str(discover_config["min_repeat_span"]),
                "--min-read-length",
                str(discover_config["min_read_length"]),
                "--kmer-size",
                str(discover_config["kmer_size"]),
                "--top-periods",
                str(discover_config["top_periods"]),
                "--min-seed-occurrences",
                str(discover_config["min_seed_occurrences"]),
                "--min-spacing-support",
                str(discover_config["min_spacing_support"]),
                "--max-pairs-per-kmer",
                str(discover_config["max_pairs_per_kmer"]),
                "--max-reads",
                str(config["read_count"]),
                "--progress-every",
                str(discover_config["progress_every"]),
            ],
        ),
        (
            "quantify",
            [
                *cli,
                "quantify",
                "--reads",
                str(simulated / "reads.fa"),
                "--catalog",
                str(discover / "monomers.fa"),
                "--genome-size",
                str(source_length),
                "--haploid-depth",
                f"{haploid_depth:.6f}",
                "--outdir",
                str(quantify),
            ],
        ),
        (
            "locate",
            [
                *cli,
                "locate",
                "--assembly",
                str(simulated / "assembly.fa"),
                "--catalog",
                str(discover / "monomers.fa"),
                "--copy-number",
                str(quantify / "copy_number.tsv"),
                "--window-size",
                str(locate_config["window_size"]),
                "--step-size",
                str(locate_config["step_size"]),
                "--k",
                str(locate_config["k"]),
                "--outdir",
                str(locate),
            ],
        ),
        (
            "probe",
            [
                *cli,
                "probe",
                "--catalog",
                str(discover / "monomers.fa"),
                "--assembly",
                str(simulated / "assembly.fa"),
                "--copy-number",
                str(quantify / "copy_number.tsv"),
                "--arrays",
                str(locate / "arrays.bed"),
                "--min-len",
                str(probe_config["min_len"]),
                "--max-len",
                str(probe_config["max_len"]),
                "--outdir",
                str(probe),
            ],
        ),
        ("validate", [*cli, "validate", "--project", str(scale_dir)]),
    ]


def run_command(name: str, command: list[str], logs_dir: Path) -> CommandResult:
    start = time.perf_counter()
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    runtime = time.perf_counter() - start
    (logs_dir / f"{name}.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (logs_dir / f"{name}.stderr.log").write_text(completed.stderr, encoding="utf-8")
    (logs_dir / f"{name}.command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")
    notes = "peak_memory_not_recorded_python_stdlib"
    if completed.stderr.strip():
        notes = f"{notes};stderr={completed.stderr.strip().splitlines()[-1]}"
    return CommandResult(name=name, runtime_seconds=runtime, exit_status=completed.returncode, notes=notes)


def summary_row(
    config: dict[str, Any],
    scale: str,
    command_result: CommandResult,
    total_read_bp: int,
    output_validated: bool,
    recovered_family_count: int,
    discover_metrics: dict[str, str],
) -> dict[str, str]:
    return {
        "benchmark_id": str(config["benchmark_id"]),
        "scale": scale,
        "seed": str(config["seed"]),
        "read_count": str(config["read_count"]),
        "read_length": str(config["read_length"]),
        "total_read_bp": str(total_read_bp),
        "monomer_lengths": ",".join(str(value) for value in config["monomer_lengths"]),
        "command": command_result.name,
        "runtime_seconds": f"{command_result.runtime_seconds:.4f}",
        "exit_status": str(command_result.exit_status),
        "output_validated": str(output_validated).lower(),
        "recovered_family_count": str(recovered_family_count),
        "processed_reads": discover_metrics.get("processed_reads", ""),
        "processed_bases": discover_metrics.get("processed_bases", ""),
        "candidate_reads": discover_metrics.get("candidate_reads", ""),
        "candidates_per_mb": discover_metrics.get("candidates_per_mb", ""),
        "reads_per_second": discover_metrics.get("reads_per_second", ""),
        "mb_per_second": discover_metrics.get("mb_per_second", ""),
        "peak_memory_mb": "",
        "algorithm_mode": "spacing_prefilter" if command_result.name == "discover" else "",
        "notes": command_result.notes,
    }


def parse_discover_metrics(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    progress_lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if "progress processed_reads=" in line
    ]
    if not progress_lines:
        return {}
    values = dict(re.findall(r"([a-z_]+)=([^ ]+)", progress_lines[-1]))
    processed_bases = float(values.get("processed_bases", "0"))
    candidate_reads = float(values.get("candidate_reads", "0"))
    values["candidates_per_mb"] = (
        f"{candidate_reads / (processed_bases / 1_000_000):.6f}"
        if processed_bases > 0
        else "0.000000"
    )
    return values


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def append_summary_note(path: Path, note: str) -> None:
    with path.open("at", encoding="utf-8") as handle:
        handle.write(f"# {note}\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_recovered_families(path: Path) -> int:
    return len(read_tsv(path))


def build_accuracy_rows(benchmark_id: str, scale_dir: Path) -> list[dict[str, str]]:
    truth_rows = read_tsv(scale_dir / "simulated" / "truth_copy_number.tsv")
    family_rows = read_tsv(scale_dir / "discover" / "families.tsv")
    copy_rows = read_tsv(scale_dir / "quantify" / "copy_number.tsv")
    comparison_rows = read_tsv(scale_dir / "locate" / "assembly_vs_read_cn.tsv")
    copy_by_family = {row["family_id"]: row for row in copy_rows if "family_id" in row}
    status_by_family = {row["family_id"]: row.get("status", "") for row in comparison_rows if "family_id" in row}
    recovered = [
        (int(row["monomer_length_bp"]), row["family_id"])
        for row in family_rows
        if row.get("monomer_length_bp", "").isdigit()
    ]

    rows = []
    for truth in truth_rows:
        expected_length = int(truth["monomer_length_bp"])
        expected_bp = float(truth["read_repeat_bp"])
        closest = min(recovered, key=lambda item: abs(item[0] - expected_length), default=None)
        if closest is None:
            rows.append(
                {
                    "benchmark_id": benchmark_id,
                    "expected_monomer_length": str(expected_length),
                    "recovered_closest_length": "",
                    "length_error_bp": "",
                    "expected_read_copy_bp": f"{expected_bp:.4f}",
                    "estimated_read_copy_bp": "",
                    "copy_number_relative_error": "",
                    "locate_status": "",
                    "notes": "no_recovered_family",
                }
            )
            continue
        recovered_length, family_id = closest
        estimated_bp = float(copy_by_family.get(family_id, {}).get("estimated_bp", "0") or 0)
        relative_error = abs(estimated_bp - expected_bp) / expected_bp if expected_bp > 0 else 0.0
        rows.append(
            {
                "benchmark_id": benchmark_id,
                "expected_monomer_length": str(expected_length),
                "recovered_closest_length": str(recovered_length),
                "length_error_bp": str(abs(recovered_length - expected_length)),
                "expected_read_copy_bp": f"{expected_bp:.4f}",
                "estimated_read_copy_bp": f"{estimated_bp:.4f}",
                "copy_number_relative_error": f"{relative_error:.6f}",
                "locate_status": status_by_family.get(family_id, ""),
                "notes": "truth_used_for_benchmark_evaluation_only",
            }
        )
    return rows


def missing_expected_outputs(scale_dir: Path, expected_files: list[str]) -> list[str]:
    return [path for path in expected_files if not (scale_dir / path).is_file()]


if __name__ == "__main__":
    raise SystemExit(main())
