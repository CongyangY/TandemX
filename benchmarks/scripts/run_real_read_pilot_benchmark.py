#!/usr/bin/env python3
"""Benchmark TandemX discover and validate on bounded real-read subsets."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path


FIELDS = [
    "input_file",
    "max_reads",
    "processed_reads",
    "processed_bases",
    "runtime_seconds",
    "reads_per_second",
    "mb_per_second",
    "candidate_reads",
    "candidate_rate",
    "recovered_family_count",
    "peak_memory_mb",
    "output_validated",
    "exit_status",
    "command",
    "notes",
]

PROGRESS_PATTERN = re.compile(
    r"processed_reads=(?P<reads>\d+).*?processed_bases=(?P<bases>\d+).*?"
    r"candidate_reads=(?P<candidates>\d+).*?reads_per_second=(?P<reads_s>[\d.]+).*?"
    r"mb_per_second=(?P<mb_s>[\d.]+)"
)


def parse_positive_int_list(value: str) -> list[int]:
    values = [int(part.strip()) for part in value.split(",") if part.strip()]
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("--max-reads requires positive comma-separated integers")
    return values


def parse_progress(log_path: Path) -> dict[str, str]:
    matches = PROGRESS_PATTERN.findall(log_path.read_text(encoding="utf-8"))
    if not matches:
        return {}
    reads, bases, candidates, reads_s, mb_s = matches[-1]
    return {
        "processed_reads": reads,
        "processed_bases": bases,
        "candidate_reads": candidates,
        "reads_per_second": reads_s,
        "mb_per_second": mb_s,
    }


def count_families(path: Path) -> int:
    if not path.exists():
        return 0
    return max(0, len(path.read_text(encoding="utf-8").splitlines()) - 1)


def run_logged(command: list[str], log_path: Path) -> tuple[int, float]:
    started = time.perf_counter()
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    runtime = time.perf_counter() - started
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        result.stdout + ("\n" if result.stdout and result.stderr else "") + result.stderr,
        encoding="utf-8",
    )
    return result.returncode, runtime


def run_scale(args: argparse.Namespace, max_reads: int) -> dict[str, str | int]:
    scale_dir = args.outdir / f"discover_{max_reads}"
    command = [
        sys.executable,
        "-m",
        "tandemx.cli",
        "discover",
        "--reads",
        str(args.reads),
        "--outdir",
        str(scale_dir),
        "--max-reads",
        str(max_reads),
        "--min-period",
        str(args.min_period),
        "--max-period",
        str(args.max_period),
        "--top-periods",
        str(args.top_periods),
        "--min-spacing-support",
        str(args.min_spacing_support),
        "--min-support-reads",
        str(args.min_support_reads),
        "--progress-every",
        str(args.progress_every),
    ]
    if args.max_read_bases is not None:
        command.extend(["--max-read-bases", str(args.max_read_bases)])

    exit_status, runtime = run_logged(command, args.outdir / "logs" / f"discover_{max_reads}.log")
    metrics = parse_progress(scale_dir / "run.log") if (scale_dir / "run.log").exists() else {}

    validated = False
    if exit_status == 0:
        validate_command = [
            sys.executable,
            "-m",
            "tandemx.cli",
            "validate",
            "--project",
            str(scale_dir),
        ]
        validate_status, _ = run_logged(
            validate_command,
            args.outdir / "logs" / f"validate_{max_reads}.log",
        )
        validated = validate_status == 0

    processed_reads = int(metrics.get("processed_reads", "0"))
    candidate_reads = int(metrics.get("candidate_reads", "0"))
    return {
        "input_file": str(args.reads),
        "max_reads": max_reads,
        "processed_reads": processed_reads,
        "processed_bases": metrics.get("processed_bases", "0"),
        "runtime_seconds": f"{runtime:.4f}",
        "reads_per_second": metrics.get("reads_per_second", "0"),
        "mb_per_second": metrics.get("mb_per_second", "0"),
        "candidate_reads": candidate_reads,
        "candidate_rate": f"{candidate_reads / processed_reads:.6f}" if processed_reads else "0",
        "recovered_family_count": count_families(scale_dir / "families.tsv"),
        "peak_memory_mb": "NA",
        "output_validated": str(validated).lower(),
        "exit_status": exit_status,
        "command": shlex.join(command),
        "notes": "peak_memory_not_portably_available",
    }


def write_summary(path: Path, rows: list[dict[str, str | int]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--max-reads", required=True, type=parse_positive_int_list)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--max-read-bases", type=int)
    parser.add_argument("--min-period", type=int, default=50)
    parser.add_argument("--max-period", type=int, default=1000)
    parser.add_argument("--top-periods", type=int, default=3)
    parser.add_argument("--min-spacing-support", type=int, default=2)
    parser.add_argument("--min-support-reads", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args()

    if not args.reads.is_file():
        parser.error(f"reads file does not exist: {args.reads}")
    if args.max_read_bases is not None and args.max_read_bases <= 0:
        parser.error("--max-read-bases must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []
    summary_path = args.outdir / "tmpfq_benchmark_summary.tsv"
    for max_reads in args.max_reads:
        rows.append(run_scale(args, max_reads))
        write_summary(summary_path, rows)
    return 0 if all(int(row["exit_status"]) == 0 for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
