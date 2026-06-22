#!/usr/bin/env python3
"""Run a truth-aware read-level comparison of TandemX, TRF and TideHunter."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import shutil
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml

DNA = "ACGT"
RAW_FIELDS = [
    "benchmark_id", "dataset_id", "tool", "repetition", "runtime_seconds", "exit_status",
    "read_count", "read_length", "total_bases", "monomer_length", "error_rate",
    "positive_reads", "predicted_reads", "correct_reads", "false_positive_reads",
    "recall", "precision", "false_positive_rate", "monomer_length_mae_bp", "notes",
]


@dataclass(frozen=True)
class TruthRecord:
    read_id: str
    is_positive: bool
    monomer_length: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def mutate_substitutions(sequence: str, error_rate: float, rng: random.Random) -> str:
    bases = list(sequence)
    for index, base in enumerate(bases):
        if rng.random() < error_rate:
            bases[index] = rng.choice([candidate for candidate in DNA if candidate != base])
    return "".join(bases)


def write_dataset(dataset: dict[str, object], seed: int, outdir: Path) -> tuple[Path, Path]:
    rng = random.Random(seed)
    read_count = int(dataset["read_count"])
    read_length = int(dataset["read_length"])
    monomer_length = int(dataset["monomer_length"])
    positive_count = round(read_count * float(dataset["positive_fraction"]))
    monomer = random_sequence(rng, monomer_length)
    records: list[tuple[str, str]] = []
    truth: list[TruthRecord] = []
    for index in range(read_count):
        positive = index < positive_count
        read_id = f"{'repeat' if positive else 'background'}_{index + 1:06d}"
        if positive:
            sequence = (monomer * math.ceil(read_length / monomer_length))[:read_length]
            sequence = mutate_substitutions(sequence, float(dataset["error_rate"]), rng)
        else:
            sequence = random_sequence(rng, read_length)
        records.append((read_id, sequence))
        truth.append(TruthRecord(read_id, positive, monomer_length if positive else 0))
    rng.shuffle(records)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = outdir / "reads.fa"
    with fasta.open("w", encoding="utf-8") as handle:
        for read_id, sequence in records:
            handle.write(f">{read_id}\n{sequence}\n")
    truth_path = outdir / "truth.tsv"
    with truth_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["read_id", "truth_positive", "monomer_length_bp"])
        for row in truth:
            writer.writerow([row.read_id, str(row.is_positive).lower(), row.monomer_length])
    (outdir / "monomer.fa").write_text(f">truth_monomer\n{monomer}\n", encoding="utf-8")
    return fasta, truth_path


def read_truth(path: Path) -> dict[str, TruthRecord]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["read_id"]: TruthRecord(row["read_id"], row["truth_positive"] == "true", int(row["monomer_length_bp"]))
            for row in csv.DictReader(handle, delimiter="\t")
        }


def parse_tandemx(path: Path) -> dict[str, list[int]]:
    predictions: dict[str, list[int]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            predictions.setdefault(row["read_id"], []).append(int(row["period_bp"]))
    return predictions


def parse_tidehunter(path: Path) -> dict[str, list[int]]:
    predictions: dict[str, list[int]] = {}
    if not path.exists():
        return predictions
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 7 or fields[0].startswith("#"):
                continue
            try:
                period = int(fields[6])
            except ValueError:
                continue
            predictions.setdefault(fields[0], []).append(period)
    return predictions


def parse_trf(path: Path) -> dict[str, list[int]]:
    predictions: dict[str, list[int]] = {}
    current_read: str | None = None
    if not path.exists():
        return predictions
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("@"):
                current_read = stripped[1:].split()[0]
                continue
            fields = stripped.split()
            if current_read is None or len(fields) < 3:
                continue
            try:
                int(fields[0]); int(fields[1]); period = int(fields[2])
            except ValueError:
                continue
            predictions.setdefault(current_read, []).append(period)
    return predictions


def score_predictions(predictions: dict[str, list[int]], truth: dict[str, TruthRecord]) -> dict[str, float | int]:
    positive_ids = {read_id for read_id, record in truth.items() if record.is_positive}
    background_ids = set(truth) - positive_ids
    predicted_ids = set(predictions)
    correct_ids: set[str] = set()
    absolute_errors: list[int] = []
    for read_id in positive_ids:
        periods = predictions.get(read_id, [])
        if not periods:
            continue
        expected = truth[read_id].monomer_length
        error = min(abs(period - expected) for period in periods)
        absolute_errors.append(error)
        if error <= max(2, round(expected * 0.02)):
            correct_ids.add(read_id)
    false_positive_ids = predicted_ids & background_ids
    return {
        "positive_reads": len(positive_ids), "predicted_reads": len(predicted_ids),
        "correct_reads": len(correct_ids), "false_positive_reads": len(false_positive_ids),
        "recall": len(correct_ids) / len(positive_ids) if positive_ids else 0.0,
        "precision": len(correct_ids) / len(predicted_ids) if predicted_ids else 0.0,
        "false_positive_rate": len(false_positive_ids) / len(background_ids) if background_ids else 0.0,
        "monomer_length_mae_bp": statistics.fmean(absolute_errors) if absolute_errors else math.nan,
    }


def run_logged(command: list[str], stdout_path: Path, stderr_path: Path) -> tuple[int, float]:
    started = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    return completed.returncode, time.perf_counter() - started


def command_for_tool(tool: str, executable: Path, reads: Path, run_dir: Path, min_period: int, max_period: int, read_length: int) -> tuple[list[str], Path]:
    if tool == "tandemx":
        output = run_dir / "tandemx" / "candidate_reads.tsv"
        return ([str(executable), "discover", "--reads", str(reads), "--outdir", str(run_dir / "tandemx"),
                 "--min-period", str(min_period), "--max-period", str(max_period), "--min-support-reads", "1",
                 "--min-repeat-span", str(max(200, read_length // 2)), "--min-read-length", str(max(100, read_length // 2)),
                 "--kmer-size", "11", "--top-periods", "5", "--min-seed-occurrences", "2",
                 "--min-spacing-support", "2", "--max-pairs-per-kmer", "100", "--progress-every", "1000000",
                 "--kmer-backend", "rust"], output)
    if tool == "trf":
        output = run_dir / "trf.tsv"
        return ([str(executable), str(reads), "2", "7", "7", "80", "10", "50", str(max_period), "-ngs", "-h"], output)
    if tool == "tidehunter":
        output = run_dir / "tidehunter.tsv"
        return ([str(executable), "-t", "1", "-f", "2", "-p", str(min_period), "-P", str(max_period),
                 "-c", "2", "-o", str(output), str(reads)], output)
    raise ValueError(f"Unsupported tool: {tool}")


def parse_output(tool: str, path: Path) -> dict[str, list[int]]:
    return {"tandemx": parse_tandemx, "trf": parse_trf, "tidehunter": parse_tidehunter}[tool](path)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def aggregate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        groups.setdefault((str(row["dataset_id"]), str(row["tool"])), []).append(row)
    output: list[dict[str, object]] = []
    for (dataset_id, tool), group in sorted(groups.items()):
        reference = group[0]
        runtime = statistics.median(float(row["runtime_seconds"]) for row in group)
        output.append({
            "dataset_id": dataset_id, "tool": tool, "read_count": reference["read_count"],
            "total_bases": reference["total_bases"], "monomer_length": reference["monomer_length"],
            "error_rate": reference["error_rate"], "median_runtime_seconds": runtime,
            "median_mb_per_second": (float(reference["total_bases"]) / 1_000_000) / runtime,
            "recall": statistics.fmean(float(row["recall"]) for row in group),
            "precision": statistics.fmean(float(row["precision"]) for row in group),
            "false_positive_rate": statistics.fmean(float(row["false_positive_rate"]) for row in group),
            "monomer_length_mae_bp": statistics.fmean(float(row["monomer_length_mae_bp"]) for row in group),
            "successful_runs": sum(int(row["exit_status"]) == 0 for row in group),
        })
    return output


def main() -> int:
    args = parse_args(); config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    root = Path.cwd(); outdir = args.outdir.resolve(); outdir.mkdir(parents=True, exist_ok=True)
    tools = {name: (root / path).resolve() for name, path in config["tools"].items()}
    for name, executable in tools.items():
        if not executable.is_file(): raise SystemExit(f"Missing {name} executable: {executable}")
    min_period, max_period = (int(value) for value in config["period_range"])
    raw_rows: list[dict[str, object]] = []
    for dataset_index, dataset in enumerate(config["datasets"]):
        dataset_id = str(dataset["id"]); dataset_dir = outdir / "datasets" / dataset_id
        reads, truth_path = write_dataset(dataset, int(config["seed"]) + dataset_index, dataset_dir)
        truth = read_truth(truth_path)
        for tool, executable in tools.items():
            for repetition in range(1, int(config["repetitions"]) + 1):
                run_dir = outdir / "runs" / dataset_id / tool / f"rep_{repetition}"
                if run_dir.exists(): shutil.rmtree(run_dir)
                run_dir.mkdir(parents=True)
                command, output = command_for_tool(tool, executable, reads, run_dir, min_period, max_period, int(dataset["read_length"]))
                stdout_path = output if tool == "trf" else run_dir / "stdout.log"
                status, runtime = run_logged(command, stdout_path, run_dir / "stderr.log")
                predictions = parse_output(tool, output) if status == 0 and output.exists() else {}
                raw_rows.append({
                    "benchmark_id": config["benchmark_id"], "dataset_id": dataset_id, "tool": tool,
                    "repetition": repetition, "runtime_seconds": f"{runtime:.6f}", "exit_status": status,
                    "read_count": dataset["read_count"], "read_length": dataset["read_length"],
                    "total_bases": int(dataset["read_count"]) * int(dataset["read_length"]),
                    "monomer_length": dataset["monomer_length"], "error_rate": dataset["error_rate"],
                    **score_predictions(predictions, truth), "notes": "single_thread_same_fasta_substitution_errors_only",
                })
                print(f"{dataset_id}\t{tool}\trep={repetition}\t{runtime:.3f}s\texit={status}", flush=True)
    summary = aggregate(raw_rows)
    write_tsv(outdir / "raw_runs.tsv", raw_rows, RAW_FIELDS)
    write_tsv(outdir / "summary.tsv", summary, list(summary[0]) if summary else [])
    (outdir / "environment.json").write_text(json.dumps({
        "benchmark_id": config["benchmark_id"], "config": str(args.config.resolve()),
        "tools": {name: str(path) for name, path in tools.items()}, "threads": 1,
        "accuracy_tolerance": "max(2 bp, 2% of truth period)",
    }, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
