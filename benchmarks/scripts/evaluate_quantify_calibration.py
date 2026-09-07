"""Run a paired development ablation of public quantify depth calibration."""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from tandemx.discover.mvp import read_fasta
from tandemx.quantify.mvp import (
    family_kmer_membership,
    monomer_kmer_counts,
    read_monomer_fasta,
)
from tandemx.utils.kmers import canonical_kmer, is_low_complexity_kmer


METHODS = (
    "baseline_total_bases",
    "oracle_error_survival",
    "empirical_controls",
    "empirical_controls_plus_oracle_error",
)


def select_single_copy_controls(
    genome_path: Path,
    truth_path: Path,
    catalogue_path: Path,
    *,
    k: int,
    desired: int,
    stride: int,
    candidate_multiplier: int,
) -> tuple[list[str], dict[str, object]]:
    """Select and whole-genome verify canonical controls outside planted arrays."""
    if min(k, desired, stride, candidate_multiplier) < 1:
        raise ValueError("Control-selection sizes must be positive")
    records = list(read_fasta(genome_path))
    if len(records) != 1:
        raise ValueError("Calibration genomes must contain exactly one sequence")
    genome = records[0].sequence.upper()
    if set(genome) - set("ACGT"):
        raise ValueError("Calibration genome must contain only ACGT")
    truth = read_table(truth_path, {"start", "end"})
    intervals = sorted((int(row["start"]), int(row["end"])) for row in truth)
    catalogue = list(read_monomer_fasta(catalogue_path))
    membership = family_kmer_membership(catalogue, k)
    repeat_targets = {
        word
        for monomer in catalogue
        for word in monomer_kmer_counts(monomer.sequence, k)
        if len(membership[word]) == 1
    }
    candidates: list[str] = []
    candidate_set: set[str] = set()
    interval_index = 0
    limit = desired * candidate_multiplier
    for start in range(0, len(genome) - k + 1, stride):
        while interval_index < len(intervals) and intervals[interval_index][1] <= start:
            interval_index += 1
        if interval_index < len(intervals):
            array_start, array_end = intervals[interval_index]
            if start < array_end and start + k > array_start:
                continue
        word = genome[start : start + k]
        word = canonical_kmer(word)
        if (
            word in candidate_set
            or word in repeat_targets
            or is_low_complexity_kmer(word)
        ):
            continue
        candidates.append(word)
        candidate_set.add(word)
        if len(candidates) >= limit:
            break
    if len(candidates) < desired:
        raise ValueError("Insufficient candidate background k-mers")
    counts: Counter[str] = Counter()
    for start in range(len(genome) - k + 1):
        canonical = canonical_kmer(genome[start : start + k])
        if canonical in candidate_set:
            counts[canonical] += 1
    selected = [word for word in candidates if counts[word] == 1][:desired]
    if len(selected) != desired:
        raise ValueError(
            f"Only {len(selected)} of {desired} requested controls are genome-unique"
        )
    receipt = {
        "k": k,
        "requested_controls": desired,
        "selected_controls": len(selected),
        "candidate_controls": len(candidates),
        "selection_stride_bp": stride,
        "excluded_planted_array_count": len(intervals),
        "genome_sha256": digest_file(genome_path),
        "truth_sha256": digest_file(truth_path),
        "catalogue_sha256": digest_file(catalogue_path),
        "criterion": "canonical_whole_genome_occurrence_equals_1;outside_planted_array_intervals;not_repeat_diagnostic;not_low_complexity",
        "warning": "simulation_truth_assisted_control_selection_not_available_for_blind_real_data",
    }
    return selected, receipt


def summarize(metrics: list[dict[str, object]], executions: list[dict[str, object]]) -> list[dict[str, object]]:
    by_execution = {
        (row["seed"], row["condition_id"], row["method"]): row for row in executions
    }
    groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in metrics:
        groups[(row["method"], row["coverage"], row["error_model"])].append(row)
    output: list[dict[str, object]] = []
    for (method, coverage, error_model), rows in sorted(groups.items()):
        run_rows = [
            by_execution[key]
            for key in {
                (row["seed"], row["condition_id"], row["method"]) for row in rows
            }
        ]
        interval_rows = [row for row in rows if row["interval_contains_truth"] is not None]
        output.append(
            {
                "method": method,
                "coverage": coverage,
                "error_model": error_model,
                "independent_genomes": len({row["seed"] for row in rows}),
                "family_conditions": len(rows),
                "mean_signed_relative_error": statistics.mean(float(row["signed_relative_error"]) for row in rows),
                "mean_absolute_relative_error": statistics.mean(float(row["absolute_relative_error"]) for row in rows),
                "mean_estimator_minus_sampling_oracle": statistics.mean(float(row["estimator_minus_sampling_oracle"]) for row in rows),
                "diagnostic_spread_truth_coverage": statistics.mean(bool(row["interval_contains_truth"]) for row in interval_rows) if interval_rows else None,
                "mean_diagnostic_spread_relative_width": statistics.mean(float(row["interval_relative_width"]) for row in interval_rows) if interval_rows else None,
                "median_runtime_seconds": statistics.median(float(row["runtime_seconds"]) for row in run_rows),
                "median_peak_rss_mib": statistics.median(float(row["peak_rss_mib"]) for row in run_rows),
                "warning": "dependent_families_within_three_genomes;diagnostic_spread_not_sampling_CI;single_resource_execution_per_condition",
            }
        )
    return output


def worker(config_path: Path, datasets: list[Path], outdir: Path) -> None:
    config = yaml.safe_load(config_path.read_text())
    if config.get("split") != "development" or tuple(config.get("methods", ())) != METHODS:
        raise ValueError("Require the frozen development method matrix")
    seeds = list(config["seeds"])
    receipts = [json.loads((path / "generation_receipt.json").read_text()) for path in datasets]
    if [receipt["seed"] for receipt in receipts] != seeds:
        raise ValueError("Dataset order/seeds differ from the frozen config")
    if any(
        not receipt["complete"] or receipt["heldout_used"] or receipt["split"] != "development"
        for receipt in receipts
    ):
        raise ValueError("Calibration requires complete development datasets")
    if any(
        len(receipt["conditions_completed"]) != int(config["condition_count_per_seed"])
        for receipt in receipts
    ):
        raise ValueError("Condition count differs from the frozen config")

    executions: list[dict[str, object]] = []
    metrics: list[dict[str, object]] = []
    for dataset, generation in zip(datasets, receipts):
        genome_dir = dataset / "genome"
        genome_receipt = json.loads((genome_dir / "manifest.json").read_text())
        for name in ("genome.fa", "catalogue.fa", "truth_copy_number.tsv"):
            if digest_file(genome_dir / name) != genome_receipt["files"][name]:
                raise ValueError(f"Genome input hash differs: {name}")
        controls, control_receipt = select_single_copy_controls(
            genome_dir / "genome.fa",
            genome_dir / "truth_copy_number.tsv",
            genome_dir / "catalogue.fa",
            k=int(config["k"]),
            desired=int(config["control_kmer_count"]),
            stride=int(config["control_stride_bp"]),
            candidate_multiplier=int(config["control_candidate_multiplier"]),
        )
        control_dir = outdir / "controls" / f"s{generation['seed']}"
        control_dir.mkdir(parents=True)
        controls_path = control_dir / "single_copy_kmers.tsv"
        controls_path.write_text(
            "kmer\texpected_copy_number\n"
            + "".join(f"{word}\t1\n" for word in controls)
        )
        control_receipt["controls_sha256"] = digest_file(controls_path)
        (control_dir / "receipt.json").write_text(json.dumps(control_receipt, indent=2) + "\n")
        truth = {
            row["family_id"]: row
            for row in read_table(genome_dir / "truth_copy_number.tsv")
        }
        if len(truth) != generation["family_count"]:
            raise ValueError("Truth family count differs from generation receipt")
        for condition in generation["conditions_completed"]:
            read_dir = dataset / "reads" / condition["condition_id"]
            sampling_path = read_dir / "manifest.json"
            if digest_file(sampling_path) != condition["manifest_sha256"]:
                raise ValueError("Read manifest hash differs from generation receipt")
            sampling = json.loads(sampling_path.read_text())
            reads_path = read_dir / "reads.fa"
            if digest_file(reads_path) != sampling["files"]["reads.fa"]:
                raise ValueError("Read hash differs from sampling receipt")
            total_error = sum(
                float(sampling[name])
                for name in ("substitution_rate", "insertion_rate", "deletion_rate")
            )
            for method in METHODS:
                run_dir = (
                    outdir
                    / "runs"
                    / f"s{generation['seed']}"
                    / condition["condition_id"]
                    / method
                )
                run_dir.mkdir(parents=True)
                command = [
                    sys.executable,
                    "-m",
                    "tandemx.cli",
                    "quantify",
                    "--reads",
                    str(reads_path),
                    "--catalog",
                    str(genome_dir / "catalogue.fa"),
                    "--genome-size",
                    str(genome_receipt["genome_bp"]),
                    "--k",
                    str(config["k"]),
                    "--kmer-backend",
                    "rust",
                    "--no-progress",
                    "--outdir",
                    str(run_dir / "output"),
                ]
                if "controls" in method:
                    command.extend(["--single-copy-kmers", str(controls_path)])
                if "error" in method:
                    command.extend(["--read-error-rate", str(total_error)])
                (run_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                measured = run_process(
                    command,
                    run_dir / "stdout.log",
                    run_dir / "stderr.log",
                    float(config["timeout_seconds"]),
                )
                execution = {
                    "seed": generation["seed"],
                    "condition_id": condition["condition_id"],
                    "coverage": condition["coverage"],
                    "error_model": condition["label"],
                    "method": method,
                    "total_error_rate_supplied": total_error if "error" in method else None,
                    **measured,
                }
                execution["status"] = (
                    "ok" if execution["exit_code"] == 0 and not execution["timed_out"] else "failed"
                )
                executions.append(execution)
                write_table(outdir / "executions.tsv", executions, list(executions[0]))
                (run_dir / "receipt.json").write_text(json.dumps(execution, indent=2) + "\n")
                if execution["status"] != "ok":
                    continue
                estimates = {
                    row["family_id"]: row
                    for row in read_table(run_dir / "output" / "copy_number.tsv")
                }
                if estimates.keys() != truth.keys():
                    raise ValueError("Quantify output family set differs from truth")
                for family_id, estimate in estimates.items():
                    truth_row = truth[family_id]
                    copies = float(truth_row["copies"])
                    value = float(estimate["estimated_copy_number"])
                    low = float(estimate["copy_number_interval_low"])
                    high = float(estimate["copy_number_interval_high"])
                    oracle = (
                        float(sampling["sampled_repeat_bp"].get(family_id, 0))
                        / float(truth_row["period"])
                        / float(sampling["actual_source_coverage"])
                    )
                    metrics.append(
                        {
                            "seed": generation["seed"],
                            "condition_id": condition["condition_id"],
                            "coverage": condition["coverage"],
                            "error_model": condition["label"],
                            "method": method,
                            "family_id": family_id,
                            "period": truth_row["period"],
                            "unit_substitution_rate": truth_row["unit_substitution_rate"],
                            "truth_copies": copies,
                            "estimated_copy_number": value,
                            "signed_relative_error": (value - copies) / copies,
                            "absolute_relative_error": abs(value - copies) / copies,
                            "sampling_oracle_copy_estimate": oracle,
                            "estimator_minus_sampling_oracle": (value - oracle) / copies,
                            "interval_low": low,
                            "interval_high": high,
                            "interval_contains_truth": low <= copies <= high,
                            "interval_relative_width": (high - low) / copies,
                            "normalization_method": estimate["normalization_method"],
                            "kmer_survival_probability": estimate["kmer_survival_probability"],
                            "control_mean_depth": estimate.get("single_copy_control_mean_depth", "NA"),
                            "control_zero_fraction": estimate.get("single_copy_control_zero_fraction", "NA"),
                        }
                    )
                write_table(outdir / "metrics.tsv", metrics, list(metrics[0]))
    summaries = summarize(metrics, executions)
    write_table(outdir / "summary.tsv", summaries, list(summaries[0]))
    validation = {
        "complete": len(executions)
        == len(datasets) * int(config["condition_count_per_seed"]) * len(METHODS)
        and all(row["status"] == "ok" for row in executions),
        "executions": len(executions),
        "successful_executions": sum(row["status"] == "ok" for row in executions),
        "family_conditions": len(metrics),
        "independent_genomes": len(datasets),
        "heldout_used": False,
        "warning": "development_ablation;oracle_error_rate_is_not_blind;diagnostic_spread_not_sampling_CI",
    }
    (outdir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["complete"]:
        raise RuntimeError("Calibration matrix incomplete; preserve failed executions")


def run(config_path: Path, datasets: list[Path], outdir: Path) -> None:
    config_path = config_path.resolve()
    datasets = [path.resolve() for path in datasets]
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    outdir.mkdir(parents=True)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir / "source_snapshot"
    provenance = source_manifest(root, snapshot)
    script_target = snapshot / Path(__file__).relative_to(root)
    script_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(__file__, script_target)
    config_target = outdir / "frozen_config.yaml"
    shutil.copyfile(config_path, config_target)
    provenance.update(
        {
            "script_sha256": digest_file(script_target),
            "config_sha256": digest_file(config_target),
            "datasets": {
                str(path): digest_file(path / "generation_receipt.json") for path in datasets
            },
            "scope": "development_known_catalogue_public_quantify_depth_calibration_ablation",
            "resource_note": "single_execution_per_method_condition_not_publication_timing",
        }
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    command = [
        sys.executable,
        str(script_target),
        "--worker",
        "--config",
        str(config_target),
        "--datasets",
        *map(str, datasets),
        "--outdir",
        str(outdir),
    ]
    measured = run_process(
        command,
        outdir / "stdout.log",
        outdir / "stderr.log",
        14400,
        {**os.environ, "PYTHONPATH": str(snapshot)},
        snapshot,
    )
    (outdir / "execution.json").write_text(
        json.dumps({"command": command, **measured}, indent=2) + "\n"
    )
    if measured["exit_code"] != 0 or measured["timed_out"]:
        raise RuntimeError("Quantify calibration failed; inspect preserved logs")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--datasets", required=True, nargs="+", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    (worker if args.worker else run)(args.config, args.datasets, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
