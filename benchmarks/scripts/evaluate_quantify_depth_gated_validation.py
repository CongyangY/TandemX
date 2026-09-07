"""Run the frozen untouched-genome validation of depth-gated quantification."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from benchmarks.scripts.evaluate_quantify_calibration import select_single_copy_controls


METHODS = ("baseline_total_bases", "empirical_controls")
CANDIDATE = "depth_gated_controls"


def _gate(name: str, observed: float, operator: str, threshold: float) -> dict[str, object]:
    comparisons = {
        "<=": observed <= threshold,
        ">=": observed >= threshold,
        ">": observed > threshold,
    }
    if operator not in comparisons:
        raise ValueError(f"Unsupported gate operator: {operator}")
    return {
        "name": name,
        "observed": observed,
        "operator": operator,
        "threshold": threshold,
        "passed": comparisons[operator],
    }


def build_candidate(
    metrics: list[dict[str, Any]], threshold: float
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Apply the frozen condition-level depth gate without fitting validation data."""
    if threshold < 0:
        raise ValueError("Control-depth threshold must be nonnegative")
    keys = {
        (
            str(row["seed"]),
            str(row["condition_id"]),
            str(row["family_id"]),
            str(row["method"]),
        ): row
        for row in metrics
    }
    if len(keys) != len(metrics):
        raise ValueError("Duplicate validation family-method keys")
    base_keys = {key[:3] for key in keys if key[3] == METHODS[0]}
    expected = {(*key, method) for key in base_keys for method in METHODS}
    if not base_keys or set(keys) != expected:
        raise ValueError("Validation methods do not have one-to-one family keys")

    candidate_rows: list[dict[str, Any]] = []
    paired_rows: list[dict[str, Any]] = []
    selection_by_condition: dict[tuple[str, str], str] = {}
    for key in sorted(base_keys):
        baseline = keys[(*key, METHODS[0])]
        controls = keys[(*key, METHODS[1])]
        control_depth = float(controls["control_mean_depth"])
        selected_method = METHODS[1] if control_depth >= threshold else METHODS[0]
        condition_key = key[:2]
        previous = selection_by_condition.setdefault(condition_key, selected_method)
        if previous != selected_method:
            raise ValueError("Candidate selection differs among families in one condition")
        selected = dict(controls if selected_method == METHODS[1] else baseline)
        selected["method"] = CANDIDATE
        selected["candidate_selected_method"] = selected_method
        selected["candidate_control_mean_depth"] = control_depth
        candidate_rows.append(selected)
        baseline_error = float(baseline["absolute_relative_error"])
        candidate_error = float(selected["absolute_relative_error"])
        delta = baseline_error - candidate_error
        paired_rows.append(
            {
                "seed": key[0],
                "condition_id": key[1],
                "coverage": baseline["coverage"],
                "error_model": baseline["error_model"],
                "family_id": key[2],
                "baseline_absolute_relative_error": baseline_error,
                "candidate_absolute_relative_error": candidate_error,
                "baseline_minus_candidate_absolute_error": delta,
                "outcome": (
                    "improved"
                    if delta > 1e-12
                    else "worse"
                    if delta < -1e-12
                    else "equal"
                ),
                "candidate_selected_method": selected_method,
                "control_mean_depth": control_depth,
            }
        )
    return candidate_rows, paired_rows


def summarize(
    metrics: list[dict[str, Any]], group_field: str | None
) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in metrics:
        stratum = "all" if group_field is None else str(row[group_field])
        groups[(str(row["method"]), stratum)].append(row)
    output: list[dict[str, object]] = []
    for (method, stratum), rows in sorted(groups.items()):
        output.append(
            {
                "method": method,
                "stratum": stratum,
                "family_conditions": len(rows),
                "mean_absolute_relative_error": statistics.mean(
                    float(row["absolute_relative_error"]) for row in rows
                ),
                "mean_signed_relative_error": statistics.mean(
                    float(row["signed_relative_error"]) for row in rows
                ),
            }
        )
    return output


def evaluate_gates(
    metrics: list[dict[str, Any]],
    paired: list[dict[str, Any]],
    executions: list[dict[str, Any]],
    limits: dict[str, Any],
) -> dict[str, Any]:
    aggregate = {(row["method"], row["stratum"]): row for row in summarize(metrics, None)}
    by_seed = {(row["method"], row["stratum"]): row for row in summarize(metrics, "seed")}
    by_coverage = {
        (row["method"], row["stratum"]): row for row in summarize(metrics, "coverage")
    }
    baseline_mare = float(aggregate[(METHODS[0], "all")]["mean_absolute_relative_error"])
    candidate_mare = float(aggregate[(CANDIDATE, "all")]["mean_absolute_relative_error"])
    seed_reductions = [
        float(by_seed[(METHODS[0], seed)]["mean_absolute_relative_error"])
        - float(by_seed[(CANDIDATE, seed)]["mean_absolute_relative_error"])
        for seed in sorted({str(row["seed"]) for row in metrics}, key=int)
    ]
    coverage_reductions = [
        float(by_coverage[(METHODS[0], coverage)]["mean_absolute_relative_error"])
        - float(by_coverage[(CANDIDATE, coverage)]["mean_absolute_relative_error"])
        for coverage in sorted(
            {str(row["coverage"]) for row in metrics}, key=float
        )
    ]
    improved = sum(row["outcome"] == "improved" for row in paired)
    equal = sum(row["outcome"] == "equal" for row in paired)
    selection = {
        (str(row["seed"]), str(row["condition_id"])): str(
            row["candidate_selected_method"]
        )
        for row in paired
    }
    selected_counts = {
        method: sum(value == method for value in selection.values()) for method in METHODS
    }
    failed_executions = sum(row.get("status") != "ok" for row in executions)
    gates = [
        _gate(
            "failed_executions",
            float(failed_executions),
            "<=",
            float(limits["failed_executions_max"]),
        ),
        _gate(
            "aggregate_mare_reduction",
            baseline_mare - candidate_mare,
            ">=",
            float(limits["aggregate_mare_reduction_min"]),
        ),
        _gate(
            "candidate_mare",
            candidate_mare,
            "<=",
            float(limits["candidate_mare_max"]),
        ),
        _gate(
            "minimum_seed_mare_reduction",
            min(seed_reductions),
            ">",
            float(limits["minimum_seed_mare_reduction_strictly_greater_than"]),
        ),
        _gate(
            "minimum_coverage_mare_reduction",
            min(coverage_reductions),
            ">=",
            float(limits["minimum_coverage_mare_reduction_min"]),
        ),
        _gate(
            "paired_improved_fraction",
            improved / len(paired),
            ">=",
            float(limits["paired_improved_fraction_min"]),
        ),
        _gate(
            "paired_nonworse_fraction",
            (improved + equal) / len(paired),
            ">=",
            float(limits["paired_nonworse_fraction_min"]),
        ),
        _gate(
            "baseline_selected_condition_count",
            float(selected_counts[METHODS[0]]),
            ">=",
            float(limits["selected_condition_count_per_branch_min"]),
        ),
        _gate(
            "controls_selected_condition_count",
            float(selected_counts[METHODS[1]]),
            ">=",
            float(limits["selected_condition_count_per_branch_min"]),
        ),
    ]
    return {
        "status": "passed" if all(row["passed"] for row in gates) else "failed",
        "gates": gates,
        "failed_gate_names": [row["name"] for row in gates if not row["passed"]],
        "baseline_mean_absolute_relative_error": baseline_mare,
        "candidate_mean_absolute_relative_error": candidate_mare,
        "aggregate_mare_reduction": baseline_mare - candidate_mare,
        "minimum_seed_mare_reduction": min(seed_reductions),
        "minimum_coverage_mare_reduction": min(coverage_reductions),
        "paired_improved": improved,
        "paired_equal": equal,
        "paired_worse": len(paired) - improved - equal,
        "selected_condition_counts": selected_counts,
        "warning": "frozen_validation_same_factorial_process;simulation_truth_assisted_controls;conditional_known_catalogue_scope",
    }


def worker(config_path: Path, datasets: list[Path], outdir: Path) -> None:
    config = yaml.safe_load(config_path.read_text())
    if (
        config.get("split") != "validation"
        or tuple(config.get("methods", ())) != METHODS
        or config.get("candidate", {}).get("rule")
        != "empirical_controls_when_control_mean_depth_at_least_threshold_else_total_bases"
    ):
        raise ValueError("Require the frozen depth-gated validation configuration")
    seeds = [int(value) for value in config["seeds"]]
    forbidden = {int(value) for value in config["forbidden_development_seeds"]}
    if set(seeds) & forbidden or len(set(seeds)) != len(seeds):
        raise ValueError("Validation seeds overlap development or are duplicated")
    generations = [json.loads((path / "generation_receipt.json").read_text()) for path in datasets]
    if [int(row["seed"]) for row in generations] != seeds:
        raise ValueError("Dataset order/seeds differ from frozen validation config")
    if any(
        not row["complete"]
        or row["split"] != "validation"
        or row.get("validation_used") is not True
        or row.get("heldout_used") is not False
        or row.get("config_sha256")
        != config["dataset_generation"]["config_sha256"]
        or row.get("histogram_sha256")
        != config["dataset_generation"]["length_histogram_sha256"]
        or row.get("source_sha256")
        != config["dataset_generation"]["source_sha256"]
        or len(row["conditions_completed"]) != int(config["condition_count_per_seed"])
        for row in generations
    ):
        raise ValueError("Validation datasets are incomplete or have wrong split semantics")

    executions: list[dict[str, Any]] = []
    raw_metrics: list[dict[str, Any]] = []
    for dataset, generation in zip(datasets, generations):
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
        (control_dir / "receipt.json").write_text(
            json.dumps(control_receipt, indent=2) + "\n"
        )
        truth = {
            row["family_id"]: row
            for row in read_table(genome_dir / "truth_copy_number.tsv")
        }
        if len(truth) != int(generation["family_count"]):
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
                if method == METHODS[1]:
                    command.extend(["--single-copy-kmers", str(controls_path)])
                (run_dir / "command.json").write_text(
                    json.dumps(command, indent=2) + "\n"
                )
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
                    **measured,
                }
                execution["status"] = (
                    "ok"
                    if execution["exit_code"] == 0 and not execution["timed_out"]
                    else "failed"
                )
                executions.append(execution)
                write_table(outdir / "executions.tsv", executions, list(executions[0]))
                (run_dir / "receipt.json").write_text(
                    json.dumps(execution, indent=2) + "\n"
                )
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
                    raw_metrics.append(
                        {
                            "seed": generation["seed"],
                            "condition_id": condition["condition_id"],
                            "coverage": condition["coverage"],
                            "error_model": condition["label"],
                            "method": method,
                            "family_id": family_id,
                            "period": truth_row["period"],
                            "unit_substitution_rate": truth_row[
                                "unit_substitution_rate"
                            ],
                            "truth_copies": copies,
                            "estimated_copy_number": value,
                            "signed_relative_error": (value - copies) / copies,
                            "absolute_relative_error": abs(value - copies) / copies,
                            "normalization_method": estimate["normalization_method"],
                            "control_mean_depth": estimate.get(
                                "single_copy_control_mean_depth", "NA"
                            ),
                            "control_zero_fraction": estimate.get(
                                "single_copy_control_zero_fraction", "NA"
                            ),
                        }
                    )
                write_table(
                    outdir / "raw_metrics.tsv", raw_metrics, list(raw_metrics[0])
                )
    expected_executions = len(datasets) * int(config["condition_count_per_seed"]) * len(METHODS)
    complete = len(executions) == expected_executions and all(
        row["status"] == "ok" for row in executions
    )
    if complete:
        threshold = float(config["candidate"]["control_mean_depth_threshold"])
        candidate, paired = build_candidate(raw_metrics, threshold)
        metrics = [*raw_metrics, *candidate]
        metric_fields = list(
            dict.fromkeys(field for row in metrics for field in row)
        )
        write_table(outdir / "metrics.tsv", metrics, metric_fields)
        write_table(outdir / "paired.tsv", paired, list(paired[0]))
        for name, field in (("aggregate_summary.tsv", None), ("seed_summary.tsv", "seed"), ("coverage_summary.tsv", "coverage")):
            rows = summarize(metrics, field)
            write_table(outdir / name, rows, list(rows[0]))
        gate_result = evaluate_gates(
            metrics, paired, executions, config["acceptance_gates"]
        )
    else:
        gate_result = {
            "status": "failed",
            "gates": [],
            "failed_gate_names": ["incomplete_execution_matrix"],
            "warning": "preserved_failed_execution_matrix;scientific_metrics_not_completed",
        }
    gate_result.update(
        {
            "benchmark_id": config["benchmark_id"],
            "config_sha256": digest_file(config_path),
            "candidate_rule": config["candidate"],
            "validation_seeds_consumed": seeds,
        }
    )
    (outdir / "gate_results.json").write_text(
        json.dumps(gate_result, indent=2) + "\n"
    )
    validation = {
        "complete": complete,
        "executions": len(executions),
        "successful_executions": sum(row["status"] == "ok" for row in executions),
        "raw_family_conditions": len(raw_metrics),
        "candidate_family_conditions": len(candidate) if complete else 0,
        "independent_genomes": len(datasets),
        "validation_used": True,
        "heldout_used": False,
        "promotion_status": gate_result["status"],
        "warning": "frozen_validation_same_factorial_process;conditions_within_genome_are_dependent",
    }
    (outdir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")


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
    helpers = (
        Path(__file__).resolve(),
        root / "benchmarks" / "scripts" / "evaluate_quantify_calibration.py",
    )
    helper_hashes = {}
    for helper in helpers:
        relative = helper.relative_to(root)
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(helper, target)
        helper_hashes[relative.as_posix()] = digest_file(target)
    config_target = outdir / "frozen_config.yaml"
    shutil.copyfile(config_path, config_target)
    provenance.update(
        {
            "helper_hashes": helper_hashes,
            "config_sha256": digest_file(config_target),
            "datasets": {
                str(path): digest_file(path / "generation_receipt.json")
                for path in datasets
            },
            "scope": "frozen_untouched_genome_depth_gated_quantify_validation",
            "resource_note": "single_execution_per_method_condition_not_publication_timing",
        }
    )
    (outdir / "environment.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    script_target = snapshot / Path(__file__).resolve().relative_to(root)
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
        raise RuntimeError("Depth-gated validation failed; inspect preserved logs")
    if not (outdir / "validation.json").is_file():
        raise RuntimeError("Depth-gated validation did not produce a receipt")


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
