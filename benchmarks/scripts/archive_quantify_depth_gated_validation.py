"""Validate and compact the frozen depth-gated quantify validation evidence."""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml

from benchmarks.challenge.schema import digest_file, write_table


ROOT_FILES = (
    "environment.json",
    "execution.json",
    "executions.tsv",
    "raw_metrics.tsv",
    "metrics.tsv",
    "paired.tsv",
    "aggregate_summary.tsv",
    "seed_summary.tsv",
    "coverage_summary.tsv",
    "gate_results.json",
    "validation.json",
    "frozen_config.yaml",
    "stdout.log",
    "stderr.log",
)
METHODS = ("baseline_total_bases", "empirical_controls", "depth_gated_controls")


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _copy(source: Path, destination: Path, relative: Path) -> dict[str, object]:
    if not source.is_file():
        raise ValueError(f"Required archive source is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    source_hash = digest_file(source)
    if destination.stat().st_size != source.stat().st_size or digest_file(destination) != source_hash:
        raise OSError(f"Archive copy differs from source: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source),
        "sha256": source_hash,
        "bytes": destination.stat().st_size,
    }


def recompute_decision(
    metrics: list[dict[str, str]],
    paired: list[dict[str, str]],
    executions: list[dict[str, str]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Recompute the principal result and descriptive resource summaries."""
    by_method: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        by_method[row["method"]].append(row)
    if set(by_method) != set(METHODS) or any(len(by_method[name]) != 1485 for name in METHODS):
        raise ValueError("Validation metrics do not contain 1,485 rows for each method")

    aggregate = {
        method: statistics.mean(
            float(row["absolute_relative_error"]) for row in by_method[method]
        )
        for method in METHODS
    }
    bias = {
        method: statistics.mean(
            float(row["signed_relative_error"]) for row in by_method[method]
        )
        for method in METHODS
    }
    seed_reductions: dict[str, float] = {}
    for seed in sorted({row["seed"] for row in metrics}, key=int):
        values = {
            method: statistics.mean(
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if row["seed"] == seed
            )
            for method in METHODS
        }
        seed_reductions[seed] = values["baseline_total_bases"] - values["depth_gated_controls"]
    coverage_reductions: dict[str, float] = {}
    for coverage in sorted({row["coverage"] for row in metrics}, key=float):
        values = {
            method: statistics.mean(
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if row["coverage"] == coverage
            )
            for method in METHODS
        }
        coverage_reductions[coverage] = (
            values["baseline_total_bases"] - values["depth_gated_controls"]
        )
    outcomes = {
        name: sum(row["outcome"] == name for row in paired)
        for name in ("improved", "equal", "worse")
    }
    selections = {
        (row["seed"], row["condition_id"]): row["candidate_selected_method"]
        for row in paired
    }
    if len(selections) != 27:
        raise ValueError("Depth-gated selection does not contain 27 unique conditions")
    selection_counts = {
        method: sum(value == method for value in selections.values())
        for method in ("baseline_total_bases", "empirical_controls")
    }

    resource_rows: list[dict[str, object]] = []
    for method in ("baseline_total_bases", "empirical_controls"):
        method_rows = [row for row in executions if row["method"] == method]
        for coverage in ("1", "5", "20"):
            rows = [row for row in method_rows if row["coverage"] == coverage]
            if len(rows) != 9:
                raise ValueError(f"Expected nine resource rows for {method} at {coverage}x")
            runtimes = [float(row["runtime_seconds"]) for row in rows]
            memories = [float(row["peak_rss_mib"]) for row in rows]
            resource_rows.append(
                {
                    "method": method,
                    "coverage": coverage,
                    "executions": len(rows),
                    "median_runtime_seconds": statistics.median(runtimes),
                    "minimum_runtime_seconds": min(runtimes),
                    "maximum_runtime_seconds": max(runtimes),
                    "median_peak_rss_mib": statistics.median(memories),
                    "minimum_peak_rss_mib": min(memories),
                    "maximum_peak_rss_mib": max(memories),
                }
            )
    decision = {
        "status": "passed_frozen_simulation_gates",
        "aggregate_mean_absolute_relative_error": aggregate,
        "aggregate_mean_signed_relative_error": bias,
        "aggregate_mare_reduction": aggregate["baseline_total_bases"]
        - aggregate["depth_gated_controls"],
        "candidate_minus_empirical_controls_mare": aggregate["depth_gated_controls"]
        - aggregate["empirical_controls"],
        "per_seed_mare_reduction": seed_reductions,
        "per_coverage_mare_reduction": coverage_reductions,
        "paired_outcomes": outcomes,
        "selected_condition_counts": selection_counts,
        "interpretation": (
            "the frozen gate improves total-bases normalization on every validation genome; "
            "ungated empirical controls retain slightly lower aggregate MARE in this split"
        ),
        "warning": (
            "simulation_truth_assisted_controls;known_catalogue;conditions_within_genome_dependent;"
            "single_execution_resources_not_publication_timing;not_biological_truth"
        ),
    }
    return decision, resource_rows


def _assert_close(observed: float, expected: float, label: str) -> None:
    if not math.isclose(observed, expected, rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"Recomputed {label} differs from gate receipt")


def validate_run(
    run: Path, input_audit: Path
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    dict[str, object],
    list[dict[str, object]],
]:
    missing = [name for name in ROOT_FILES if not (run / name).is_file()]
    if missing:
        raise ValueError(f"Validation run lacks required files: {', '.join(missing)}")
    validation = json.loads((run / "validation.json").read_text())
    expected_validation = {
        "complete": True,
        "executions": 54,
        "successful_executions": 54,
        "raw_family_conditions": 2970,
        "candidate_family_conditions": 1485,
        "independent_genomes": 3,
        "validation_used": True,
        "heldout_used": False,
        "promotion_status": "passed",
        "warning": "frozen_validation_same_factorial_process;conditions_within_genome_are_dependent",
    }
    if validation != expected_validation:
        raise ValueError("Validation receipt differs from the completed frozen matrix")
    executions = _read(run / "executions.tsv")
    raw_metrics = _read(run / "raw_metrics.tsv")
    metrics = _read(run / "metrics.tsv")
    paired = _read(run / "paired.tsv")
    if len(executions) != 54 or any(row["status"] != "ok" for row in executions):
        raise ValueError("Execution matrix is incomplete or contains failures")
    if len({(row["seed"], row["condition_id"], row["method"]) for row in executions}) != 54:
        raise ValueError("Execution matrix contains duplicate keys")
    if len(raw_metrics) != 2970 or len(metrics) != 4455 or len(paired) != 1485:
        raise ValueError("Scientific table row counts differ from the frozen matrix")

    config = yaml.safe_load((run / "frozen_config.yaml").read_text())
    if config.get("benchmark_id") != "quantify_depth_gated_validation_v1":
        raise ValueError("Unexpected validation benchmark ID")
    environment = json.loads((run / "environment.json").read_text())
    if environment.get("config_sha256") != digest_file(run / "frozen_config.yaml"):
        raise ValueError("Environment/config hash mismatch")
    snapshot = Path(environment["source_snapshot"])
    for relative, expected in environment["file_hashes"].items():
        source = snapshot / relative
        if not source.is_file() or digest_file(source) != expected:
            raise ValueError(f"Source snapshot hash mismatch: {relative}")
    for relative, expected in environment["helper_hashes"].items():
        source = snapshot / relative
        if not source.is_file() or digest_file(source) != expected:
            raise ValueError(f"Validation helper hash mismatch: {relative}")

    audit = json.loads(input_audit.read_text())
    if not audit.get("complete") or audit.get("status") != "passed" or audit.get("dataset_count") != 3:
        raise ValueError("Input dataset audit did not pass")
    audit_receipts = {
        row["dataset"]: row["generation_receipt_sha256"] for row in audit["datasets"]
    }
    if audit_receipts != environment["datasets"]:
        raise ValueError("Input audit dataset receipts differ from the validation environment")

    decision, resource_rows = recompute_decision(metrics, paired, executions)
    gates = json.loads((run / "gate_results.json").read_text())
    if gates.get("status") != "passed" or gates.get("failed_gate_names"):
        raise ValueError("Frozen scientific gates did not pass")
    _assert_close(
        float(decision["aggregate_mare_reduction"]),
        float(gates["aggregate_mare_reduction"]),
        "aggregate MARE reduction",
    )
    _assert_close(
        min(float(value) for value in decision["per_seed_mare_reduction"].values()),
        float(gates["minimum_seed_mare_reduction"]),
        "minimum seed MARE reduction",
    )
    _assert_close(
        min(float(value) for value in decision["per_coverage_mare_reduction"].values()),
        float(gates["minimum_coverage_mare_reduction"]),
        "minimum coverage MARE reduction",
    )
    if decision["paired_outcomes"] != {
        "improved": gates["paired_improved"],
        "equal": gates["paired_equal"],
        "worse": gates["paired_worse"],
    } or decision["selected_condition_counts"] != gates["selected_condition_counts"]:
        raise ValueError("Recomputed paired outcomes or branch selections differ")
    return validation, environment, executions, metrics, paired, decision, resource_rows


def archive(run: Path, input_audit: Path, outdir: Path) -> dict[str, object]:
    run = run.resolve()
    input_audit = input_audit.resolve()
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    (
        validation,
        environment,
        executions,
        metrics,
        paired,
        decision,
        resource_rows,
    ) = validate_run(run, input_audit)
    outdir.mkdir(parents=True)
    manifest: list[dict[str, object]] = []
    for name in ROOT_FILES:
        relative = Path("run") / name
        manifest.append(_copy(run / name, outdir / relative, relative))
    manifest.append(_copy(input_audit, outdir / "input_audit.json", Path("input_audit.json")))

    snapshot = Path(environment["source_snapshot"])
    for relative_text in environment["helper_hashes"]:
        relative = Path("source_snapshot") / relative_text
        manifest.append(_copy(snapshot / relative_text, outdir / relative, relative))

    config = yaml.safe_load((run / "frozen_config.yaml").read_text())
    dataset_by_seed = {}
    for dataset_text, expected_hash in environment["datasets"].items():
        dataset = Path(dataset_text)
        receipt = json.loads((dataset / "generation_receipt.json").read_text())
        if digest_file(dataset / "generation_receipt.json") != expected_hash:
            raise ValueError(f"Generation receipt changed: {dataset}")
        dataset_by_seed[int(receipt["seed"])] = (dataset, receipt)
    if list(sorted(dataset_by_seed)) != [int(seed) for seed in config["seeds"]]:
        raise ValueError("Validation datasets differ from frozen seeds")
    for seed, (dataset, receipt) in sorted(dataset_by_seed.items()):
        sources = [
            Path("generation_receipt.json"),
            Path("run_config.json"),
            Path("length_histogram.tsv"),
            Path("genome/manifest.json"),
            *[
                Path("reads") / row["condition_id"] / "manifest.json"
                for row in receipt["conditions_completed"]
            ],
        ]
        for source_relative in sources:
            relative = Path("datasets") / f"s{seed}" / source_relative
            manifest.append(_copy(dataset / source_relative, outdir / relative, relative))

    for seed in config["seeds"]:
        for name in ("single_copy_kmers.tsv", "receipt.json"):
            source_relative = Path("controls") / f"s{seed}" / name
            manifest.append(_copy(run / source_relative, outdir / source_relative, source_relative))

    artifact_rows: list[dict[str, object]] = []
    for execution in sorted(
        executions, key=lambda row: (int(row["seed"]), row["condition_id"], row["method"])
    ):
        relative_run = (
            Path("runs")
            / f"s{execution['seed']}"
            / execution["condition_id"]
            / execution["method"]
        )
        source = run / relative_run
        paths = {
            "command": source / "command.json",
            "receipt": source / "receipt.json",
            "copy_number": source / "output/copy_number.tsv",
            "run_config": source / "output/run_config.yaml",
            "run_log": source / "output/run.log",
            "stdout": source / "stdout.log",
            "stderr": source / "stderr.log",
        }
        if any(not path.is_file() for path in paths.values()):
            raise ValueError(f"Execution artifacts are incomplete: {relative_run}")
        receipt = json.loads(paths["receipt"].read_text())
        if receipt.get("status") != "ok" or int(receipt.get("exit_code", -1)) != 0:
            raise ValueError(f"Execution receipt is not successful: {relative_run}")
        artifact_rows.append(
            {
                "seed": execution["seed"],
                "condition_id": execution["condition_id"],
                "method": execution["method"],
                "command_json": json.dumps(
                    json.loads(paths["command"].read_text()), separators=(",", ":")
                ),
                **{f"{name}_sha256": digest_file(path) for name, path in paths.items()},
            }
        )
    write_table(outdir / "execution_artifacts.tsv", artifact_rows, list(artifact_rows[0]))
    write_table(outdir / "resource_summary.tsv", resource_rows, list(resource_rows[0]))
    (outdir / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    archive_summary = {
        "benchmark_id": config["benchmark_id"],
        "validation": validation,
        "execution_artifact_rows": len(artifact_rows),
        "metric_rows": len(metrics),
        "paired_rows": len(paired),
        "dataset_receipts": len(dataset_by_seed),
        "decision_status": decision["status"],
        "manifest_entries_excluding_self": 0,
        "warning": decision["warning"],
    }
    (outdir / "README.md").write_text(
        "# Frozen depth-gated quantify validation evidence\n\n"
        "This compact archive retains the complete 54-command, three-genome "
        "validation matrix, 4,455 metric rows, 1,485 paired family outcomes, "
        "all compact dataset/config/manifests, the independent 93-payload input "
        "audit, controls, source hashes and per-execution artifact hashes. The "
        "predeclared rule passed all nine simulation gates and improved total-bases "
        "normalization on every validation genome. Ungated empirical controls had "
        "slightly lower aggregate MARE in this split. Control selection used "
        "simulation truth, catalogues were supplied, conditions within a genome "
        "are dependent, and one execution per method-condition is not a final "
        "timing distribution or biological validation. Seeds 6401--6403 are "
        "consumed and must not be used for tuning.\n"
    )
    for name in (
        "execution_artifacts.tsv",
        "resource_summary.tsv",
        "decision.json",
        "README.md",
    ):
        path = outdir / name
        manifest.append(
            {
                "file": name,
                "source": "derived_after_validation",
                "sha256": digest_file(path),
                "bytes": path.stat().st_size,
            }
        )
    archive_summary["manifest_entries_excluding_self"] = len(manifest) + 1
    (outdir / "archive_summary.json").write_text(json.dumps(archive_summary, indent=2) + "\n")
    manifest.append(
        {
            "file": "archive_summary.json",
            "source": "derived_after_validation",
            "sha256": digest_file(outdir / "archive_summary.json"),
            "bytes": (outdir / "archive_summary.json").stat().st_size,
        }
    )
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return archive_summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--input-audit", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        archive(args.run, args.input_audit, args.outdir)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
