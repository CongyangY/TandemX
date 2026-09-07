"""Archive and independently summarize a complete quantify calibration run."""
from __future__ import annotations

import argparse
import csv
import json
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
    "metrics.tsv",
    "summary.tsv",
    "validation.json",
    "frozen_config.yaml",
    "stdout.log",
    "stderr.log",
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _copy(source: Path, destination: Path, relative: Path) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    source_digest = digest_file(source)
    if destination.stat().st_size != source.stat().st_size:
        raise OSError(f"Archive copy size differs: {source}")
    if digest_file(destination) != source_digest:
        raise OSError(f"Archive copy hash differs: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source),
        "sha256": source_digest,
        "bytes": destination.stat().st_size,
    }


def build_decision(
    metrics: list[dict[str, Any]],
    executions: list[dict[str, Any]],
    *,
    control_depth_threshold: float = 2.0,
) -> dict[str, object]:
    """Summarize development effects and a post-hoc candidate without promotion."""
    by_method: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in metrics:
        by_method[str(row["method"])].append(row)
    required = {
        "baseline_total_bases",
        "oracle_error_survival",
        "empirical_controls",
        "empirical_controls_plus_oracle_error",
    }
    if set(by_method) != required:
        raise ValueError("Calibration methods differ from the frozen four-method matrix")
    key_fields = ("seed", "condition_id", "family_id")
    keyed = {
        (str(row["seed"]), str(row["condition_id"]), str(row["family_id"]), str(row["method"])): row
        for row in metrics
    }
    if len(keyed) != len(metrics):
        raise ValueError("Duplicate family-method calibration keys")
    base_keys = {
        key[:3] for key in keyed if key[3] == "baseline_total_bases"
    }
    if any({(*key, method) for key in base_keys} != {
        candidate for candidate in keyed if candidate[3] == method
    } for method in required):
        raise ValueError("Methods do not contain one-to-one family-condition keys")

    aggregates: list[dict[str, object]] = []
    per_seed: list[dict[str, object]] = []
    for method in sorted(required):
        rows = by_method[method]
        run_rows = [row for row in executions if str(row["method"]) == method]
        aggregates.append(
            {
                "method": method,
                "family_conditions": len(rows),
                "mean_absolute_relative_error": statistics.mean(
                    float(row["absolute_relative_error"]) for row in rows
                ),
                "mean_signed_relative_error": statistics.mean(
                    float(row["signed_relative_error"]) for row in rows
                ),
                "mean_estimator_minus_sampling_oracle": statistics.mean(
                    float(row["estimator_minus_sampling_oracle"]) for row in rows
                ),
                "median_runtime_seconds": statistics.median(
                    float(row["runtime_seconds"]) for row in run_rows
                ),
                "median_peak_rss_mib": statistics.median(
                    float(row["peak_rss_mib"]) for row in run_rows
                ),
            }
        )
        for seed in sorted({str(row["seed"]) for row in rows}):
            seed_rows = [row for row in rows if str(row["seed"]) == seed]
            per_seed.append(
                {
                    "method": method,
                    "seed": seed,
                    "family_conditions": len(seed_rows),
                    "mean_absolute_relative_error": statistics.mean(
                        float(row["absolute_relative_error"]) for row in seed_rows
                    ),
                    "mean_signed_relative_error": statistics.mean(
                        float(row["signed_relative_error"]) for row in seed_rows
                    ),
                }
            )

    paired: list[dict[str, object]] = []
    for method in sorted(required - {"baseline_total_bases"}):
        deltas = [
            float(keyed[(*key, "baseline_total_bases")]["absolute_relative_error"])
            - float(keyed[(*key, method)]["absolute_relative_error"])
            for key in base_keys
        ]
        paired.append(
            {
                "method": method,
                "improved": sum(value > 1e-12 for value in deltas),
                "equal": sum(abs(value) <= 1e-12 for value in deltas),
                "worse": sum(value < -1e-12 for value in deltas),
                "mean_absolute_error_reduction": statistics.mean(deltas),
            }
        )

    control_oracle_identical = all(
        str(keyed[(*key, "empirical_controls")]["estimated_copy_number"])
        == str(keyed[(*key, "empirical_controls_plus_oracle_error")]["estimated_copy_number"])
        for key in base_keys
    )
    candidate_errors: list[float] = []
    candidate_by_seed: dict[str, list[float]] = defaultdict(list)
    control_rows_used = 0
    for key in base_keys:
        baseline = keyed[(*key, "baseline_total_bases")]
        control = keyed[(*key, "empirical_controls")]
        use_control = float(control["control_mean_depth"]) >= control_depth_threshold
        selected = control if use_control else baseline
        error = float(selected["absolute_relative_error"])
        candidate_errors.append(error)
        candidate_by_seed[key[0]].append(error)
        control_rows_used += int(use_control)
    baseline_seed = {
        str(row["seed"]): float(row["mean_absolute_relative_error"])
        for row in per_seed
        if row["method"] == "baseline_total_bases"
    }
    candidate_seed = {
        seed: statistics.mean(values) for seed, values in sorted(candidate_by_seed.items())
    }
    seed_improvements = {
        seed: baseline_seed[seed] - candidate_seed[seed] for seed in candidate_seed
    }
    aggregate_by_method = {str(row["method"]): row for row in aggregates}
    return {
        "status": "development_candidate_selected_not_promoted",
        "aggregates": aggregates,
        "per_seed": per_seed,
        "paired_against_baseline": paired,
        "control_plus_oracle_estimates_identical": control_oracle_identical,
        "posthoc_candidate": {
            "rule": "empirical_controls_when_control_mean_depth_at_least_threshold_else_total_bases",
            "control_mean_depth_threshold": control_depth_threshold,
            "mean_absolute_relative_error": statistics.mean(candidate_errors),
            "baseline_mean_absolute_relative_error": aggregate_by_method[
                "baseline_total_bases"
            ]["mean_absolute_relative_error"],
            "control_family_rows_used": control_rows_used,
            "total_family_rows": len(base_keys),
            "per_seed_mean_absolute_relative_error": candidate_seed,
            "per_seed_improvement_over_baseline": seed_improvements,
            "all_seed_improvements_positive": all(value > 0 for value in seed_improvements.values()),
        },
        "promotion": {
            "passed": False,
            "reason": "development_only;posthoc_rule_requires_frozen_untouched_genomes;simulation_truth_assisted_controls_not_blind_real_controls",
        },
    }


def validate_run(run: Path) -> tuple[dict[str, object], list[dict[str, str]], list[dict[str, str]]]:
    missing = [name for name in ROOT_FILES if not (run / name).is_file()]
    if missing:
        raise ValueError(f"Calibration run lacks required files: {', '.join(missing)}")
    validation = json.loads((run / "validation.json").read_text())
    executions = _read(run / "executions.tsv")
    metrics = _read(run / "metrics.tsv")
    summary = _read(run / "summary.tsv")
    if validation != {
        "complete": True,
        "executions": 108,
        "successful_executions": 108,
        "family_conditions": 5940,
        "independent_genomes": 3,
        "heldout_used": False,
        "warning": "development_ablation;oracle_error_rate_is_not_blind;diagnostic_spread_not_sampling_CI",
    }:
        raise ValueError("Calibration validation receipt differs from the frozen matrix")
    if len(executions) != 108 or any(row["status"] != "ok" for row in executions):
        raise ValueError("Calibration executions are incomplete or failed")
    if len(metrics) != 5940 or len(summary) != 36:
        raise ValueError("Calibration result row counts differ from the frozen matrix")
    config = yaml.safe_load((run / "frozen_config.yaml").read_text())
    if config.get("benchmark_id") != "quantify_calibration_development_v1":
        raise ValueError("Unexpected calibration benchmark ID")
    environment = json.loads((run / "environment.json").read_text())
    if environment.get("config_sha256") != digest_file(run / "frozen_config.yaml"):
        raise ValueError("Environment/config hash mismatch")
    snapshot = Path(str(environment["source_snapshot"]))
    for relative, expected in environment["file_hashes"].items():
        source = snapshot / relative
        if not source.is_file() or digest_file(source) != expected:
            raise ValueError(f"Source snapshot hash mismatch: {relative}")
    if environment.get("script_sha256") != digest_file(
        snapshot / "benchmarks/scripts/evaluate_quantify_calibration.py"
    ):
        raise ValueError("Runner snapshot hash mismatch")
    build_decision(metrics, executions)
    return validation, executions, metrics


def archive(run: Path, outdir: Path) -> dict[str, object]:
    run = run.resolve()
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    validation, executions, metrics = validate_run(run)
    outdir.mkdir(parents=True)
    manifest: list[dict[str, object]] = []
    for name in ROOT_FILES:
        relative = Path("run") / name
        manifest.append(_copy(run / name, outdir / relative, relative))
    for seed in (6301, 6302, 6303):
        for name in ("single_copy_kmers.tsv", "receipt.json"):
            relative = Path("controls") / f"s{seed}" / name
            manifest.append(_copy(run / relative, outdir / relative, relative))

    command_rows: list[dict[str, object]] = []
    by_execution = {
        (row["seed"], row["condition_id"], row["method"]): row for row in executions
    }
    for key, execution in sorted(by_execution.items()):
        seed, condition, method = key
        source = run / "runs" / f"s{seed}" / condition / method
        command_path = source / "command.json"
        receipt_path = source / "receipt.json"
        output = source / "output"
        required = (
            command_path,
            receipt_path,
            output / "copy_number.tsv",
            output / "run_config.yaml",
            output / "run.log",
            source / "stdout.log",
            source / "stderr.log",
        )
        if any(not path.is_file() for path in required):
            raise ValueError(f"Execution artifacts incomplete: {key}")
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("status") != "ok" or int(receipt.get("exit_code", -1)) != 0:
            raise ValueError(f"Execution receipt disagrees with matrix: {key}")
        command_rows.append(
            {
                "seed": seed,
                "condition_id": condition,
                "method": method,
                "command_json": json.dumps(json.loads(command_path.read_text()), separators=(",", ":")),
                "command_sha256": digest_file(command_path),
                "receipt_sha256": digest_file(receipt_path),
                "copy_number_sha256": digest_file(output / "copy_number.tsv"),
                "run_config_sha256": digest_file(output / "run_config.yaml"),
                "run_log_sha256": digest_file(output / "run.log"),
                "stdout_sha256": digest_file(source / "stdout.log"),
                "stderr_sha256": digest_file(source / "stderr.log"),
            }
        )
    write_table(outdir / "execution_artifacts.tsv", command_rows, list(command_rows[0]))
    decision = build_decision(metrics, executions)
    (outdir / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    archive_summary = {
        "benchmark_id": "quantify_calibration_development_v1",
        "validation": validation,
        "execution_artifact_rows": len(command_rows),
        "decision_status": decision["status"],
        "warning": "development_only;not_heldout;oracle_input_not_blind;diagnostic_spread_not_sampling_CI",
    }
    (outdir / "archive_summary.json").write_text(json.dumps(archive_summary, indent=2) + "\n")
    (outdir / "README.md").write_text(
        "# Quantify calibration development evidence\n\n"
        "This compact archive retains the complete three-genome, 108-execution "
        "public-command matrix, 5,940 family rows, controls, source/config hashes, "
        "per-execution artifact hashes and a recomputed development decision. The "
        "oracle error rate is not blind, control selection uses simulation truth, "
        "and the diagnostic spread is not a sampling confidence interval. The "
        "post-hoc depth-gated control rule is a candidate for a future frozen "
        "held-out test; it is not promoted here.\n"
    )
    for name in ("execution_artifacts.tsv", "decision.json", "archive_summary.json", "README.md"):
        path = outdir / name
        manifest.append(
            {
                "file": name,
                "source": "derived_after_validation",
                "sha256": digest_file(path),
                "bytes": path.stat().st_size,
            }
        )
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return archive_summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        archive(args.run, args.outdir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
