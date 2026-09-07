"""Validate and compactly archive an exact-output quantify calibration replay."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
from pathlib import Path
from typing import Any

from benchmarks.challenge.schema import digest_file, write_table


METHODS = (
    "baseline_total_bases",
    "oracle_error_survival",
    "empirical_controls",
    "empirical_controls_plus_oracle_error",
)
RUN_FILES = (
    "environment.json",
    "execution.json",
    "executions.tsv",
    "summary.tsv",
    "validation.json",
    "frozen_config.yaml",
    "stdout.log",
    "stderr.log",
)
RESOURCE_SUMMARY_FIELDS = {"median_runtime_seconds", "median_peak_rss_mib"}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _copy(source: Path, destination: Path, relative: Path) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    expected = digest_file(source)
    if destination.stat().st_size != source.stat().st_size or digest_file(destination) != expected:
        raise OSError(f"Replay archive copy differs: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source.resolve()),
        "sha256": expected,
        "bytes": destination.stat().st_size,
    }


def _validate_source_snapshot(run: Path, environment: dict[str, Any]) -> None:
    hashes = environment.get("file_hashes")
    snapshot = Path(str(environment.get("source_snapshot", "")))
    if not isinstance(hashes, dict) or not hashes or not snapshot.is_dir():
        raise ValueError(f"Missing source snapshot provenance: {run}")
    expected_digest = hashlib.sha256(
        json.dumps(hashes, sort_keys=True).encode()
    ).hexdigest()
    if expected_digest != environment.get("source_digest"):
        raise ValueError(f"Source digest differs from file manifest: {run}")
    for relative, expected in hashes.items():
        source = snapshot / relative
        if not source.is_file() or digest_file(source) != expected:
            raise ValueError(f"Source snapshot hash differs: {run}/{relative}")


def _validate_run(run: Path) -> dict[str, Any]:
    run = run.resolve()
    required = [run / name for name in (*RUN_FILES, "metrics.tsv")]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Calibration run is incomplete: {', '.join(missing)}")
    environment = _json(run / "environment.json")
    execution = _json(run / "execution.json")
    validation = _json(run / "validation.json")
    executions = _read(run / "executions.tsv")
    metrics = _read(run / "metrics.tsv")
    summaries = _read(run / "summary.tsv")
    if (
        validation.get("complete") is not True
        or execution.get("exit_code") != 0
        or execution.get("timed_out") is not False
        or int(validation.get("executions", -1)) != len(executions)
        or int(validation.get("successful_executions", -1)) != len(executions)
        or int(validation.get("family_conditions", -1)) != len(metrics)
        or not executions
        or not metrics
        or any(row.get("status") != "ok" for row in executions)
        or {row.get("method") for row in executions} != set(METHODS)
    ):
        raise ValueError(f"Calibration run validation is incomplete: {run}")
    execution_keys = {
        (row["seed"], row["condition_id"], row["method"]) for row in executions
    }
    if len(execution_keys) != len(executions):
        raise ValueError(f"Duplicate calibration executions: {run}")
    metric_keys = {
        (row["seed"], row["condition_id"], row["method"], row["family_id"])
        for row in metrics
    }
    if len(metric_keys) != len(metrics):
        raise ValueError(f"Duplicate calibration metric rows: {run}")
    if min(
        float(execution["runtime_seconds"]),
        float(execution["peak_rss_mib"]),
        *(float(row["runtime_seconds"]) for row in executions),
        *(float(row["peak_rss_mib"]) for row in executions),
    ) <= 0:
        raise ValueError(f"Calibration resource measurements must be positive: {run}")
    _validate_source_snapshot(run, environment)
    return {
        "path": run,
        "environment": environment,
        "execution": execution,
        "validation": validation,
        "executions": executions,
        "metrics": metrics,
        "summaries": summaries,
        "execution_keys": execution_keys,
        "metric_keys": metric_keys,
    }


def compare_runs(baseline: Path, replay: Path) -> dict[str, Any]:
    """Require exact scientific parity and summarize descriptive resources."""
    old = _validate_run(baseline)
    new = _validate_run(replay)
    if old["validation"] != new["validation"]:
        raise ValueError("Replay validation receipt differs from baseline")
    if old["execution_keys"] != new["execution_keys"]:
        raise ValueError("Replay execution keys differ from baseline")
    if old["metric_keys"] != new["metric_keys"]:
        raise ValueError("Replay metric keys differ from baseline")
    if digest_file(old["path"] / "frozen_config.yaml") != digest_file(
        new["path"] / "frozen_config.yaml"
    ):
        raise ValueError("Replay frozen configuration differs from baseline")
    baseline_metrics_hash = digest_file(old["path"] / "metrics.tsv")
    replay_metrics_hash = digest_file(new["path"] / "metrics.tsv")
    if baseline_metrics_hash != replay_metrics_hash:
        raise ValueError("Replay scientific metrics are not byte-identical")

    old_summaries = {
        (row["method"], row["coverage"], row["error_model"]): row
        for row in old["summaries"]
    }
    new_summaries = {
        (row["method"], row["coverage"], row["error_model"]): row
        for row in new["summaries"]
    }
    if old_summaries.keys() != new_summaries.keys():
        raise ValueError("Replay summary strata differ from baseline")
    for key in old_summaries:
        old_row = old_summaries[key]
        new_row = new_summaries[key]
        if old_row.keys() != new_row.keys() or any(
            old_row[field] != new_row[field]
            for field in old_row
            if field not in RESOURCE_SUMMARY_FIELDS
        ):
            raise ValueError(f"Replay non-resource summary differs: {key}")

    hash_rows: list[dict[str, object]] = [
        {
            "product": "metrics.tsv",
            "baseline_sha256": baseline_metrics_hash,
            "replay_sha256": replay_metrics_hash,
            "byte_identical": True,
        }
    ]
    seeds = sorted({row["seed"] for row in old["executions"]}, key=int)
    for seed in seeds:
        relative = Path("controls") / f"s{seed}" / "single_copy_kmers.tsv"
        old_path = old["path"] / relative
        new_path = new["path"] / relative
        if not old_path.is_file() or not new_path.is_file():
            raise ValueError(f"Missing control panel: seed {seed}")
        old_hash, new_hash = digest_file(old_path), digest_file(new_path)
        if old_hash != new_hash:
            raise ValueError(f"Replay control panel differs: seed {seed}")
        hash_rows.append(
            {
                "product": relative.as_posix(),
                "baseline_sha256": old_hash,
                "replay_sha256": new_hash,
                "byte_identical": True,
            }
        )
    for seed, condition, method in sorted(old["execution_keys"]):
        relative = (
            Path("runs")
            / f"s{seed}"
            / condition
            / method
            / "output"
            / "copy_number.tsv"
        )
        old_path = old["path"] / relative
        new_path = new["path"] / relative
        if not old_path.is_file() or not new_path.is_file():
            raise ValueError(f"Missing copy-number product: {relative}")
        old_hash, new_hash = digest_file(old_path), digest_file(new_path)
        if old_hash != new_hash:
            raise ValueError(f"Replay copy-number product differs: {relative}")
        hash_rows.append(
            {
                "product": relative.as_posix(),
                "baseline_sha256": old_hash,
                "replay_sha256": new_hash,
                "byte_identical": True,
            }
        )

    old_executions = {
        (row["seed"], row["condition_id"], row["method"]): row
        for row in old["executions"]
    }
    new_executions = {
        (row["seed"], row["condition_id"], row["method"]): row
        for row in new["executions"]
    }
    resource_rows: list[dict[str, object]] = []
    for method in METHODS:
        keys = sorted(key for key in old_executions if key[2] == method)
        old_runtime = [float(old_executions[key]["runtime_seconds"]) for key in keys]
        new_runtime = [float(new_executions[key]["runtime_seconds"]) for key in keys]
        old_rss = [float(old_executions[key]["peak_rss_mib"]) for key in keys]
        new_rss = [float(new_executions[key]["peak_rss_mib"]) for key in keys]
        old_runtime_median = statistics.median(old_runtime)
        new_runtime_median = statistics.median(new_runtime)
        old_rss_median = statistics.median(old_rss)
        new_rss_median = statistics.median(new_rss)
        resource_rows.append(
            {
                "method": method,
                "paired_executions": len(keys),
                "baseline_median_runtime_seconds": old_runtime_median,
                "replay_median_runtime_seconds": new_runtime_median,
                "median_runtime_change_percent": 100
                * (new_runtime_median / old_runtime_median - 1),
                "median_runtime_speedup": old_runtime_median / new_runtime_median,
                "paired_runtime_geometric_mean_ratio": math.exp(
                    statistics.mean(
                        math.log(new_value / old_value)
                        for old_value, new_value in zip(old_runtime, new_runtime)
                    )
                ),
                "baseline_median_peak_rss_mib": old_rss_median,
                "replay_median_peak_rss_mib": new_rss_median,
                "median_peak_rss_change_percent": 100
                * (new_rss_median / old_rss_median - 1),
                "warning": "same_machine_single_replay;descriptive_engineering_resources_not_publication_timing",
            }
        )
    baseline_driver_runtime = float(old["execution"]["runtime_seconds"])
    replay_driver_runtime = float(new["execution"]["runtime_seconds"])
    comparison = {
        "status": "exact_scientific_output_parity",
        "baseline": str(old["path"]),
        "replay": str(new["path"]),
        "baseline_git_head": old["environment"].get("git_head"),
        "replay_git_head": new["environment"].get("git_head"),
        "validation_identical": True,
        "frozen_config_identical": True,
        "metrics_byte_identical": True,
        "metric_rows_compared": len(old["metrics"]),
        "nonresource_summary_cell_differences": 0,
        "copy_number_files_compared": len(old["execution_keys"]),
        "copy_number_file_differences": 0,
        "control_panels_compared": len(seeds),
        "control_panel_differences": 0,
        "exact_hash_rows": len(hash_rows),
        "baseline_metrics_sha256": baseline_metrics_hash,
        "replay_metrics_sha256": replay_metrics_hash,
        "baseline_driver_runtime_seconds": baseline_driver_runtime,
        "replay_driver_runtime_seconds": replay_driver_runtime,
        "driver_runtime_change_percent": 100
        * (replay_driver_runtime / baseline_driver_runtime - 1),
        "driver_speedup": baseline_driver_runtime / replay_driver_runtime,
        "baseline_driver_peak_rss_mib": float(old["execution"]["peak_rss_mib"]),
        "replay_driver_peak_rss_mib": float(new["execution"]["peak_rss_mib"]),
        "warning": "development_matrix_replay;same_machine_single_replay;engineering_not_publication_timing",
    }
    return {
        "comparison": comparison,
        "resources": resource_rows,
        "hashes": hash_rows,
        "baseline": old,
        "replay": new,
    }


def archive(baseline: Path, replay: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    checked = compare_runs(baseline, replay)
    outdir.mkdir(parents=True)
    (outdir / "comparison.json").write_text(
        json.dumps(checked["comparison"], indent=2) + "\n"
    )
    write_table(
        outdir / "resource_comparison.tsv",
        checked["resources"],
        list(checked["resources"][0]),
    )
    write_table(
        outdir / "exact_output_hashes.tsv",
        checked["hashes"],
        list(checked["hashes"][0]),
    )
    (outdir / "README.md").write_text(
        "# Quantify FASTA-survival replay\n\n"
        "This compact archive verifies that the committed all-ACGT FASTA shortcut "
        "preserved the frozen development matrix exactly: all scientific metric "
        "rows, per-execution copy-number products, control panels and non-resource "
        "summary fields agree. Resource values are one historical run and one "
        "same-machine replay, so they support an engineering regression check "
        "rather than a publication timing distribution.\n"
    )
    manifest: list[dict[str, object]] = []
    replay_root = checked["replay"]["path"]
    for name in RUN_FILES:
        relative = Path("replay") / name
        manifest.append(_copy(replay_root / name, outdir / relative, relative))
    for seed in sorted(
        {row["seed"] for row in checked["replay"]["executions"]}, key=int
    ):
        source = replay_root / "controls" / f"s{seed}" / "receipt.json"
        relative = Path("replay") / "controls" / f"s{seed}" / "receipt.json"
        manifest.append(_copy(source, outdir / relative, relative))
    for name in (
        "comparison.json",
        "resource_comparison.tsv",
        "exact_output_hashes.tsv",
        "README.md",
    ):
        path = outdir / name
        manifest.append(
            {
                "file": name,
                "source": "derived_after_exact_output_validation",
                "sha256": digest_file(path),
                "bytes": path.stat().st_size,
            }
        )
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return checked["comparison"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        archive(args.baseline, args.replay, args.outdir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
