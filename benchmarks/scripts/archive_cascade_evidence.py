"""Archive split cascade benchmarks after strict provenance and parity checks."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import shutil
from typing import Iterable

from benchmarks.challenge.schema import digest_file, write_table


ROOT_FILES = ("environment.json", "run_config.yaml", "raw_runs.tsv", "summary.tsv", "validation.json", "run.log")
HIGHER_IS_BETTER = (
    "array_recall", "array_precision", "array_f1", "read_detection_recall",
    "read_detection_precision", "sequence_family_recall",
)
LOWER_IS_BETTER = (
    "negative_read_call_rate", "matched_period_mae_bp", "matched_boundary_mae_bp",
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: str | None) -> float | None:
    if value in (None, "", "NA", "nan"):
        return None
    return float(value)


def _copy(source: Path, target: Path, relative: Path) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    expected = digest_file(source)
    if source.stat().st_size != target.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Cascade archive copy differs: {source}")
    return {"file": relative.as_posix(), "source": str(source), "sha256": expected, "bytes": target.stat().st_size}


def validate_run(path: Path, repetitions: int = 3) -> dict[str, object]:
    path = path.resolve()
    missing = [name for name in ROOT_FILES if not (path / name).is_file()]
    if missing:
        raise ValueError(f"Missing cascade evidence in {path}: {', '.join(missing)}")
    environment = json.loads((path / "environment.json").read_text())
    validation = json.loads((path / "validation.json").read_text())
    raw = _read_tsv(path / "raw_runs.tsv")
    summary = _read_tsv(path / "summary.tsv")
    if not summary:
        raise ValueError(f"Empty cascade summary: {path}")
    if validation.get("complete") is not True or int(validation.get("failed_runs", -1)) != 0:
        raise ValueError(f"Incomplete cascade run: {path}")
    expected_runs = repetitions * len(summary)
    if len(raw) != expected_runs or int(validation.get("total_runs", -1)) != expected_runs:
        raise ValueError(f"Cascade raw/summary row-count mismatch: {path}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        grouped[row["scenario"]].append(row)
        if row["status"] != "ok" or row["exit_code"] != "0" or row["timed_out"] != "False":
            raise ValueError(f"Failed cascade repetition: {row['scenario']}")
    if set(grouped) != {row["scenario"] for row in summary}:
        raise ValueError(f"Cascade raw/summary scenarios differ: {path}")
    for scenario, rows in grouped.items():
        if {int(row["repetition"]) for row in rows} != set(range(1, repetitions + 1)):
            raise ValueError(f"Incomplete cascade repetitions: {scenario}")
        for field in ("prediction_sha256", "catalog_sha256"):
            if len({row[field] for row in rows}) != 1:
                raise ValueError(f"Non-deterministic cascade {field}: {scenario}")
    for row in summary:
        if int(row["successful_runs"]) != repetitions or int(row["attempted_runs"]) != repetitions or row["deterministic"] != "True":
            raise ValueError(f"Invalid cascade summary status: {row['scenario']}")
    source_digest = environment.get("source_digest")
    if not isinstance(source_digest, str) or len(source_digest) != 64:
        raise ValueError(f"Missing cascade source digest: {path}")
    return {"path": path, "environment": environment, "raw": raw, "summary": summary}


def compare_to_baseline(baseline: Path, runs: Iterable[dict[str, object]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    baseline_rows = _read_tsv(baseline.resolve())
    baseline_by_scenario = {row["scenario"]: row for row in baseline_rows}
    cascade_rows: list[dict[str, str]] = []
    source_digests: set[str] = set()
    for run in runs:
        cascade_rows.extend(run["summary"])  # type: ignore[arg-type]
        source_digests.add(run["environment"]["source_digest"])  # type: ignore[index]
    cascade_by_scenario = {row["scenario"]: row for row in cascade_rows}
    if len(cascade_by_scenario) != len(cascade_rows):
        raise ValueError("Cascade run parts contain duplicate scenarios")
    if set(cascade_by_scenario) != set(baseline_by_scenario):
        raise ValueError("Cascade and baseline scenario sets differ")
    if len(source_digests) != 1:
        raise ValueError("Cascade run parts use different source snapshots")

    comparison: list[dict[str, object]] = []
    runtime_ratios: list[float] = []
    rss_ratios: list[float] = []
    accuracy_changes = 0
    for scenario in sorted(cascade_by_scenario):
        old = baseline_by_scenario[scenario]
        new = cascade_by_scenario[scenario]
        runtime_ratio = float(new["median_runtime_seconds"]) / float(old["median_runtime_seconds"])
        rss_ratio = float(new["median_peak_rss_mib"]) / float(old["median_peak_rss_mib"])
        runtime_ratios.append(runtime_ratio)
        rss_ratios.append(rss_ratio)
        changes: list[str] = []
        for metric in HIGHER_IS_BETTER + LOWER_IS_BETTER:
            before, after = _number(old.get(metric)), _number(new.get(metric))
            if before is None or after is None:
                if before != after:
                    raise ValueError(f"Cascade metric availability differs: {scenario}/{metric}")
                continue
            delta = after - before
            if metric in HIGHER_IS_BETTER and delta < -1e-12:
                raise ValueError(f"Cascade accuracy regressed: {scenario}/{metric}")
            if metric in LOWER_IS_BETTER and delta > 1e-12:
                raise ValueError(f"Cascade error metric regressed: {scenario}/{metric}")
            if abs(delta) > 1e-12:
                changes.append(f"{metric}:{before:.12g}->{after:.12g}")
                accuracy_changes += 1
        comparison.append({
            "scenario": scenario,
            "baseline_runtime_seconds": float(old["median_runtime_seconds"]),
            "cascade_runtime_seconds": float(new["median_runtime_seconds"]),
            "runtime_ratio": runtime_ratio,
            "baseline_peak_rss_mib": float(old["median_peak_rss_mib"]),
            "cascade_peak_rss_mib": float(new["median_peak_rss_mib"]),
            "peak_rss_ratio": rss_ratio,
            "accuracy_change": ";".join(changes) if changes else "none",
        })
    aggregate = {
        "scenario_count": len(comparison),
        "cascade_source_digest": next(iter(source_digests)),
        "all_common_accuracy_metrics_noninferior": True,
        "changed_accuracy_or_error_cells": accuracy_changes,
        "runtime_geometric_mean_ratio": math.exp(sum(math.log(x) for x in runtime_ratios) / len(runtime_ratios)),
        "runtime_faster_scenarios": sum(x < 1 for x in runtime_ratios),
        "peak_rss_geometric_mean_ratio": math.exp(sum(math.log(x) for x in rss_ratios) / len(rss_ratios)),
        "peak_rss_lower_scenarios": sum(x < 1 for x in rss_ratios),
        "warning": "development_seed_only;baseline_has_one_timing_run_per_scenario;not_independent_heldout_evidence",
    }
    return comparison, aggregate


def archive(baseline: Path, run_paths: list[Path], outdir: Path) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    runs = [validate_run(path) for path in run_paths]
    comparison, aggregate = compare_to_baseline(baseline, runs)
    outdir.mkdir(parents=True)
    combined_raw = [row for run in runs for row in run["raw"]]  # type: ignore[misc]
    combined_summary = [row for run in runs for row in run["summary"]]  # type: ignore[misc]
    write_table(outdir / "combined_raw_runs.tsv", combined_raw, list(combined_raw[0]))
    write_table(outdir / "combined_summary.tsv", combined_summary, list(combined_summary[0]))
    write_table(outdir / "comparison.tsv", comparison, list(comparison[0]))
    (outdir / "aggregate.json").write_text(json.dumps(aggregate, indent=2) + "\n")
    (outdir / "README.md").write_text(
        "# Cascade native-screen development evidence\n\n"
        "This compact archive combines non-overlapping development-seed runs, checks all three "
        "repetitions and output hashes, and compares them with the existing elastic baseline. "
        "It is development evidence; held-out seeds were not used.\n"
    )
    manifest: list[dict[str, object]] = []
    for index, run in enumerate(runs, start=1):
        for name in ROOT_FILES:
            relative = Path(f"run_part_{index}") / name
            manifest.append(_copy(run["path"] / name, outdir / relative, relative))  # type: ignore[operator]
    for name in ("combined_raw_runs.tsv", "combined_summary.tsv", "comparison.tsv", "aggregate.json", "README.md"):
        path = outdir / name
        manifest.append({"file": name, "source": "derived_after_validation", "sha256": digest_file(path), "bytes": path.stat().st_size})
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return aggregate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--run", action="append", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    archive(args.baseline_summary, args.run, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
