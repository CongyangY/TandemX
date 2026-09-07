"""Compact and independently summarize cascade gap-free development evidence."""
from __future__ import annotations

import argparse
import csv
import json
import math
import pstats
import shutil
from pathlib import Path
from typing import Any

from benchmarks.challenge.schema import digest_file, write_table


ROOT_FILES = (
    "environment.json",
    "run_config.yaml",
    "validation.json",
    "raw_runs.tsv",
    "summary.tsv",
    "run.log",
)
SCIENTIFIC_FIELDS = (
    "array_recall",
    "array_precision",
    "array_f1",
    "read_detection_recall",
    "read_detection_precision",
    "negative_read_call_rate",
    "sequence_family_recall",
    "matched_period_mae_bp",
    "matched_boundary_mae_bp",
    "cyclic_monomer_recall",
    "homologous_consensus_fraction",
    "mean_best_cyclic_edit_similarity",
    "base_union_recall",
    "base_union_precision",
    "base_union_f1",
    "duplicate_bp_fraction",
)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _copy(source: Path, destination: Path, relative: Path) -> dict[str, object]:
    if not source.is_file():
        raise ValueError(f"Required source is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    checksum = digest_file(source)
    if digest_file(destination) != checksum or destination.stat().st_size != source.stat().st_size:
        raise OSError(f"Copied evidence differs from source: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source),
        "sha256": checksum,
        "bytes": destination.stat().st_size,
    }


def runtime_geomean_ratio(summary: list[dict[str, str]]) -> float:
    by_key = {(row["scenario"], row["tool"]): row for row in summary}
    scenarios = sorted({row["scenario"] for row in summary})
    ratios = [
        float(by_key[(scenario, "tandemx")]["median_runtime_seconds"])
        / float(by_key[(scenario, "tidehunter")]["median_runtime_seconds"])
        for scenario in scenarios
    ]
    if not ratios or any(value <= 0 for value in ratios):
        raise ValueError("Runtime ratios must be complete and positive")
    return math.exp(sum(math.log(value) for value in ratios) / len(ratios))


def _validate_run(run: Path, config_hash: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    missing = [name for name in ROOT_FILES if not (run / name).is_file()]
    if missing:
        raise ValueError(f"Development run is incomplete: {', '.join(missing)}")
    validation = json.loads((run / "validation.json").read_text())
    if validation != {
        "complete": True,
        "successful_runs": 96,
        "failed_runs": 0,
        "total_runs": 96,
    }:
        raise ValueError("Development run receipt is not the complete 96-run matrix")
    environment = json.loads((run / "environment.json").read_text())
    if environment.get("config_sha256") != config_hash or environment.get("split") != "development":
        raise ValueError("Development run config hash or split differs")
    raw = _rows(run / "raw_runs.tsv")
    summary = _rows(run / "summary.tsv")
    if len(raw) != 96 or len(summary) != 32 or any(row["status"] != "ok" for row in raw):
        raise ValueError("Development matrix rows are incomplete or failed")
    return raw, summary


def _dataset_hashes(run: Path) -> dict[str, str]:
    return {
        path.parent.name: digest_file(path)
        for path in sorted((run / "datasets").glob("*/manifest.json"))
    }


def _profile_rows(profiles: Path) -> tuple[list[dict[str, object]], dict[str, str]]:
    rows: list[dict[str, object]] = []
    hashes: dict[str, str] = {}
    for path in sorted(profiles.glob("*.prof")):
        hashes[path.name] = digest_file(path)
        stats = pstats.Stats(str(path))
        ranked = sorted(stats.stats.items(), key=lambda item: item[1][3], reverse=True)
        rank = 0
        for (filename, line, function), (primitive, calls, self_time, cumulative, _callers) in ranked:
            if "/tandemx/" not in filename and "tandemx._rust_core" not in function:
                continue
            rank += 1
            rows.append(
                {
                    "profile": path.stem,
                    "rank": rank,
                    "filename": filename,
                    "line": line,
                    "function": function,
                    "primitive_calls": primitive,
                    "total_calls": calls,
                    "self_seconds": self_time,
                    "cumulative_seconds": cumulative,
                }
            )
            if rank == 20:
                break
    if not rows:
        raise ValueError("No TandemX profile rows were found")
    return rows, hashes


def archive(
    baseline: Path,
    rejected: Path,
    candidate: Path,
    audit: Path,
    profiles: Path,
    config: Path,
    outdir: Path,
) -> dict[str, Any]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    config_hash = digest_file(config)
    runs = {"baseline": baseline, "rejected_intermediate": rejected, "candidate": candidate}
    parsed = {name: _validate_run(path, config_hash) for name, path in runs.items()}
    dataset_hashes = {name: _dataset_hashes(path) for name, path in runs.items()}
    if len(dataset_hashes["baseline"]) != 16 or len({json.dumps(value, sort_keys=True) for value in dataset_hashes.values()}) != 1:
        raise ValueError("Development runs do not use identical 16-dataset manifests")

    outdir.mkdir(parents=True)
    manifest: list[dict[str, object]] = []
    for label, run in runs.items():
        for name in ROOT_FILES:
            relative = Path(label) / name
            manifest.append(_copy(run / name, outdir / relative, relative))
    manifest.append(_copy(config, outdir / "development_config.yaml", Path("development_config.yaml")))
    for name in ("receipt.json", "screen_features.tsv", "threshold_grid.tsv", "scenario_summary.tsv"):
        relative = Path("fast_path_audit") / name
        manifest.append(_copy(audit / name, outdir / relative, relative))

    summaries = {name: value[1] for name, value in parsed.items()}
    ratios = {name: runtime_geomean_ratio(rows) for name, rows in summaries.items()}
    baseline_by = {(row["scenario"], row["tool"]): row for row in summaries["baseline"]}
    candidate_by = {(row["scenario"], row["tool"]): row for row in summaries["candidate"]}
    comparison_rows: list[dict[str, object]] = []
    metric_differences: list[dict[str, object]] = []
    for key in sorted(baseline_by):
        old, new = baseline_by[key], candidate_by[key]
        if key[1] == "tandemx":
            comparison_rows.append(
                {
                    "scenario": key[0],
                    "baseline_runtime_seconds": old["median_runtime_seconds"],
                    "candidate_runtime_seconds": new["median_runtime_seconds"],
                    "speedup": float(old["median_runtime_seconds"])
                    / float(new["median_runtime_seconds"]),
                    "baseline_peak_rss_mib": old["median_peak_rss_mib"],
                    "candidate_peak_rss_mib": new["median_peak_rss_mib"],
                }
            )
        for field in SCIENTIFIC_FIELDS:
            if old[field] != new[field]:
                metric_differences.append(
                    {
                        "scenario": key[0],
                        "tool": key[1],
                        "metric": field,
                        "baseline": old[field],
                        "candidate": new[field],
                    }
                )

    branch_rows: list[dict[str, object]] = []
    for scenario in sorted({row["scenario"] for row in summaries["candidate"]}):
        path = candidate / "runs" / f"{scenario}_s1201" / "tandemx" / "rep1" / "discover" / "candidate_reads.tsv"
        counts = {"cascade_gap_free": 0, "cascade_elastic": 0}
        for row in _rows(path):
            for branch in counts:
                if branch in row["warning"].split(";"):
                    counts[branch] += 1
        branch_rows.append({"scenario": scenario, **counts})

    profile_rows, profile_hashes = _profile_rows(profiles)
    write_table(outdir / "scenario_comparison.tsv", comparison_rows, list(comparison_rows[0]))
    write_table(outdir / "metric_differences.tsv", metric_differences, list(metric_differences[0]))
    write_table(outdir / "branch_counts.tsv", branch_rows, list(branch_rows[0]))
    write_table(outdir / "profile_top.tsv", profile_rows, list(profile_rows[0]))
    decision = {
        "status": "development_passed",
        "runtime_geomean_ratio": ratios,
        "candidate_gate": "tandemx_divided_by_tidehunter_runtime_geomean_ratio_le_2",
        "candidate_gate_passed": ratios["candidate"] <= 2.0,
        "candidate_metric_difference_rows": len(metric_differences),
        "profile_sha256": profile_hashes,
        "interpretation": (
            "the guarded 30%-span gap-free path passed the development runtime gate; "
            "the unguarded intermediate indel boundary regression is retained"
        ),
        "warning": "development_seed_1201_used_for_selection;not_validation;timing_repetitions_not_biological_replicates",
    }
    (outdir / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    (outdir / "README.md").write_text(
        "# Cascade gap-free development evidence\n\n"
        "Seed 1201 is development data. The baseline, rejected intermediate and "
        "guarded candidate each completed 96 runs on identical inputs. The "
        "intermediate reached a 1.933 runtime ratio but worsened 0.1% indel "
        "boundary MAE to 14.221 bp. Requiring >=95% valid shifted columns and "
        "<=2% unit-span residual retained 403 truth-scored safe development "
        "acceptances and restored boundary metrics. The final candidate ratio "
        "was 1.989. Validation seed 2201 must be run once only after the frozen "
        "source and configuration pass hosted CI.\n"
    )
    for name in (
        "scenario_comparison.tsv",
        "metric_differences.tsv",
        "branch_counts.tsv",
        "profile_top.tsv",
        "decision.json",
        "README.md",
    ):
        path = outdir / name
        manifest.append(
            {
                "file": name,
                "source": "derived_from_validated_development_inputs",
                "sha256": digest_file(path),
                "bytes": path.stat().st_size,
            }
        )
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--rejected", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--profiles", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(archive(args.baseline, args.rejected, args.candidate, args.audit, args.profiles, args.config, args.outdir), indent=2))


if __name__ == "__main__":
    main()
