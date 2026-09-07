"""Audit observable cascade-screen thresholds on development datasets.

This analysis uses planted truth only to score candidate rules. The proposed
runtime rules themselves use screen interval length, shifted identity and local
composition; they never receive truth labels.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Iterable

from benchmarks.challenge.schema import digest_file, read_table, write_table
from tandemx.discover.rust_backend import scan_read_for_periods
from tandemx.io.sequences import read_sequence_records


FEATURE_FIELDS = [
    "dataset_id",
    "scenario",
    "seed",
    "read_id",
    "read_length",
    "truth_array_count",
    "screen_status",
    "candidate_period_count",
    "best_period",
    "periodicity_score",
    "repeat_start",
    "repeat_end",
    "repeat_span",
    "repeat_span_fraction",
    "unit_span_residual_fraction",
    "shifted_identity",
    "valid_pair_fraction",
    "composition_adjusted_identity",
    "maximum_truth_iou",
    "minimum_truth_period_error_bp",
]


def parse_grid(text: str) -> list[float]:
    values = [float(value) for value in text.split(",")]
    if not values or any(not math.isfinite(value) or value < 0 or value > 1 for value in values):
        raise ValueError("threshold grids must contain finite values in [0,1]")
    return values


def interval_iou(left: tuple[int, int], right: tuple[int, int]) -> float:
    intersection = max(0, min(left[1], right[1]) - max(left[0], right[0]))
    union = max(left[1], right[1]) - min(left[0], right[0])
    return intersection / union if union else 0.0


def shifted_evidence(
    sequence: str, period: int, start: int, end: int
) -> tuple[float, float, float]:
    start = max(0, start)
    end = min(len(sequence), end)
    if period <= 0 or end - start <= period:
        return 0.0, 0.0, 0.0
    pairs = [
        (left, left + period)
        for left in range(start, end - period)
        if sequence[left] in "ACGT" and sequence[left + period] in "ACGT"
    ]
    if not pairs:
        return 0.0, 0.0, 0.0
    valid_fraction = len(pairs) / (end - start - period)
    identity = sum(sequence[left] == sequence[right] for left, right in pairs) / len(pairs)
    local = sequence[start:end]
    counts = {base: local.count(base) for base in "ACGT"}
    total = sum(counts.values())
    if total == 0:
        return identity, valid_fraction, 0.0
    chance = sum((count / total) ** 2 for count in counts.values())
    adjusted = (identity - chance) / (1 - chance) if chance < 1 else 0.0
    return identity, valid_fraction, adjusted


def feature_rows(dataset_root: Path) -> tuple[list[dict[str, object]], dict[str, str]]:
    rows: list[dict[str, object]] = []
    manifests: dict[str, str] = {}
    datasets = sorted(path for path in dataset_root.iterdir() if path.is_dir())
    if not datasets:
        raise ValueError(f"No datasets found in {dataset_root}")
    for dataset in datasets:
        manifest_path = dataset / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifests[dataset.name] = digest_file(manifest_path)
        scenario = manifest["scenario"]
        truth_by_read: dict[str, list[dict[str, str]]] = {}
        for truth in read_table(dataset / "truth_arrays.tsv"):
            truth_by_read.setdefault(truth["read_id"], []).append(truth)
        for record in read_sequence_records(dataset / "reads.fa"):
            screen = scan_read_for_periods(
                record.sequence,
                k=11,
                min_period=30,
                max_period=1000,
                top_periods=5,
                min_seed_occurrences=2,
                min_spacing_support=2,
                max_pairs_per_kmer=100,
                min_repeat_span=100,
            )
            identity, valid_fraction, adjusted = shifted_evidence(
                record.sequence,
                screen.best_period,
                screen.repeat_start,
                screen.repeat_end,
            )
            truths = truth_by_read.get(record.id, [])
            ious = [
                interval_iou(
                    (screen.repeat_start, screen.repeat_end),
                    (int(truth["start"]), int(truth["end"])),
                )
                for truth in truths
            ]
            period_errors = [abs(screen.best_period - int(truth["period"])) for truth in truths]
            span = max(0, screen.repeat_end - screen.repeat_start)
            remainder = span % screen.best_period if screen.best_period else 0
            residual = (
                min(remainder, screen.best_period - remainder) / screen.best_period
                if screen.best_period
                else 1.0
            )
            rows.append(
                {
                    "dataset_id": dataset.name,
                    "scenario": scenario["name"],
                    "seed": manifest["seed"],
                    "read_id": record.id,
                    "read_length": len(record.sequence),
                    "truth_array_count": len(truths),
                    "screen_status": screen.status,
                    "candidate_period_count": len(screen.candidate_periods),
                    "best_period": screen.best_period,
                    "periodicity_score": screen.periodicity_score,
                    "repeat_start": screen.repeat_start,
                    "repeat_end": screen.repeat_end,
                    "repeat_span": span,
                    "repeat_span_fraction": span / len(record.sequence),
                    "unit_span_residual_fraction": residual,
                    "shifted_identity": identity,
                    "valid_pair_fraction": valid_fraction,
                    "composition_adjusted_identity": adjusted,
                    "maximum_truth_iou": max(ious, default=0.0),
                    "minimum_truth_period_error_bp": min(period_errors, default="NA"),
                }
            )
    return rows, manifests


def threshold_rows(
    features: Iterable[dict[str, object]],
    span_fractions: list[float],
    identities: list[float],
    *,
    minimum_span: int = 100,
    minimum_adjusted_identity: float = 0.7,
    minimum_valid_pair_fraction: float = 0.95,
    maximum_unit_span_residual_fraction: float = 0.02,
) -> list[dict[str, object]]:
    features = list(features)
    result: list[dict[str, object]] = []
    single_positive = sum(int(row["truth_array_count"]) == 1 for row in features)
    for span_fraction in span_fractions:
        for identity_threshold in identities:
            accepted = [
                row
                for row in features
                if int(row["best_period"]) > 0
                and int(row["repeat_span"])
                >= max(minimum_span, math.ceil(span_fraction * int(row["read_length"])))
                and float(row["shifted_identity"]) >= identity_threshold
                and float(row["valid_pair_fraction"])
                >= minimum_valid_pair_fraction
                and float(row["composition_adjusted_identity"])
                >= minimum_adjusted_identity
                and float(row["unit_span_residual_fraction"])
                <= maximum_unit_span_residual_fraction
            ]
            correct_single = sum(
                int(row["truth_array_count"]) == 1
                and float(row["maximum_truth_iou"]) >= 0.5
                and int(row["minimum_truth_period_error_bp"]) <= max(3, round(0.05 * int(row["best_period"])))
                for row in accepted
            )
            negative = sum(int(row["truth_array_count"]) == 0 for row in accepted)
            multi = sum(int(row["truth_array_count"]) > 1 for row in accepted)
            wrong_single = sum(int(row["truth_array_count"]) == 1 for row in accepted) - correct_single
            result.append(
                {
                    "minimum_span_fraction": span_fraction,
                    "minimum_shifted_identity": identity_threshold,
                    "minimum_composition_adjusted_identity": minimum_adjusted_identity,
                    "minimum_valid_pair_fraction": minimum_valid_pair_fraction,
                    "maximum_unit_span_residual_fraction": maximum_unit_span_residual_fraction,
                    "accepted_reads": len(accepted),
                    "accepted_fraction": len(accepted) / len(features),
                    "correct_single_array_acceptances": correct_single,
                    "correct_single_array_fraction": correct_single / single_positive,
                    "negative_acceptances": negative,
                    "wrong_single_array_acceptances": wrong_single,
                    "multi_array_acceptances": multi,
                    "development_safe": negative == 0 and wrong_single == 0 and multi == 0,
                }
            )
    return result


def scenario_rows(features: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in features:
        grouped.setdefault(str(row["scenario"]), []).append(row)
    result = []
    for scenario, rows in sorted(grouped.items()):
        positives = [row for row in rows if int(row["truth_array_count"]) > 0]
        result.append(
            {
                "scenario": scenario,
                "reads": len(rows),
                "truth_positive_reads": len(positives),
                "median_screen_span_fraction_positive": statistics.median(
                    float(row["repeat_span_fraction"]) for row in positives
                )
                if positives
                else "NA",
                "median_shifted_identity_positive": statistics.median(
                    float(row["shifted_identity"]) for row in positives
                )
                if positives
                else "NA",
                "maximum_adjusted_identity_negative": max(
                    (
                        float(row["composition_adjusted_identity"])
                        for row in rows
                        if int(row["truth_array_count"]) == 0
                    ),
                    default="NA",
                ),
            }
        )
    return result


def run(
    dataset_root: Path,
    outdir: Path,
    span_fractions: list[float],
    identities: list[float],
) -> dict[str, object]:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"Choose a new empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    features, manifests = feature_rows(dataset_root)
    thresholds = threshold_rows(features, span_fractions, identities)
    scenarios = scenario_rows(features)
    write_table(outdir / "screen_features.tsv", features, FEATURE_FIELDS)
    write_table(outdir / "threshold_grid.tsv", thresholds, list(thresholds[0]))
    write_table(outdir / "scenario_summary.tsv", scenarios, list(scenarios[0]))
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    receipt = {
        "complete": True,
        "dataset_root": str(dataset_root.resolve()),
        "dataset_manifests": manifests,
        "git_head": revision,
        "feature_rows": len(features),
        "threshold_rows": len(thresholds),
        "scenario_rows": len(scenarios),
        "span_fractions": span_fractions,
        "identity_thresholds": identities,
        "minimum_composition_adjusted_identity": 0.7,
        "minimum_valid_pair_fraction": 0.95,
        "maximum_unit_span_residual_fraction": 0.02,
        "warning": "development_truth_used_for_scoring;rules_use_observed_screen_features_only;not_validation",
    }
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--span-fractions", default="0,0.1,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--identity-thresholds", default="0.95,0.975,0.99,0.995")
    args = parser.parse_args()
    receipt = run(
        args.dataset_root,
        args.outdir,
        parse_grid(args.span_fractions),
        parse_grid(args.identity_thresholds),
    )
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
