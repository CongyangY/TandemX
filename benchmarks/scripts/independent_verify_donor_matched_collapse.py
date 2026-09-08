#!/usr/bin/env python3
"""Independently verify a donor-matched collapse evaluation artifact.

This verifier deliberately uses only the Python standard library and does not
import TandemX's evaluator, BED reader, or copy-number reader.  It recomputes
the family universe, interval unions, fate labels, confusion matrix, confidence
intervals, and missing-bp statistics directly from the frozen inputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, median
from typing import Any, Iterable, Sequence


FAMILY_NUMERIC_FIELDS = (
    "read_estimated_bp",
    "old_assembly_bp",
    "new_assembly_bp",
    "sensitivity_assembly_bp",
    "old_read_ratio",
    "new_read_ratio",
    "old_new_ratio",
    "predicted_missing_bp",
    "observed_gain_bp",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_copy_number(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "family_id" not in reader.fieldnames:
            raise ValueError(f"missing family_id field: {path}")
        bp_field = "estimated_bp" if "estimated_bp" in reader.fieldnames else "estimated_repeat_bp"
        if bp_field not in reader.fieldnames:
            raise ValueError(f"missing estimated bp field: {path}")
        for line_number, row in enumerate(reader, start=2):
            family_id = row.get("family_id", "")
            if not family_id:
                raise ValueError(f"empty family_id at {path}:{line_number}")
            try:
                value = float(row[bp_field])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid estimated bp at {path}:{line_number}") from exc
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"non-finite or negative estimated bp at {path}:{line_number}")
            values[family_id] = values.get(family_id, 0.0) + value
    return values


def read_bed_unions(path: Path) -> dict[str, float]:
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED row has fewer than four fields at {path}:{line_number}")
            chrom, family_id = fields[0], fields[3]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"non-integer BED coordinate at {path}:{line_number}") from exc
            if not chrom or not family_id or start < 0 or end <= start:
                raise ValueError(f"invalid BED interval at {path}:{line_number}")
            intervals.setdefault((family_id, chrom), []).append((start, end))
    totals: dict[str, float] = {}
    for (family_id, _chrom), group in intervals.items():
        merged_bp = 0
        current_start = current_end = None
        for start, end in sorted(group):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= current_end:
                current_end = max(current_end, end)
            else:
                merged_bp += current_end - current_start
                current_start, current_end = start, end
        if current_start is not None:
            merged_bp += current_end - current_start
        totals[family_id] = totals.get(family_id, 0.0) + merged_bp
    return totals


def read_family_metrics(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "family_id", *FAMILY_NUMERIC_FIELDS, "reference_state", "prediction_state",
            "eligibility", "fate", "outcome",
        }
        if not reader.fieldnames or required - set(reader.fieldnames):
            raise ValueError(f"missing family metric fields: {sorted(required - set(reader.fieldnames or []))}")
        rows: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            family_id = row["family_id"]
            if not family_id or family_id in rows:
                raise ValueError(f"empty or duplicate family_id at {path}:{line_number}")
            rows[family_id] = row
    return rows


def ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator > 0 else None


def reference_state(old_bp: float, new_bp: float, collapse: float, overexpansion: float) -> str:
    value = old_bp / new_bp
    if value < collapse:
        return "reference_collapse"
    if value > overexpansion:
        return "old_exceeds_new"
    return "reference_retained"


def prediction_state(old_bp: float, read_bp: float, collapse: float, overexpansion: float) -> str:
    value = old_bp / read_bp
    if value < collapse:
        return "predicted_collapse"
    if value > overexpansion:
        return "predicted_old_excess"
    return "predicted_retained"


def expected_family_rows(
    read_bp: dict[str, float],
    old_bp: dict[str, float],
    new_bp: dict[str, float],
    sensitivity_bp: dict[str, float],
    *,
    collapse: float,
    overexpansion: float,
    primary_min_new_bp: int,
) -> dict[str, dict[str, Any]]:
    expected: dict[str, dict[str, Any]] = {}
    for family_id in sorted(set(read_bp) | set(old_bp) | set(new_bp) | set(sensitivity_bp)):
        read_value = read_bp.get(family_id, 0.0)
        old_value = old_bp.get(family_id, 0.0)
        new_value = new_bp.get(family_id, 0.0)
        sensitivity_value = sensitivity_bp.get(family_id, 0.0)
        ref_state = (
            reference_state(old_value, new_value, collapse, overexpansion)
            if new_value > 0 else "unresolved"
        )
        pred_state = (
            prediction_state(old_value, read_value, collapse, overexpansion)
            if read_value > 0 else "unresolved"
        )
        eligibility, fate, outcome = "eligible", "evaluated", "NA"
        if new_value < primary_min_new_bp:
            eligibility, fate = "not_source_eligible", "new_reference_below_primary_min_bp"
        elif read_value <= 0:
            eligibility, fate = "technical_failure", "missing_positive_read_estimate"
        else:
            truth_positive = ref_state == "reference_collapse"
            predicted_positive = pred_state == "predicted_collapse"
            outcome = (
                "TP" if truth_positive and predicted_positive else
                "FN" if truth_positive else
                "FP" if predicted_positive else "TN"
            )
        expected[family_id] = {
            "family_id": family_id,
            "read_estimated_bp": read_value,
            "old_assembly_bp": old_value,
            "new_assembly_bp": new_value,
            "sensitivity_assembly_bp": sensitivity_value,
            "old_read_ratio": ratio(old_value, read_value),
            "new_read_ratio": ratio(new_value, read_value),
            "old_new_ratio": ratio(old_value, new_value),
            "predicted_missing_bp": max(read_value - old_value, 0.0),
            "observed_gain_bp": max(new_value - old_value, 0.0),
            "reference_state": ref_state,
            "prediction_state": pred_state,
            "eligibility": eligibility,
            "fate": fate,
            "outcome": outcome,
        }
    return expected


def ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and values[order[end]] == values[order[index]]:
            end += 1
        rank = (index + 1 + end) / 2
        for original_index in order[index:end]:
            result[original_index] = rank
        index = end
    return result


def pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) < 2 or len(left) != len(right):
        return None
    left_mean, right_mean = mean(left), mean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    denominator = math.sqrt(
        sum((x - left_mean) ** 2 for x in left) * sum((y - right_mean) ** 2 for y in right)
    )
    return numerator / denominator if denominator else None


def divide(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> list[float] | None:
    if total <= 0:
        return None
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    half = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


def expected_summary(
    rows: Iterable[dict[str, Any]],
    *,
    collapse: float,
    overexpansion: float,
    primary_min_new_bp: int,
    sensitivity_thresholds: Sequence[int],
) -> dict[str, Any]:
    all_rows = list(rows)
    primary = [row for row in all_rows if row["eligibility"] == "eligible"]
    confusion = {label: sum(row["outcome"] == label for row in primary) for label in ("TP", "FN", "FP", "TN")}
    tp, fn, fp, tn = (confusion[label] for label in ("TP", "FN", "FP", "TN"))
    sensitivity = divide(tp, tp + fn)
    specificity = divide(tn, tn + fp)
    predicted = [float(row["predicted_missing_bp"]) for row in primary]
    observed = [float(row["observed_gain_bp"]) for row in primary]
    errors = [abs(left - right) for left, right in zip(predicted, observed)]
    threshold_rows = []
    for threshold in sensitivity_thresholds:
        selected = [
            row for row in all_rows
            if float(row["new_assembly_bp"]) >= threshold and float(row["read_estimated_bp"]) > 0
        ]
        counts = {label: 0 for label in ("TP", "FN", "FP", "TN")}
        for row in selected:
            truth_positive = row["reference_state"] == "reference_collapse"
            predicted_positive = row["prediction_state"] == "predicted_collapse"
            label = (
                "TP" if truth_positive and predicted_positive else
                "FN" if truth_positive else
                "FP" if predicted_positive else "TN"
            )
            counts[label] += 1
        threshold_rows.append({"min_new_bp": threshold, "families": len(selected), **counts})
    mcc_denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "interpretation": "donor_matched_retrospective_reference_proxy_not_absolute_biological_truth",
        "eligibility_rule": "new_assembly_bp_at_least_primary_min_new_bp_independent_of_read_estimate",
        "config": {
            "collapse_threshold": collapse,
            "overexpansion_threshold": overexpansion,
            "primary_min_new_bp": primary_min_new_bp,
            "min_new_bp_sensitivity": list(sensitivity_thresholds),
        },
        "all_family_rows": len(all_rows),
        "eligible_family_rows": len(primary),
        "not_source_eligible_rows": sum(row["eligibility"] == "not_source_eligible" for row in all_rows),
        "technical_failure_rows": sum(row["eligibility"] == "technical_failure" for row in all_rows),
        "confusion": confusion,
        "sensitivity": sensitivity,
        "sensitivity_wilson95": wilson(tp, tp + fn),
        "false_positive_rate": divide(fp, fp + tn),
        "false_positive_rate_wilson95": wilson(fp, fp + tn),
        "precision": divide(tp, tp + fp),
        "precision_wilson95": wilson(tp, tp + fp),
        "specificity": specificity,
        "balanced_accuracy": mean([sensitivity, specificity]) if sensitivity is not None and specificity is not None else None,
        "matthews_correlation": (tp * tn - fp * fn) / mcc_denominator if mcc_denominator else None,
        "predicted_missing_bp_total": sum(predicted),
        "observed_gain_bp_total": sum(observed),
        "missing_bp_mean_absolute_error": mean(errors) if errors else None,
        "missing_bp_median_absolute_error": median(errors) if errors else None,
        "missing_bp_pearson": pearson(predicted, observed),
        "missing_bp_spearman": pearson(ranks(predicted), ranks(observed)) if predicted else None,
        "min_new_bp_sensitivity": threshold_rows,
    }


def close(left: Any, right: Any, tolerance: float = 1e-9) -> bool:
    if left is None or right is None:
        return left is None and right is None
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(close(x, y, tolerance) for x, y in zip(left, right))
    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(close(left[key], right[key], tolerance) for key in left)
    return left == right


def parsed_numeric(value: str) -> float | None:
    return None if value == "NA" else float(value)


def verify(
    *,
    config_path: Path,
    copy_number_path: Path,
    old_arrays_path: Path,
    new_arrays_path: Path,
    sensitivity_arrays_path: Path,
    family_metrics_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    config_payload = json.loads(config_path.read_text(encoding="utf-8"))
    config = config_payload["evaluation"]
    collapse = float(config["collapse_threshold"])
    overexpansion = float(config["overexpansion_threshold"])
    primary_min_new_bp = int(config["primary_min_new_bp"])
    sensitivity_thresholds = [int(value) for value in config["min_new_bp_sensitivity"]]
    expected_rows = expected_family_rows(
        read_copy_number(copy_number_path),
        read_bed_unions(old_arrays_path),
        read_bed_unions(new_arrays_path),
        read_bed_unions(sensitivity_arrays_path),
        collapse=collapse,
        overexpansion=overexpansion,
        primary_min_new_bp=primary_min_new_bp,
    )
    observed_rows = read_family_metrics(family_metrics_path)
    failures: list[str] = []
    if set(expected_rows) != set(observed_rows):
        failures.append("family_universe_mismatch")
    for family_id in sorted(set(expected_rows) & set(observed_rows)):
        expected, observed = expected_rows[family_id], observed_rows[family_id]
        for field in FAMILY_NUMERIC_FIELDS:
            if not close(expected[field], parsed_numeric(observed[field])):
                failures.append(f"family_numeric_mismatch:{family_id}:{field}")
        for field in ("reference_state", "prediction_state", "eligibility", "fate", "outcome"):
            if expected[field] != observed[field]:
                failures.append(f"family_label_mismatch:{family_id}:{field}")
    expected_summary_payload = expected_summary(
        expected_rows.values(),
        collapse=collapse,
        overexpansion=overexpansion,
        primary_min_new_bp=primary_min_new_bp,
        sensitivity_thresholds=sensitivity_thresholds,
    )
    observed_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    for key, value in expected_summary_payload.items():
        if key not in observed_summary or not close(value, observed_summary[key]):
            failures.append(f"summary_mismatch:{key}")
    paths = {
        "config": config_path,
        "copy_number": copy_number_path,
        "old_arrays": old_arrays_path,
        "new_arrays": new_arrays_path,
        "sensitivity_arrays": sensitivity_arrays_path,
        "family_metrics": family_metrics_path,
        "summary": summary_path,
    }
    return {
        "schema_version": 1,
        "verification_method": "independent_standard_library_recomputation_without_tandemx_evaluator_imports",
        "verification_passed": not failures,
        "family_rows_recomputed": len(expected_rows),
        "checks": {
            "family_universe_match": "family_universe_mismatch" not in failures,
            "family_values_and_labels_match": not any(item.startswith("family_") and item != "family_universe_mismatch" for item in failures),
            "summary_statistics_match": not any(item.startswith("summary_mismatch:") for item in failures),
        },
        "failures": failures,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--copy-number", required=True, type=Path)
    parser.add_argument("--old-arrays", required=True, type=Path)
    parser.add_argument("--new-arrays", required=True, type=Path)
    parser.add_argument("--sensitivity-arrays", required=True, type=Path)
    parser.add_argument("--family-metrics", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(f"refusing to overwrite receipt: {args.receipt}")
    result = verify(
        config_path=args.config,
        copy_number_path=args.copy_number,
        old_arrays_path=args.old_arrays,
        new_arrays_path=args.new_arrays,
        sensitivity_arrays_path=args.sensitivity_arrays,
        family_metrics_path=args.family_metrics,
        summary_path=args.summary,
    )
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if result["verification_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
