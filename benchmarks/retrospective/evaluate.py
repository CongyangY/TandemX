"""Evaluate frozen read-based collapse predictions against a newer assembly.

The newer assembly is treated as a retrospective reference proxy, not absolute
biological truth.  Source eligibility is deliberately independent of the read
estimate under evaluation: a family enters the primary denominator solely when
the newer assembly contains at least the predeclared amount of localized array
sequence.  Every discovered or localized family remains in the output table.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from statistics import mean, median
from typing import Iterable, Sequence

from benchmarks.challenge.schema import write_table
from tandemx.compare.mvp import read_arrays_bed, read_copy_number, union_interval_length


FAMILY_FIELDS = [
    "family_id",
    "read_estimated_bp",
    "old_assembly_bp",
    "new_assembly_bp",
    "sensitivity_assembly_bp",
    "old_read_ratio",
    "new_read_ratio",
    "old_new_ratio",
    "predicted_missing_bp",
    "observed_gain_bp",
    "reference_state",
    "prediction_state",
    "eligibility",
    "fate",
    "outcome",
]


@dataclass(frozen=True)
class EvaluationConfig:
    collapse_threshold: float = 0.6
    overexpansion_threshold: float = 1.5
    primary_min_new_bp: int = 15_000
    min_new_bp_sensitivity: tuple[int, ...] = (5_000, 15_000, 50_000)

    def validate(self) -> None:
        if not 0 < self.collapse_threshold < 1:
            raise ValueError("collapse_threshold must be in (0, 1)")
        if self.overexpansion_threshold <= 1:
            raise ValueError("overexpansion_threshold must be greater than 1")
        if self.primary_min_new_bp <= 0:
            raise ValueError("primary_min_new_bp must be positive")
        if not self.min_new_bp_sensitivity or any(value <= 0 for value in self.min_new_bp_sensitivity):
            raise ValueError("min_new_bp_sensitivity must contain positive values")
        if self.primary_min_new_bp not in self.min_new_bp_sensitivity:
            raise ValueError("primary_min_new_bp must be included in min_new_bp_sensitivity")


def assembly_bp_by_family(path: Path) -> dict[str, float]:
    """Return chromosome-aware interval-union bp for each family."""
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for row in read_arrays_bed(path):
        grouped.setdefault((row.family_id, row.chrom), []).append((row.start, row.end))
    result: dict[str, float] = {}
    for (family_id, _chrom), intervals in grouped.items():
        result[family_id] = result.get(family_id, 0.0) + union_interval_length(intervals)
    return result


def _ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator > 0 else None


def _reference_state(old_bp: float, new_bp: float, config: EvaluationConfig) -> str:
    ratio = old_bp / new_bp
    if ratio < config.collapse_threshold:
        return "reference_collapse"
    if ratio > config.overexpansion_threshold:
        return "old_exceeds_new"
    return "reference_retained"


def _prediction_state(old_bp: float, read_bp: float, config: EvaluationConfig) -> str:
    ratio = old_bp / read_bp
    if ratio < config.collapse_threshold:
        return "predicted_collapse"
    if ratio > config.overexpansion_threshold:
        return "predicted_old_excess"
    return "predicted_retained"


def evaluate_families(
    read_bp: dict[str, float],
    old_bp: dict[str, float],
    new_bp: dict[str, float],
    sensitivity_bp: dict[str, float] | None,
    config: EvaluationConfig,
) -> list[dict[str, object]]:
    """Evaluate every family without dropping adverse or ineligible rows."""
    config.validate()
    sensitivity_bp = sensitivity_bp or {}
    family_ids = sorted(set(read_bp) | set(old_bp) | set(new_bp) | set(sensitivity_bp))
    rows: list[dict[str, object]] = []
    for family_id in family_ids:
        read_value = read_bp.get(family_id, 0.0)
        old_value = old_bp.get(family_id, 0.0)
        new_value = new_bp.get(family_id, 0.0)
        sensitivity_value = sensitivity_bp.get(family_id, 0.0)
        reference_state = _reference_state(old_value, new_value, config) if new_value > 0 else "unresolved"
        prediction_state = _prediction_state(old_value, read_value, config) if read_value > 0 else "unresolved"
        eligibility = "eligible"
        fate = "evaluated"
        outcome = "NA"
        if new_value < config.primary_min_new_bp:
            eligibility = "not_source_eligible"
            fate = "new_reference_below_primary_min_bp"
        else:
            if read_value <= 0:
                eligibility = "technical_failure"
                fate = "missing_positive_read_estimate"
            else:
                truth_positive = reference_state == "reference_collapse"
                predicted_positive = prediction_state == "predicted_collapse"
                outcome = (
                    "TP" if truth_positive and predicted_positive else
                    "FN" if truth_positive else
                    "FP" if predicted_positive else
                    "TN"
                )
        rows.append(
            {
                "family_id": family_id,
                "read_estimated_bp": read_value,
                "old_assembly_bp": old_value,
                "new_assembly_bp": new_value,
                "sensitivity_assembly_bp": sensitivity_value,
                "old_read_ratio": _ratio(old_value, read_value),
                "new_read_ratio": _ratio(new_value, read_value),
                "old_new_ratio": _ratio(old_value, new_value),
                "predicted_missing_bp": max(read_value - old_value, 0.0),
                "observed_gain_bp": max(new_value - old_value, 0.0),
                "reference_state": reference_state,
                "prediction_state": prediction_state,
                "eligibility": eligibility,
                "fate": fate,
                "outcome": outcome,
            }
        )
    return rows


def _rank(values: Sequence[float]) -> list[float]:
    ordered = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[index]]:
            end += 1
        average_rank = (index + 1 + end) / 2
        for ordered_index in ordered[index:end]:
            ranks[ordered_index] = average_rank
        index = end
    return ranks


def _pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) < 2 or len(left) != len(right):
        return None
    left_mean = mean(left)
    right_mean = mean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    denominator = math.sqrt(left_ss * right_ss)
    return numerator / denominator if denominator else None


def _safe_divide(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def _wilson(successes: int, total: int, z: float = 1.959963984540054) -> list[float] | None:
    if total <= 0:
        return None
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    half = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


def summarize_rows(rows: Sequence[dict[str, object]], config: EvaluationConfig) -> dict[str, object]:
    config.validate()
    primary = [row for row in rows if row["eligibility"] == "eligible"]
    confusion = {name: sum(row["outcome"] == name for row in primary) for name in ("TP", "FN", "FP", "TN")}
    tp, fn, fp, tn = (confusion[name] for name in ("TP", "FN", "FP", "TN"))
    mcc_denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    sensitivity = _safe_divide(tp, tp + fn)
    specificity = _safe_divide(tn, tn + fp)
    predicted = [float(row["predicted_missing_bp"]) for row in primary]
    observed = [float(row["observed_gain_bp"]) for row in primary]
    absolute_errors = [abs(left - right) for left, right in zip(predicted, observed)]
    sensitivity_rows = []
    for threshold in config.min_new_bp_sensitivity:
        selected = [row for row in rows if float(row["new_assembly_bp"]) >= threshold and float(row["read_estimated_bp"]) > 0]
        counts = {name: 0 for name in ("TP", "FN", "FP", "TN")}
        for row in selected:
            truth_positive = row["reference_state"] == "reference_collapse"
            predicted_positive = row["prediction_state"] == "predicted_collapse"
            outcome = (
                "TP" if truth_positive and predicted_positive else
                "FN" if truth_positive else
                "FP" if predicted_positive else
                "TN"
            )
            counts[outcome] += 1
        sensitivity_rows.append({"min_new_bp": threshold, "families": len(selected), **counts})
    return {
        "schema_version": 1,
        "interpretation": "donor_matched_retrospective_reference_proxy_not_absolute_biological_truth",
        "eligibility_rule": "new_assembly_bp_at_least_primary_min_new_bp_independent_of_read_estimate",
        "config": asdict(config),
        "all_family_rows": len(rows),
        "eligible_family_rows": len(primary),
        "not_source_eligible_rows": sum(row["eligibility"] == "not_source_eligible" for row in rows),
        "technical_failure_rows": sum(row["eligibility"] == "technical_failure" for row in rows),
        "confusion": confusion,
        "sensitivity": sensitivity,
        "sensitivity_wilson95": _wilson(tp, tp + fn),
        "false_positive_rate": _safe_divide(fp, fp + tn),
        "false_positive_rate_wilson95": _wilson(fp, fp + tn),
        "precision": _safe_divide(tp, tp + fp),
        "precision_wilson95": _wilson(tp, tp + fp),
        "specificity": specificity,
        "balanced_accuracy": mean([sensitivity, specificity]) if sensitivity is not None and specificity is not None else None,
        "matthews_correlation": (tp * tn - fp * fn) / mcc_denominator if mcc_denominator else None,
        "predicted_missing_bp_total": sum(predicted),
        "observed_gain_bp_total": sum(observed),
        "missing_bp_mean_absolute_error": mean(absolute_errors) if absolute_errors else None,
        "missing_bp_median_absolute_error": median(absolute_errors) if absolute_errors else None,
        "missing_bp_pearson": _pearson(predicted, observed),
        "missing_bp_spearman": _pearson(_rank(predicted), _rank(observed)) if predicted else None,
        "min_new_bp_sensitivity": sensitivity_rows,
    }


def run_evaluation(
    copy_number: Path,
    old_arrays: Path,
    new_arrays: Path,
    outdir: Path,
    config: EvaluationConfig,
    sensitivity_arrays: Path | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    read_bp = read_copy_number(copy_number)
    rows = evaluate_families(
        read_bp,
        assembly_bp_by_family(old_arrays),
        assembly_bp_by_family(new_arrays),
        assembly_bp_by_family(sensitivity_arrays) if sensitivity_arrays else None,
        config,
    )
    summary = summarize_rows(rows, config)
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir / "family_metrics.tsv", rows, FAMILY_FIELDS)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return rows, summary
