import pytest

from benchmarks.scripts.plot_ey15_donor_validation import (
    annotation_class_recall,
    metric_rows,
)


def test_metric_rows_preserve_wilson_intervals_and_normalization() -> None:
    summary = {
        "sensitivity": 0.8,
        "sensitivity_wilson95": [0.5, 0.95],
        "precision": 0.75,
        "precision_wilson95": [0.45, 0.92],
        "false_positive_rate": 0.1,
        "false_positive_rate_wilson95": [0.02, 0.35],
    }
    rows = metric_rows(summary, "depth107")
    assert [row["metric"] for row in rows] == [
        "Sensitivity",
        "Precision",
        "False-positive rate",
    ]
    assert rows[0] == {
        "normalization": "depth107",
        "metric": "Sensitivity",
        "value": 0.8,
        "lower": 0.5,
        "upper": 0.95,
    }


def test_annotation_class_recall_uses_the_documented_column() -> None:
    rows = [
        {"annotation_class": "centromere", "annotation_bp_recall": "0.03175"},
        {"annotation_class": "5S_rDNA", "annotation_bp_recall": "0.48134"},
    ]
    assert annotation_class_recall(rows, "centromere") == pytest.approx(0.03175)


def test_annotation_class_recall_requires_one_matching_row() -> None:
    with pytest.raises(ValueError, match="expected one annotation-class row"):
        annotation_class_recall([], "centromere")
