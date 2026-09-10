"""Independent unified endpoint tests; toy only, never a holdout benchmark."""

import pytest

from benchmarks.challenge.schema import ArrayRecord
from benchmarks.challenge.unified_metrics import interval_metrics


LENGTHS = {"positive": 200, "negative": 100}
TRUTH = [ArrayRecord("positive", 20, 120, 120, family_id="F1")]


def test_union_overlap_and_duplicate_predictions() -> None:
    predictions = [
        ArrayRecord("positive", 20, 80, 120, family_id="native1"),
        ArrayRecord("positive", 50, 120, 120, family_id="native1"),
    ]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {"native1": "F1"})
    assert result["global"]["predicted_union_bp"] == 100
    assert result["global"]["true_positive_bp"] == 100
    assert result["per_family"]["F1"]["predicted_bp"] == 100
    assert result["per_family"]["F1"]["false_positive_bp"] == 0


def test_cross_family_overlap_is_ambiguous_but_global_precision_keeps_union() -> None:
    predictions = [
        ArrayRecord("positive", 20, 120, 120, family_id="native1"),
        ArrayRecord("positive", 70, 170, 120, family_id="native2"),
    ]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {"native1": "F1", "native2": "F2"})
    assert result["global"]["predicted_union_bp"] == 150
    assert result["global"]["true_positive_bp"] == 100
    assert result["global"]["base_union_precision"] == pytest.approx(2 / 3)
    assert result["per_family"]["F1"]["ambiguous_bp"] == 50
    assert result["per_family"]["F1"]["true_positive_bp"] == 50
    assert result["per_family"]["F2"]["true_positive_bp"] == 0


def test_negative_read_prediction_is_reported() -> None:
    predictions = [ArrayRecord("negative", 10, 40, 120, family_id="native1")]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {"native1": "F1"})
    assert result["global"]["negative_read_predicted_bp"] == 30
    assert result["global"]["base_union_precision"] == 0.0
    assert result["per_family"]["F1"]["false_positive_bp"] == 30


def test_unknown_native_mapping_stays_in_global_union_and_is_unassigned() -> None:
    predictions = [ArrayRecord("positive", 20, 120, 120, family_id="native_unknown")]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {})
    assert result["global"]["predicted_union_bp"] == 100
    assert result["global"]["base_union_precision"] == 1.0
    assert result["global"]["unassigned_native_bp"] == 100
    assert result["per_family"]["F1"]["true_positive_bp"] == 0


@pytest.mark.parametrize(
    "predictions, truths, lengths, message",
    [
        ([ArrayRecord("missing", 0, 1, 120, family_id="n")], TRUTH, LENGTHS, "Unknown read ID"),
        ([ArrayRecord("positive", 0, 201, 120, family_id="n")], TRUTH, LENGTHS, "Interval"),
        (TRUTH, TRUTH, {"positive": 50}, "Interval"),
    ],
)
def test_invalid_read_or_interval_fails(predictions, truths, lengths, message) -> None:
    with pytest.raises(ValueError, match=message):
        interval_metrics(predictions, truths, lengths, {"n": "F1"})


def test_unmatched_and_ambiguous_mapping_values_are_not_family_correct() -> None:
    predictions = [ArrayRecord("positive", 20, 120, 120, family_id="native1")]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {"native1": None})
    assert result["global"]["unassigned_native_bp"] == 100
    assert result["per_family"]["F1"]["true_positive_bp"] == 0


def test_overlapping_unmatched_native_families_count_once_per_read() -> None:
    predictions = [
        ArrayRecord("positive", 20, 100, 120, family_id="native1"),
        ArrayRecord("positive", 70, 150, 120, family_id="native2"),
    ]
    result = interval_metrics(predictions, TRUTH, LENGTHS, {})
    assert result["global"]["unassigned_native_bp"] == 130
