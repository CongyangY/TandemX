"""Mechanical metric tests; these toy cases are not holdout evaluations."""

import pytest

from benchmarks.challenge.phase_gate_metrics import score_predictions


LENGTHS = {"r1": 200}
TRUTH = [("r1", 20, 120, "F1")]
TRUTH_TWO = TRUTH + [("r1", 0, 1, "F2")]


def test_perfect_prediction_has_full_precision_recall() -> None:
    result = score_predictions({"r1": {"F1": [(20, 120)]}}, TRUTH, LENGTHS, ["F1"])
    row = result["per_family"]["F1"]
    assert row["read_truth_bp"] == 100
    assert row["predicted_unique_bp"] == 100
    assert row["true_positive_bp"] == 100
    assert row["false_positive_bp"] == 0
    assert row["recall"] == row["precision"] == 1.0
    assert row["absolute_relative_error"] == 0.0
    assert result["overall"]["mass_conserved"]


def test_wrong_family_is_false_positive_not_background() -> None:
    result = score_predictions({"r1": {"F2": [(20, 120)]}}, TRUTH_TWO, LENGTHS, ["F1", "F2"])
    assert result["per_family"]["F1"]["true_positive_bp"] == 0
    assert result["per_family"]["F2"]["false_positive_bp"] == 100
    assert result["per_family"]["F2"]["background_false_bp"] == 0


def test_prediction_on_background_is_background_false() -> None:
    result = score_predictions({"r1": {"F1": [(150, 180)]}}, TRUTH, LENGTHS, ["F1"])
    row = result["per_family"]["F1"]
    assert row["false_positive_bp"] == row["background_false_bp"] == 30
    assert row["absolute_relative_error"] == 0.7


def test_cross_family_overlap_is_ambiguous_and_not_forced_to_either() -> None:
    predictions = {"r1": {"F1": [(20, 120)], "F2": [(70, 170)]}}
    result = score_predictions(predictions, TRUTH_TWO, LENGTHS, ["F1", "F2"])
    assert result["overall"]["ambiguous_bp"] == 50
    assert result["per_family"]["F1"]["ambiguous_eligible_bp"] == 50
    assert result["per_family"]["F2"]["ambiguous_eligible_bp"] == 50
    assert result["per_family"]["F1"]["true_positive_bp"] == 50
    assert result["per_family"]["F2"]["true_positive_bp"] == 0


def test_empty_predictions_report_unassigned_bases() -> None:
    result = score_predictions({}, TRUTH, LENGTHS, ["F1"])
    assert result["overall"]["predicted_unique_bp"] == 0
    assert result["overall"]["ambiguous_bp"] == 0
    assert result["overall"]["unassigned_bp"] == 200
    assert result["per_family"]["F1"]["recall"] == 0.0


def test_zero_truth_has_undefined_relative_error_and_reports_false_bp() -> None:
    result = score_predictions({"empty": {"F1": [(5, 15)]}}, [], {"empty": 80}, ["F1"])
    assert result["per_family"]["F1"]["read_truth_bp"] == 0
    assert result["per_family"]["F1"]["absolute_relative_error"] is None
    assert result["per_family"]["F1"]["background_false_bp"] == 10
    assert result["overall"] == {
        "ambiguous_bp": 0,
        "unassigned_bp": 70,
        "predicted_unique_bp": 10,
        "read_bases": 80,
        "mass_conserved": True,
    }


def test_negative_read_length_fails_fast() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        score_predictions({}, [], {"bad": -1}, ["F1"])


@pytest.mark.parametrize(
    "predictions, message",
    [
        ({"missing": {"F1": [(0, 1)]}}, "Unknown prediction read ID"),
        ({"r1": {"F9": [(0, 1)]}}, "Unknown prediction family"),
        ({"r1": {"F1": [(0, 201)]}}, "exceeds read length"),
    ],
)
def test_invalid_prediction_coordinates_and_ids_fail_fast(predictions, message) -> None:
    with pytest.raises(ValueError, match=message):
        score_predictions(predictions, TRUTH, LENGTHS, ["F1"])
