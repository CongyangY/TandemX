from __future__ import annotations

import pytest

from benchmarks.scripts.analyze_cascade_fast_path import (
    interval_iou,
    parse_grid,
    shifted_evidence,
    threshold_rows,
)


def _feature(**updates: object) -> dict[str, object]:
    row: dict[str, object] = {
        "best_period": 4,
        "repeat_span": 40,
        "read_length": 100,
        "shifted_identity": 1.0,
        "valid_pair_fraction": 1.0,
        "unit_span_residual_fraction": 0.0,
        "composition_adjusted_identity": 1.0,
        "truth_array_count": 1,
        "maximum_truth_iou": 1.0,
        "minimum_truth_period_error_bp": 0,
    }
    row.update(updates)
    return row


def test_shifted_evidence_and_interval_iou() -> None:
    identity, valid_fraction, adjusted = shifted_evidence("ACGT" * 5, 4, 0, 20)
    assert identity == 1.0
    assert valid_fraction == 1.0
    assert adjusted == 1.0
    assert interval_iou((10, 30), (20, 40)) == pytest.approx(1 / 3)


def test_threshold_grid_marks_truth_failures_without_using_them_in_rule() -> None:
    rows = [
        _feature(),
        _feature(truth_array_count=0, composition_adjusted_identity=0.1),
        _feature(truth_array_count=2),
    ]
    result = threshold_rows(rows, [0.3], [0.995], minimum_span=1)
    assert result[0]["accepted_reads"] == 2
    assert result[0]["correct_single_array_acceptances"] == 1
    assert result[0]["negative_acceptances"] == 0
    assert result[0]["multi_array_acceptances"] == 1
    assert result[0]["development_safe"] is False


def test_parse_grid_rejects_values_outside_unit_interval() -> None:
    assert parse_grid("0,0.5,1") == [0.0, 0.5, 1.0]
    with pytest.raises(ValueError, match="in \\[0,1\\]"):
        parse_grid("1.1")
