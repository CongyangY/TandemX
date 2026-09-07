from __future__ import annotations

import pytest

from benchmarks.scripts.archive_cascade_gap_free_development import (
    runtime_geomean_ratio,
)


def test_runtime_geomean_ratio_uses_paired_scenarios() -> None:
    rows = [
        {"scenario": "a", "tool": "tandemx", "median_runtime_seconds": "2"},
        {"scenario": "a", "tool": "tidehunter", "median_runtime_seconds": "1"},
        {"scenario": "b", "tool": "tandemx", "median_runtime_seconds": "8"},
        {"scenario": "b", "tool": "tidehunter", "median_runtime_seconds": "2"},
    ]
    assert runtime_geomean_ratio(rows) == pytest.approx((2 * 4) ** 0.5)


def test_runtime_geomean_ratio_rejects_nonpositive_measurement() -> None:
    rows = [
        {"scenario": "a", "tool": "tandemx", "median_runtime_seconds": "0"},
        {"scenario": "a", "tool": "tidehunter", "median_runtime_seconds": "1"},
    ]
    with pytest.raises(ValueError, match="positive"):
        runtime_geomean_ratio(rows)
