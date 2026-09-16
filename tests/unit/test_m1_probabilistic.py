"""Focused checks for the development-only read-likelihood candidate."""
from __future__ import annotations

import pytest

from benchmarks.m1_shared_signature.probabilistic import (
    Settings, distance_matrix, infer,
)
from benchmarks.m1_shared_signature.run import SCENARIOS, make_case


def test_exact_reads_and_random_negative_rejection() -> None:
    catalogue = {"a": "AAAACCCCGGGG", "b": "CCCCGGGGAAAA",
                 "decoy": "TTTTGGGGCCCC"}
    # a/b are cyclic-equivalent; they must remain an unresolved group.
    result = infer(["AAAACCCCGGGG"] * 5 + ["ACGTACGTACGT"], catalogue,
                   Settings(max_mismatch=1))
    assert result["family"]["a"] is None
    assert result["family"]["b"] is None
    assert result["groups"]["a+b"] == 5
    assert result["family"]["decoy"] == 0
    assert result["rejected_gate"] == 1


def test_distinct_families_are_counted_without_truth_labels() -> None:
    catalogue = {"a": "AAAACCCCGGGG", "b": "TTTTGGGGCCCC"}
    result = infer(["AAAACCCCGGGG"] * 4 + ["TTTTGGGGCCCC"] * 3,
                   catalogue, Settings(max_mismatch=1))
    assert result["family"] == {"a": 4.0, "b": 3.0}
    assert result["assigned"] == 7
    assert result["rejected_unknown"] == 0


def test_empty_invalid_and_reproducible() -> None:
    catalogue = {"a": "AAAA", "b": "CCCC"}
    result = infer([], catalogue)
    assert result["family"] == {"a": 0.0, "b": 0.0}
    with pytest.raises(ValueError):
        infer(["NNNN"], catalogue)
    with pytest.raises(ValueError):
        infer(["AAAA"], catalogue, Settings(posterior_min=0.3))
    with pytest.raises(ValueError):
        distance_matrix(["AAAA"], {"a": "AAA"})
    case_catalogue, reads, _ = make_case(11, *SCENARIOS[0][1:])
    assert infer(reads, case_catalogue) == infer(reads, case_catalogue)
