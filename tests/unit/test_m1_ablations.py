"""Focused behavior tests for fixed M1 A/C/D development candidates."""
from __future__ import annotations

import pytest

from benchmarks.m1_shared_signature import ablations


def test_error_aware_design_conserves_expected_windows() -> None:
    groups, words, expected, _ = ablations.design({"a": "AAAACCCC", "b": "GGGGTTTT"})
    assert len(groups) == 2
    assert len(words) == 4 ** ablations.K
    assert list(expected.sum(axis=0)) == pytest.approx([8.0, 8.0])


def test_distinct_exact_units_and_zero_decoy() -> None:
    catalogue = {"a": "AAAACCCC", "b": "GGGGTTTT", "zero": "ACGTACGT"}
    reads = [catalogue["a"]] * 4 + [catalogue["b"]] * 3
    for implementation in (ablations.discriminative, ablations.poisson_counts,
                           ablations.em_mixture):
        result = implementation(reads, catalogue)
        assert result["family"]["a"] > result["family"]["zero"]
        assert result["family"]["b"] > result["family"]["zero"]
    mixture = ablations.em_mixture(reads, catalogue)
    assert mixture["assigned"] + mixture["rejected_ambiguous"] == len(reads)
    assert mixture["family"]["zero"] == 0


def test_identical_signatures_refuse_individual_counts() -> None:
    catalogue = {"a": "AAAACCCC", "b": "CCCCAAAA"}
    for implementation in (ablations.discriminative, ablations.poisson_counts,
                           ablations.em_mixture):
        result = implementation(["AAAACCCC"] * 5, catalogue)
        assert result["family"] == {"a": None, "b": None}
        assert result["groups"]["a+b"] > 0


def test_low_complexity_near_degenerate_design_refuses() -> None:
    a = "A" * 200
    b = "A" * 100 + "C" + "A" * 99
    for implementation in (ablations.discriminative, ablations.poisson_counts,
                           ablations.em_mixture):
        result = implementation([a], {"a": a, "b": b})
        assert result["family"] == {"a": None, "b": None}
        assert result["groups"] == {"a+b": 1.0}
        assert result["design_refusal"] == "singular_ratio_below_0.05"


def test_empty_invalid_and_reproducible() -> None:
    catalogue = {"a": "AAAACCCC"}
    for implementation in (ablations.discriminative, ablations.poisson_counts,
                           ablations.em_mixture):
        empty = implementation([], catalogue)
        assert empty["family"]["a"] == 0
        assert implementation(["AAAACCCC"], catalogue) == implementation(["AAAACCCC"], catalogue)
    with pytest.raises(ValueError):
        ablations.design({"a": "NNNNNNNN"})
    with pytest.raises(ValueError):
        ablations.discriminative(["NNNNNNNN"], catalogue)
