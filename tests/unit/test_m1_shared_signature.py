"""Regression tests for the development-only M1 tournament."""
from __future__ import annotations

import numpy as np
import pytest

from benchmarks.m1_shared_signature.model import (
    align_and_gate, circular_kmers, constrained_fit, shared_signature_fit,
)
from benchmarks.m1_shared_signature.run import SCENARIOS, evaluate, make_case


def test_rotation_invariant_signature_and_mapping() -> None:
    catalogue = {"a": "AAAACCCC", "b": "GGGGTTTT"}
    assert circular_kmers("AAAACCCC", 3) == circular_kmers("CCCCAAAA", 3)
    accepted, mapped, rejected, tied = align_and_gate(
        ["CCCCAAAA", "GGGGTTTT", "ACGTACGT"], catalogue, 0)
    assert accepted == ["CCCCAAAA", "GGGGTTTT"]
    assert mapped == {"a": 1.0, "b": 1.0}
    assert (rejected, tied) == (1, 0)


def test_identical_columns_refuse_individual_values() -> None:
    reads = ["AAAACCCC"] * 4
    family, groups, members = shared_signature_fit(
        reads, {"a": "AAAACCCC", "b": "CCCCAAAA"}, 3)
    assert family == {"a": None, "b": None}
    assert groups["a+b"] == pytest.approx(4.0)
    assert members == [["a", "b"]]


def test_near_degenerate_low_complexity_design_refuses() -> None:
    a = "A" * 200
    b = "A" * 100 + "C" + "A" * 99
    family, groups, _ = shared_signature_fit([a], {"a": a, "b": b}, 5)
    assert family == {"a": None, "b": None}
    assert groups == {"a+b": 1.0}


def test_constrained_fit_nonnegative_and_conserved() -> None:
    columns = np.array([[3.0, 0.0], [0.0, 3.0]])
    estimate = constrained_fit(columns, np.array([6.0, 12.0]), 6)
    assert estimate == pytest.approx([2.0, 4.0])
    assert np.all(constrained_fit(columns, np.zeros(2), 0) == 0)


def test_empty_invalid_and_reproducible() -> None:
    with pytest.raises(ValueError):
        align_and_gate([], {}, 0)
    with pytest.raises(ValueError):
        circular_kmers("AAAA", 5)
    with pytest.raises(ValueError):
        align_and_gate(["NNNN"], {"a": "AAAA"}, 1)
    empty_fit, empty_groups, _ = shared_signature_fit([], {"a": "AAAA"}, 2)
    assert empty_fit == {"a": 0.0}
    assert empty_groups == {}
    assert make_case(7, 3, 0.03, 0.0) == make_case(7, 3, 0.03, 0.0)


def test_case_records_negatives_and_individual_refusal() -> None:
    result = evaluate(11, SCENARIOS[-1])
    assert result["n_reads"] == 600
    assert result["accepted_by_source"]["unknown"] <= 100
    assert result["scores"]["shared_signature"]["individual_refused"] == 2
    assert result["shared_signature"]["f1"] is None
