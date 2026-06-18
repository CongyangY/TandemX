from __future__ import annotations

import random
import sys

import pytest
import tandemx

from tandemx.discover.rust_backend import (
    RustBackendUnavailable,
    rust_backend_available,
    scan_read_for_periods,
)
from tandemx.discover.spacing import (
    build_spacing_histogram,
    extract_repeated_kmer_positions,
    refine_candidate_period,
    select_candidate_periods,
)


pytestmark = pytest.mark.skipif(
    not rust_backend_available(),
    reason="TandemX Rust extension is not installed",
)


def rust_scan(sequence: str, min_period: int = 20, max_period: int = 100):
    return scan_read_for_periods(
        sequence,
        k=11,
        min_period=min_period,
        max_period=max_period,
        top_periods=3,
        min_seed_occurrences=2,
        min_spacing_support=2,
        max_pairs_per_kmer=100,
    )


def python_scan(sequence: str, min_period: int = 20, max_period: int = 100):
    positions, _ = extract_repeated_kmer_positions(sequence, 11, 2, 100)
    histogram = build_spacing_histogram(positions, min_period, max_period, 100)
    candidates = select_candidate_periods(histogram, min_period, max_period, 3, 2)
    best_period, score = refine_candidate_period(
        sequence,
        positions,
        candidates,
        min_period,
        max_period,
    )
    return candidates, best_period, score


def test_rust_detects_simple_and_non_default_repeats() -> None:
    for monomer in ("ACGTTCAGGACTAACCGTGA", "ACGTTCAGGACTAACCGTGATCGATCGATCG"):
        result = rust_scan(monomer * 20)
        assert result.best_period == len(monomer)
        assert result.periodicity_score >= 0.9
        assert result.status == "accepted"


def test_rust_matches_python_period_and_score() -> None:
    sequence = "ACGTTCAGGACTAACCGTGATCGATCGATCG" * 20
    python_candidates, python_period, python_score = python_scan(sequence)
    rust_result = rust_scan(sequence)

    assert rust_result.candidate_periods == tuple(python_candidates)
    assert rust_result.best_period == python_period
    assert rust_result.periodicity_score == pytest.approx(python_score)


def test_rust_handles_random_low_complexity_and_ambiguous_reads() -> None:
    random_sequence = "".join(random.Random(7).choices("ACGT", k=2000))
    assert rust_scan(random_sequence).status != "accepted"
    assert rust_scan("A" * 2000).candidate_periods == ()

    monomer = "ACGTTCAGGACTAACCGTGA"
    with_n = monomer * 10 + "N" + monomer * 10
    assert rust_scan(with_n).best_period == len(monomer)


def test_rust_backend_unavailable_message(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "tandemx._rust_core", None)
    monkeypatch.delattr(tandemx, "_rust_core", raising=False)
    with pytest.raises(RustBackendUnavailable, match="maturin develop"):
        rust_scan("ACGT" * 100)
