from __future__ import annotations

from tandemx.discover.mvp import is_low_complexity
from tandemx.discover.spacing import (
    build_spacing_histogram,
    extract_repeated_kmer_positions,
    is_low_complexity_kmer,
    modulo_periodicity_score,
    select_candidate_periods,
)


def test_extract_repeated_kmer_positions() -> None:
    sequence = "AACCGGTT" * 6

    positions, overflow_count = extract_repeated_kmer_positions(
        sequence,
        kmer_size=7,
        min_seed_occurrences=2,
        max_pairs_per_kmer=20,
    )

    assert positions
    assert all(len(observed) >= 2 for observed in positions.values())
    assert overflow_count == 0


def test_rolling_kmers_reset_at_ambiguous_bases() -> None:
    seed = "ACGTTCAGGAC"
    positions, overflow_count = extract_repeated_kmer_positions(
        f"{seed}N{seed}",
        kmer_size=len(seed),
        min_seed_occurrences=2,
    )

    assert list(positions.values()) == [[0, 12]]
    assert overflow_count == 0


def test_spacing_histogram_and_candidate_selection() -> None:
    positions = {
        "seed_a": [3, 24, 45, 66, 87],
        "seed_b": [10, 31, 52, 73],
    }

    histogram = build_spacing_histogram(positions, min_period=15, max_period=30, bin_size=1)
    candidates = select_candidate_periods(
        histogram,
        min_period=15,
        max_period=30,
        top_periods=3,
        min_spacing_support=3,
    )

    assert histogram[21] >= 7
    assert candidates[0] == 21


def test_modulo_periodicity_score_prefers_supported_period() -> None:
    positions = {
        "seed_a": [2, 23, 44, 65],
        "seed_b": [8, 29, 50, 71],
    }

    assert modulo_periodicity_score(positions, 21, tolerance=0) == 1.0
    assert modulo_periodicity_score(positions, 20, tolerance=0) < 1.0


def test_modulo_periodicity_score_wraps_across_period_boundary() -> None:
    positions = {"seed": [19, 20, 39, 40, 59, 60]}

    assert modulo_periodicity_score(positions, 20, tolerance=1) == 1.0


def test_low_complexity_filters_reads_and_seeds() -> None:
    assert is_low_complexity("A" * 100)
    assert is_low_complexity("AT" * 50)
    assert not is_low_complexity("ACGTACGA" * 12)
    assert is_low_complexity_kmer("AAAAAAAAAAA")
    assert is_low_complexity_kmer("ATATATATATA")
    assert not is_low_complexity_kmer("ACGTTCAGGAC")


def test_max_pairs_per_kmer_bounds_position_storage_and_histogram() -> None:
    sequence = "ACGTTCAGGAC" * 100
    positions, overflow_count = extract_repeated_kmer_positions(
        sequence,
        kmer_size=11,
        min_seed_occurrences=2,
        max_pairs_per_kmer=3,
    )

    assert overflow_count > 0
    assert all(len(observed) <= 4 for observed in positions.values())
    histogram = build_spacing_histogram(
        {"seed": list(range(0, 1000, 10))},
        min_period=10,
        max_period=100,
        max_pairs_per_kmer=3,
        bin_size=1,
    )
    assert sum(histogram.values()) <= 3
