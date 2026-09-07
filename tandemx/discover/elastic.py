"""Seed-selected, indel-aware discovery of multiple read-local repeat arrays."""
from __future__ import annotations

from math import ceil
from collections import Counter

from tandemx.discover.alignment import AlignmentHit, banded_self_align_many
from tandemx.discover.consensus import aligned_unit_consensus
from tandemx.discover.spacing import (
    best_local_periodicity_score,
    build_spacing_histogram,
    expand_candidate_periods,
    extract_repeated_kmer_positions,
)


CASCADE_GAP_FREE_MIN_READ_FRACTION = 0.30
CASCADE_GAP_FREE_MIN_SHIFTED_IDENTITY = 0.95
CASCADE_GAP_FREE_MIN_VALID_PAIR_FRACTION = 0.95
CASCADE_GAP_FREE_MAX_UNIT_RESIDUAL_FRACTION = 0.02


def _unit_span_residual_fraction(span: int, period: int) -> float:
    remainder = span % period
    return min(remainder, period - remainder) / period


def select_alignment_periods(histogram: dict[int, int], top_periods: int,
                             minimum_support: int) -> list[int]:
    """Spend the period budget on distinct alignment bands, including harmonics."""
    periods: list[int] = []
    for period, support in sorted(histogram.items(), key=lambda item: (-item[1], item[0])):
        if support < minimum_support:
            continue
        if any(abs(period - old) <= max(3, ceil(old * 0.08)) for old in periods):
            continue
        periods.append(period)
        if len(periods) == top_periods:
            break
    return periods


def suppress_overlapping_hits(hits: list[AlignmentHit]) -> list[AlignmentHit]:
    """Retain the strongest alignment per locus; do not count harmonics twice.

    At least half the shorter interval must overlap for suppression. This is
    an explicit reporting rule, not evidence that higher-order repeats are absent.
    """
    retained: list[AlignmentHit] = []
    for hit in sorted(hits, key=lambda h: (-h.alignment_score, -(h.end - h.start), h.period, h.start)):
        if any(max(0, min(hit.end, other.end) - max(hit.start, other.start)) * 2 >=
               min(hit.end - hit.start, other.end - other.start) for other in retained):
            continue
        retained.append(hit)
    return sorted(retained, key=lambda h: (h.start, h.end, h.period))


def composition_adjusted_identity(sequence: str, hit: AlignmentHit) -> float:
    """Chance-correct identity by local mononucleotide composition.

    This is a heuristic selectivity filter, not a calibrated P value. The null
    equality rate q=sum(f_b^2) accounts for compositional bias; it does not model
    dinucleotide dependencies, homology between units or selection of a local hit.
    """
    counts = Counter(sequence[hit.start:hit.end])
    total = sum(counts[b] for b in "ACGT")
    if not total:
        return 0.0
    chance = sum((counts[b] / total) ** 2 for b in "ACGT")
    return (hit.identity - chance) / (1 - chance) if chance < 1 else 0.0


def _validate_parameters(min_period: int, max_period: int, min_span: int, k: int,
                         top_periods: int, min_seed_occurrences: int,
                         min_spacing_support: int, max_pairs_per_kmer: int) -> None:
    if min_period <= 0 or max_period < min_period or min_span <= 0 or k <= 0:
        raise ValueError("Invalid elastic discovery period, span or k-mer size")
    if top_periods <= 0 or min_seed_occurrences < 2 or min_spacing_support <= 0 or max_pairs_per_kmer <= 0:
        raise ValueError("Invalid elastic discovery seed parameters")


def _seed_histogram(sequence: str, *, k: int, min_period: int, max_period: int,
                    min_seed_occurrences: int, max_pairs_per_kmer: int,
                    backend: str) -> tuple[dict[int, int], int]:
    if backend == "rust":
        from tandemx.discover.rust_backend import seed_spacing_histogram
        return seed_spacing_histogram(sequence, k, max(20, min_period), max_period,
                                      min_seed_occurrences, max_pairs_per_kmer)
    positions, overflow = extract_repeated_kmer_positions(
        sequence, k, min_seed_occurrences, max_pairs_per_kmer
    )
    return build_spacing_histogram(
        positions, max(20, min_period), max_period, max_pairs_per_kmer
    ), overflow


def _elastic_calls(sequence: str, histogram: dict[int, int], *, min_period: int,
                   max_period: int, min_span: int, top_periods: int,
                   min_spacing_support: int, backend: str
                   ) -> list[tuple[AlignmentHit, str, int]]:
    short_periods = list(range(min_period, min(19, max_period) + 1))
    periods = short_periods + select_alignment_periods(histogram, top_periods, min_spacing_support)
    hits = banded_self_align_many(sequence, periods, min_span, backend=backend)
    result = []
    eligible = [hit for hit in hits if composition_adjusted_identity(sequence, hit) >= 0.7]
    for hit in suppress_overlapping_hits(eligible):
        consensus, units = aligned_unit_consensus(sequence, hit, backend=backend)
        if not min_period <= len(consensus) <= max_period:
            continue
        result.append((hit, consensus, units))
    return result


def _dominant_clean_call(sequence: str, histogram: dict[int, int], *, min_period: int,
                         max_period: int, min_span: int, top_periods: int,
                         min_spacing_support: int, backend: str
                         ) -> tuple[AlignmentHit, str, int] | None:
    """Return a gap-free dominant array without allocating an alignment trace.

    A development audit chose a 30% read fraction and 95% shifted identity while
    retaining the composition-adjusted identity filter. The span guard keeps the
    observed two-array stress case on the elastic path. Reads outside this domain
    retain indel-aware alignment and multiple-array behavior.
    """
    peaks = list(range(min_period, min(19, max_period) + 1))
    peaks += select_alignment_periods(histogram, top_periods, min_spacing_support)
    periods = expand_candidate_periods(peaks, min_period, max_period, refinement_radius=3)
    minimum_dominant_span = max(
        min_span, ceil(CASCADE_GAP_FREE_MIN_READ_FRACTION * len(sequence))
    )
    best: tuple[tuple[int, float, int], AlignmentHit] | None = None
    for period in periods:
        identity, start, end = best_local_periodicity_score(
            sequence,
            period,
            min_span,
            acceptance_score=CASCADE_GAP_FREE_MIN_SHIFTED_IDENTITY,
        )
        span = end - start
        if (
            identity < CASCADE_GAP_FREE_MIN_SHIFTED_IDENTITY
            or span < minimum_dominant_span
            or _unit_span_residual_fraction(span, period)
            > CASCADE_GAP_FREE_MAX_UNIT_RESIDUAL_FRACTION
        ):
            continue
        pairs = tuple(
            (left, left + period)
            for left in range(start, end - period)
            if sequence[left] in "ACGT" and sequence[left + period] in "ACGT"
        )
        possible_columns = max(0, end - start - period)
        if (
            not pairs
            or len(pairs) / possible_columns
            < CASCADE_GAP_FREE_MIN_VALID_PAIR_FRACTION
        ):
            continue
        matches = sum(sequence[left] == sequence[right] for left, right in pairs)
        columns = len(pairs)
        hit = AlignmentHit(
            start, end, period, matches, columns, 0,
            2 * matches - 3 * (columns - matches), pairs,
        )
        if composition_adjusted_identity(sequence, hit) < 0.7:
            continue
        key = (span, identity, -period)
        if best is None or key > best[0]:
            best = key, hit
    if best is None:
        return None
    hit = best[1]
    consensus, units = aligned_unit_consensus(sequence, hit, backend=backend)
    if not min_period <= len(consensus) <= max_period:
        return None
    return hit, consensus, units


def _verified_gap_free_call(sequence: str, period: int, start: int, end: int,
                            *, min_period: int, max_period: int, min_span: int,
                            backend: str) -> tuple[AlignmentHit, str, int] | None:
    """Verify one native clean-path proposal without another period scan."""
    if not min_period <= period <= max_period:
        return None
    start = max(0, start)
    end = min(len(sequence), end)
    if end - start < max(
        min_span, ceil(CASCADE_GAP_FREE_MIN_READ_FRACTION * len(sequence))
    ):
        return None
    if (
        _unit_span_residual_fraction(end - start, period)
        > CASCADE_GAP_FREE_MAX_UNIT_RESIDUAL_FRACTION
    ):
        return None
    pairs = tuple(
        (left, left + period)
        for left in range(start, end - period)
        if sequence[left] in "ACGT" and sequence[left + period] in "ACGT"
    )
    possible_columns = max(0, end - start - period)
    if (
        not pairs
        or len(pairs) / possible_columns
        < CASCADE_GAP_FREE_MIN_VALID_PAIR_FRACTION
    ):
        return None
    matches = sum(sequence[left] == sequence[right] for left, right in pairs)
    columns = len(pairs)
    if matches / columns < CASCADE_GAP_FREE_MIN_SHIFTED_IDENTITY:
        return None
    hit = AlignmentHit(
        start, end, period, matches, columns, 0,
        2 * matches - 3 * (columns - matches), pairs,
    )
    if composition_adjusted_identity(sequence, hit) < 0.7:
        return None
    consensus, units = aligned_unit_consensus(sequence, hit, backend=backend)
    if not min_period <= len(consensus) <= max_period:
        return None
    return hit, consensus, units


def discover_elastic_arrays(sequence: str, *, min_period: int, max_period: int,
                            min_span: int, k: int = 11, top_periods: int = 5,
                            min_seed_occurrences: int = 2, min_spacing_support: int = 2,
                            max_pairs_per_kmer: int = 100, backend: str = "python"
                            ) -> tuple[list[tuple[AlignmentHit, str, int]], int]:
    """Return alignment evidence, consensus and unit count for each retained locus."""
    _validate_parameters(min_period, max_period, min_span, k, top_periods,
                         min_seed_occurrences, min_spacing_support, max_pairs_per_kmer)
    max_period = min(max_period, len(sequence) // 2)
    if min_period > max_period or len(sequence) < min_span:
        return [], 0
    histogram, overflow = _seed_histogram(
        sequence, k=k, min_period=min_period, max_period=max_period,
        min_seed_occurrences=min_seed_occurrences, max_pairs_per_kmer=max_pairs_per_kmer,
        backend=backend,
    )
    return _elastic_calls(
        sequence, histogram, min_period=min_period, max_period=max_period,
        min_span=min_span, top_periods=top_periods,
        min_spacing_support=min_spacing_support, backend=backend,
    ), overflow


def discover_cascade_arrays(sequence: str, *, min_period: int, max_period: int,
                            min_span: int, k: int = 11, top_periods: int = 5,
                            min_seed_occurrences: int = 2, min_spacing_support: int = 2,
                            max_pairs_per_kmer: int = 100, backend: str = "python"
                            ) -> tuple[list[tuple[AlignmentHit, str, int, str]], int]:
    """Screen, use a narrow gap-free fast path, then fall back to elastic alignment."""
    _validate_parameters(min_period, max_period, min_span, k, top_periods,
                         min_seed_occurrences, min_spacing_support, max_pairs_per_kmer)
    max_period = min(max_period, len(sequence) // 2)
    if min_period > max_period or len(sequence) < min_span:
        return [], 0
    if backend == "rust":
        from tandemx.discover.rust_backend import scan_read_for_periods

        screen = scan_read_for_periods(
            sequence, k=k, min_period=min_period, max_period=max_period,
            top_periods=top_periods, min_seed_occurrences=min_seed_occurrences,
            min_spacing_support=min_spacing_support,
            max_pairs_per_kmer=max_pairs_per_kmer, min_repeat_span=min_span,
        )
        histogram = dict(screen.spacing_support)
        overflow = screen.overflow_count
        clean = _verified_gap_free_call(
            sequence, screen.best_period, screen.repeat_start, screen.repeat_end,
            min_period=min_period, max_period=max_period, min_span=min_span,
            backend=backend,
        )
    else:
        histogram, overflow = _seed_histogram(
            sequence, k=k, min_period=min_period, max_period=max_period,
            min_seed_occurrences=min_seed_occurrences,
            max_pairs_per_kmer=max_pairs_per_kmer, backend=backend,
        )
        clean = _dominant_clean_call(
            sequence, histogram, min_period=min_period, max_period=max_period,
            min_span=min_span, top_periods=top_periods,
            min_spacing_support=min_spacing_support, backend=backend,
        )
    if clean is not None:
        return [(*clean, "cascade_gap_free")], overflow
    calls = _elastic_calls(
        sequence, histogram, min_period=min_period, max_period=max_period,
        min_span=min_span, top_periods=top_periods,
        min_spacing_support=min_spacing_support, backend=backend,
    )
    return [(*call, "cascade_elastic") for call in calls], overflow
