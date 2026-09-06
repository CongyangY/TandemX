"""Seed-selected, indel-aware discovery of multiple read-local repeat arrays."""
from __future__ import annotations

from math import ceil
from collections import Counter

from tandemx.discover.alignment import AlignmentHit, banded_self_align
from tandemx.discover.consensus import aligned_unit_consensus
from tandemx.discover.spacing import extract_repeated_kmer_positions, build_spacing_histogram


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


def discover_elastic_arrays(sequence: str, *, min_period: int, max_period: int,
                            min_span: int, k: int = 11, top_periods: int = 5,
                            min_seed_occurrences: int = 2, min_spacing_support: int = 2,
                            max_pairs_per_kmer: int = 100, backend: str = "python"
                            ) -> tuple[list[tuple[AlignmentHit, str, int]], int]:
    """Return alignment evidence, consensus and unit count for each retained locus."""
    if min_period <= 0 or max_period < min_period or min_span <= 0 or k <= 0:
        raise ValueError("Invalid elastic discovery period, span or k-mer size")
    if top_periods <= 0 or min_seed_occurrences < 2 or min_spacing_support <= 0 or max_pairs_per_kmer <= 0:
        raise ValueError("Invalid elastic discovery seed parameters")
    max_period = min(max_period, len(sequence) // 2)
    if min_period > max_period or len(sequence) < min_span:
        return [], 0
    short_periods = list(range(min_period, min(19, max_period) + 1))
    positions, overflow = extract_repeated_kmer_positions(sequence, k, min_seed_occurrences, max_pairs_per_kmer)
    histogram = build_spacing_histogram(positions, max(20, min_period), max_period, max_pairs_per_kmer)
    periods = short_periods + select_alignment_periods(histogram, top_periods, min_spacing_support)
    hits = []
    for period in periods:
        hits.extend(banded_self_align(sequence, period, min_span, backend=backend))
    result = []
    eligible = [hit for hit in hits if composition_adjusted_identity(sequence, hit) >= 0.7]
    for hit in suppress_overlapping_hits(eligible):
        consensus, units = aligned_unit_consensus(sequence, hit, backend=backend)
        if not min_period <= len(consensus) <= max_period:
            continue
        result.append((hit, consensus, units))
    return result, overflow
