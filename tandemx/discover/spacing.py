"""Bounded k-mer spacing primitives for de novo repeat discovery."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Mapping, Sequence
from math import ceil


_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BASE_CODES = {"A": 0, "C": 1, "G": 2, "T": 3}


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller forward/reverse-complement seed."""
    reverse = kmer.translate(_COMPLEMENT)[::-1]
    return min(kmer, reverse)


def is_low_complexity_kmer(kmer: str) -> bool:
    """Reject seeds with at most two bases or at least 80% of one base."""
    if not kmer:
        return True
    unique_bases = set(kmer)
    if len(unique_bases) <= 2:
        return True
    threshold = ceil(0.8 * len(kmer))
    return any(kmer.count(base) >= threshold for base in unique_bases)


def extract_repeated_kmer_positions(
    sequence: str,
    kmer_size: int,
    min_seed_occurrences: int = 2,
    max_pairs_per_kmer: int = 100,
) -> tuple[dict[int, list[int]], int]:
    """Return bounded per-read positions for repeated, non-low-complexity k-mers."""
    positions: dict[int, list[int]] = {}
    overflowed: set[int] = set()
    position_cap = max_pairs_per_kmer + 1
    sequence = sequence.upper()
    mask = (1 << (2 * kmer_size)) - 1
    reverse_shift = 2 * (kmer_size - 1)
    complexity_threshold = ceil(0.8 * kmer_size)
    forward = 0
    reverse = 0
    valid_length = 0
    base_counts = [0, 0, 0, 0]
    window_codes = [0] * kmer_size

    for index, base in enumerate(sequence):
        code = _BASE_CODES.get(base)
        if code is None:
            forward = 0
            reverse = 0
            valid_length = 0
            base_counts = [0, 0, 0, 0]
            continue

        slot = index % kmer_size
        if valid_length >= kmer_size:
            base_counts[window_codes[slot]] -= 1
        window_codes[slot] = code
        base_counts[code] += 1
        valid_length += 1
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << reverse_shift)
        if valid_length < kmer_size:
            continue
        if sum(count > 0 for count in base_counts) <= 2:
            continue
        if max(base_counts) >= complexity_threshold:
            continue

        canonical = min(forward, reverse)
        observed = positions.setdefault(canonical, [])
        if len(observed) >= position_cap:
            overflowed.add(canonical)
            continue
        observed.append(index - kmer_size + 1)
    repeated = {
        kmer: observed
        for kmer, observed in positions.items()
        if len(observed) >= min_seed_occurrences
    }
    return repeated, len(overflowed)


def build_spacing_histogram(
    repeated_positions: Mapping[object, Sequence[int]],
    min_period: int,
    max_period: int,
    max_pairs_per_kmer: int = 100,
    bin_size: int = 5,
) -> Counter[int]:
    """Count nearby seed-position spacings with a strict per-k-mer pair cap."""
    histogram: Counter[int] = Counter()
    for positions in repeated_positions.values():
        pair_count = 0
        inspected_pairs = 0
        inspection_cap = max_pairs_per_kmer * 20
        for position_gap in range(1, len(positions)):
            for left_index in range(0, len(positions) - position_gap):
                spacing = positions[left_index + position_gap] - positions[left_index]
                inspected_pairs += 1
                if min_period <= spacing <= max_period:
                    raw_bin = ((spacing + (bin_size // 2)) // bin_size) * bin_size
                    binned = min(max_period, max(min_period, raw_bin))
                    histogram[binned] += 1
                    pair_count += 1
                if pair_count >= max_pairs_per_kmer or inspected_pairs >= inspection_cap:
                    break
            if pair_count >= max_pairs_per_kmer or inspected_pairs >= inspection_cap:
                break
    return histogram


def select_candidate_periods(
    histogram: Counter[int],
    min_period: int,
    max_period: int,
    top_periods: int = 5,
    min_spacing_support: int = 2,
) -> list[int]:
    """Return the strongest supported spacing bins in deterministic order."""
    peaks = [
        (period, support)
        for period, support in histogram.items()
        if min_period <= period <= max_period and support >= min_spacing_support
    ]
    peaks.sort(key=lambda item: (-item[1], item[0]))
    return [period for period, _ in peaks[:top_periods]]


def _largest_circular_neighborhood(
    residues: list[int], period: int, tolerance: int
) -> int:
    """Return the largest neighborhood around an observed circular residue."""
    if not residues:
        return 0
    if 2 * tolerance >= period - 1:
        return len(residues)
    ordered = sorted(residues)
    doubled = ordered + [residue + period for residue in ordered]
    best = 0
    for center in ordered:
        shifted_center = center + period
        left = bisect_left(doubled, shifted_center - tolerance)
        right = bisect_right(doubled, shifted_center + tolerance)
        best = max(best, right - left)
    return best


def modulo_periodicity_score(
    repeated_positions: Mapping[object, Sequence[int]],
    period: int,
    tolerance: int = 2,
) -> float:
    """Score how consistently each repeated seed occupies one phase modulo period."""
    if period <= 0:
        return 0.0
    supported = 0
    total = 0
    for positions in repeated_positions.values():
        if len(positions) < 2:
            continue
        residues = [position % period for position in positions]
        supported += _largest_circular_neighborhood(residues, period, tolerance)
        total += len(residues)
    return supported / total if total else 0.0


def bounded_periodicity_score(
    sequence: str,
    period: int,
    max_comparisons: int = 4096,
) -> float:
    """Estimate shifted identity with a bounded number of base comparisons."""
    compared_span = len(sequence) - period
    if compared_span <= 0:
        return 0.0
    step = max(1, compared_span // max_comparisons)
    matches = 0
    valid = 0
    for index in range(0, compared_span, step):
        left = sequence[index]
        right = sequence[index + period]
        if left == "N" or right == "N":
            continue
        valid += 1
        if left == right:
            matches += 1
        if valid >= max_comparisons:
            break
    return matches / valid if valid else 0.0


def refine_candidate_period(
    sequence: str,
    repeated_positions: Mapping[object, Sequence[int]],
    candidate_periods: Sequence[int],
    min_period: int,
    max_period: int,
    refinement_radius: int = 2,
    max_seed_groups: int = 128,
) -> tuple[int, float]:
    """Refine only local variants around spacing peaks; never scan the full range."""
    local_periods = {
        candidate + offset
        for candidate in candidate_periods
        for offset in range(-refinement_radius, refinement_radius + 1)
        if min_period <= candidate + offset <= max_period
    }
    refinement_positions = dict(
        sorted(
            repeated_positions.items(),
            key=lambda item: (-len(item[1]), str(item[0])),
        )[:max_seed_groups]
    )
    best_period = 0
    best_score = 0.0
    for period in sorted(local_periods):
        seed_score = modulo_periodicity_score(refinement_positions, period)
        identity_score = bounded_periodicity_score(sequence, period)
        score = (0.20 * seed_score) + (0.80 * identity_score)
        if (score, -period) > (best_score, -best_period):
            best_period = period
            best_score = score
    return best_period, best_score
