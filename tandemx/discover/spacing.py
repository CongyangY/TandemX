"""Bounded k-mer spacing primitives for de novo repeat discovery."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Mapping, Sequence
from math import ceil, isqrt


_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BYTE_CODES = [-1] * 256
for _base, _code in ((ord("A"), 0), (ord("C"), 1), (ord("G"), 2), (ord("T"), 3)):
    _BYTE_CODES[_base] = _code
_BYTE_CODES = tuple(_BYTE_CODES)


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
    first_positions: dict[int, int] = {}
    repeated_positions: dict[int, list[int]] = {}
    overflowed: set[int] = set()
    position_cap = max_pairs_per_kmer + 1
    encoded_sequence = sequence.upper().encode("ascii")
    mask = (1 << (2 * kmer_size)) - 1
    reverse_shift = 2 * (kmer_size - 1)
    complexity_threshold = ceil(0.8 * kmer_size)
    forward = 0
    reverse = 0
    valid_length = 0
    base_counts = [0, 0, 0, 0]
    distinct_bases = 0
    window_codes = [0] * kmer_size

    for index, base in enumerate(encoded_sequence):
        code = _BYTE_CODES[base]
        if code < 0:
            forward = 0
            reverse = 0
            valid_length = 0
            base_counts = [0, 0, 0, 0]
            distinct_bases = 0
            continue

        slot = index % kmer_size
        if valid_length >= kmer_size:
            outgoing = window_codes[slot]
            base_counts[outgoing] -= 1
            if base_counts[outgoing] == 0:
                distinct_bases -= 1
        window_codes[slot] = code
        if base_counts[code] == 0:
            distinct_bases += 1
        base_counts[code] += 1
        valid_length += 1
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << reverse_shift)
        if valid_length < kmer_size:
            continue
        if max(base_counts) >= complexity_threshold:
            continue
        if distinct_bases <= 2 and _is_simple_periodic_window(window_codes, index, kmer_size):
            continue

        canonical = forward if forward < reverse else reverse
        position = index - kmer_size + 1
        observed = repeated_positions.get(canonical)
        if observed is None:
            first_position = first_positions.get(canonical)
            if first_position is None:
                first_positions[canonical] = position
                continue
            first_positions.pop(canonical)
            observed = [first_position, position]
            repeated_positions[canonical] = observed
            continue
        if len(observed) >= position_cap:
            overflowed.add(canonical)
            continue
        observed.append(position)

    if min_seed_occurrences > 2:
        repeated_positions = {
            kmer: observed
            for kmer, observed in repeated_positions.items()
            if len(observed) >= min_seed_occurrences
        }
    return repeated_positions, len(overflowed)


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
    max_comparisons: int = 1024,
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


def best_local_periodicity_score(
    sequence: str,
    period: int,
    min_repeat_span: int,
    acceptance_score: float = 0.75,
) -> tuple[float, int, int]:
    """Find the strongest local interval for one period in linear time.

    Match/mismatch weights make an interval profitable only above the requested
    identity threshold. The minimum comparison span enforces the biological
    repeat-span boundary instead of accepting short chance matches in flanks.
    """
    compared_span = len(sequence) - period
    # Require evidence spanning at least two repeat units. Otherwise a long
    # candidate can win on one accidental matching base when min_repeat_span
    # is shorter than the candidate period.
    min_compared = max(1, period, min_repeat_span - period)
    if period <= 0 or compared_span < min_compared:
        return 0.0, 0, 0

    mismatch_penalty = acceptance_score / max(1e-12, 1.0 - acceptance_score)
    prefix_score = [0.0]
    prefix_matches = [0]
    prefix_valid = [0]
    for index in range(compared_span):
        left = sequence[index]
        right = sequence[index + period]
        if left == "N" or right == "N":
            value = 0.0
            match = 0
            valid = 0
        elif left == right:
            value = 1.0
            match = 1
            valid = 1
        else:
            value = -mismatch_penalty
            match = 0
            valid = 1
        prefix_score.append(prefix_score[-1] + value)
        prefix_matches.append(prefix_matches[-1] + match)
        prefix_valid.append(prefix_valid[-1] + valid)

    best_key = (float("-inf"), 0, 0)
    best_start = 0
    best_end = 0
    minimum_prefix_value = prefix_score[0]
    minimum_prefix_index = 0
    for end in range(min_compared, compared_span + 1):
        eligible = end - min_compared
        eligible_value = prefix_score[eligible]
        if eligible_value < minimum_prefix_value:
            minimum_prefix_value = eligible_value
            minimum_prefix_index = eligible
        interval_score = prefix_score[end] - minimum_prefix_value
        interval_length = end - minimum_prefix_index
        key = (interval_score, interval_length, -minimum_prefix_index)
        if key > best_key:
            best_key = key
            best_start = minimum_prefix_index
            best_end = end

    valid = prefix_valid[best_end] - prefix_valid[best_start]
    matches = prefix_matches[best_end] - prefix_matches[best_start]
    identity = matches / valid if valid else 0.0
    repeat_end = min(len(sequence), best_end + period)
    if repeat_end - best_start < min_repeat_span or identity < acceptance_score:
        return identity, best_start, repeat_end
    return identity, best_start, repeat_end


def expand_candidate_periods(
    candidate_periods: Sequence[int],
    min_period: int,
    max_period: int,
    refinement_radius: int = 2,
) -> list[int]:
    """Add local refinements and plausible fundamental divisors of spacing peaks."""
    expanded: set[int] = set()
    for candidate in candidate_periods:
        for offset in range(-refinement_radius, refinement_radius + 1):
            local = candidate + offset
            if not min_period <= local <= max_period:
                continue
            expanded.add(local)
            for divisor in range(2, isqrt(local) + 1):
                if local % divisor:
                    continue
                for fundamental in (divisor, local // divisor):
                    if min_period <= fundamental <= max_period:
                        expanded.add(fundamental)
    return sorted(expanded)


def refine_candidate_period(
    sequence: str,
    repeated_positions: Mapping[object, Sequence[int]],
    candidate_periods: Sequence[int],
    min_period: int,
    max_period: int,
    refinement_radius: int = 2,
    max_seed_groups: int = 128,
    acceptance_score: float = 0.75,
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
            key=lambda item: (-len(item[1]), _stable_seed_key(item[0])),
        )[:max_seed_groups]
    )
    best_period = 0
    best_score = 0.0
    seed_weight = 0.20
    identity_weight = 0.80
    minimum_viable_identity = max(
        0.0,
        (acceptance_score - seed_weight) / identity_weight,
    )
    for period in sorted(local_periods):
        identity_score = bounded_periodicity_score(sequence, period)
        if identity_score < minimum_viable_identity:
            score = identity_weight * identity_score
        else:
            seed_score = modulo_periodicity_score(refinement_positions, period)
            score = (seed_weight * seed_score) + (identity_weight * identity_score)
        if (score, -period) > (best_score, -best_period):
            best_period = period
            best_score = score
    return best_period, best_score


def refine_candidate_period_with_interval(
    sequence: str,
    repeated_positions: Mapping[object, Sequence[int]],
    candidate_periods: Sequence[int],
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    refinement_radius: int = 2,
    max_seed_groups: int = 128,
    acceptance_score: float = 0.75,
) -> tuple[int, float, int, int]:
    """Refine period and return the best local 0-based half-open repeat interval."""
    local_periods = expand_candidate_periods(
        candidate_periods,
        min_period,
        max_period,
        refinement_radius,
    )
    refinement_positions = dict(
        sorted(
            repeated_positions.items(),
            key=lambda item: (-len(item[1]), _stable_seed_key(item[0])),
        )[:max_seed_groups]
    )
    best_period = 0
    best_score = 0.0
    best_start = 0
    best_end = 0
    seed_weight = 0.20
    identity_weight = 0.80
    minimum_viable_identity = max(
        0.0,
        (acceptance_score - seed_weight) / identity_weight,
    )
    for period in local_periods:
        identity_score, start, end = best_local_periodicity_score(
            sequence,
            period,
            min_repeat_span,
            acceptance_score=minimum_viable_identity,
        )
        if identity_score < minimum_viable_identity:
            score = identity_weight * identity_score
        else:
            seed_score = modulo_periodicity_score(refinement_positions, period)
            score = (seed_weight * seed_score) + (identity_weight * identity_score)
        interval_length = end - start
        best_key = (best_score, best_end - best_start, -best_period)
        key = (score, interval_length, -period)
        if key > best_key:
            best_period = period
            best_score = score
            best_start = start
            best_end = end
    return best_period, best_score, best_start, best_end


def _stable_seed_key(seed: object) -> tuple[int, int | str]:
    """Match Rust's numeric seed ordering while supporting string-keyed tests."""
    if isinstance(seed, int):
        return 0, seed
    return 1, str(seed)


def _is_simple_periodic_window(window_codes: list[int], index: int, k: int) -> bool:
    ordered = [window_codes[(index + 1 + offset) % k] for offset in range(k)]
    if len(set(ordered)) <= 1:
        return True
    return k >= 4 and all(ordered[offset] == ordered[offset % 2] for offset in range(2, k))
