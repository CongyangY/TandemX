"""Shared canonical k-mer primitives with bounded-memory rolling scans."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator


_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_BYTE_CODES = [-1] * 256
for _base, _code in (
    (ord("A"), 0),
    (ord("C"), 1),
    (ord("G"), 2),
    (ord("T"), 3),
    (ord("a"), 0),
    (ord("c"), 1),
    (ord("g"), 2),
    (ord("t"), 3),
):
    _BYTE_CODES[_base] = _code
_BYTE_CODES = tuple(_BYTE_CODES)


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def canonical_kmer(kmer: str) -> str:
    normalized = kmer.upper()
    return min(normalized, reverse_complement(normalized))


def is_low_complexity_kmer(kmer: str) -> bool:
    """Flag single-base and simple alternating k-mers without rejecting all two-base seeds."""
    normalized = kmer.upper()
    if not normalized:
        return True
    counts = Counter(normalized)
    if max(counts.values()) / len(normalized) >= 0.8:
        return True
    if len(normalized) >= 4:
        dinucleotides = {
            normalized[index : index + 2]
            for index in range(0, len(normalized) - 1, 2)
        }
        if len(dinucleotides) <= 1:
            return True
    return False


def iter_linear_canonical_kmers(sequence: str, k: int) -> Iterator[str]:
    if k <= 0:
        raise ValueError("k must be positive")
    normalized = sequence.upper()
    for index in range(0, len(normalized) - k + 1):
        kmer = normalized[index : index + k]
        if all(base in "ACGT" for base in kmer):
            yield canonical_kmer(kmer)


def iter_circular_canonical_kmers(sequence: str, k: int) -> Iterator[str]:
    """Yield one tandem-context k-mer for every monomer start phase."""
    if k <= 0:
        raise ValueError("k must be positive")
    normalized = sequence.upper()
    if not normalized:
        return
    repeat_count = (len(normalized) + k - 2) // len(normalized) + 1
    circular = (normalized * repeat_count)[: len(normalized) + k - 1]
    for index in range(len(normalized)):
        kmer = circular[index : index + k]
        if all(base in "ACGT" for base in kmer):
            yield canonical_kmer(kmer)


def circular_kmer_counts(sequence: str, k: int) -> Counter[str]:
    return Counter(iter_circular_canonical_kmers(sequence, k))


def canonical_kmer_code(kmer: str) -> int:
    """Encode an A/C/G/T k-mer using the same canonical 2-bit representation as scans."""
    normalized = kmer.upper()
    if not 1 <= len(normalized) <= 31:
        raise ValueError("encoded canonical k-mers require length in 1..=31")
    forward = 0
    reverse = 0
    for index, base in enumerate(normalized.encode("ascii")):
        code = _BYTE_CODES[base]
        if code < 0:
            raise ValueError("canonical k-mer encoding requires only A/C/G/T")
        forward = (forward << 2) | code
        reverse |= (3 - code) << (2 * index)
    return min(forward, reverse)


def iter_canonical_kmer_codes(
    sequence: str,
    k: int,
    *,
    filter_low_complexity: bool = False,
) -> Iterator[tuple[int, int]]:
    """Yield ``(0-based position, canonical 2-bit code)`` in one sequence pass."""
    if not 1 <= k <= 31:
        raise ValueError("encoded canonical k-mers require k in 1..=31")
    mask = (1 << (2 * k)) - 1
    reverse_shift = 2 * (k - 1)
    forward = 0
    reverse = 0
    valid_length = 0
    base_counts = [0, 0, 0, 0]
    window_codes = [0] * k
    distinct_bases = 0
    complexity_threshold = (4 * k + 4) // 5

    for index, base in enumerate(sequence.encode("ascii")):
        code = _BYTE_CODES[base]
        if code < 0:
            forward = 0
            reverse = 0
            valid_length = 0
            base_counts = [0, 0, 0, 0]
            distinct_bases = 0
            continue
        slot = index % k
        if valid_length >= k:
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
        if valid_length < k:
            continue
        if filter_low_complexity and (
            max(base_counts) >= complexity_threshold
            or (distinct_bases <= 2 and _is_simple_two_base_window(window_codes, index, k))
        ):
            continue
        yield index + 1 - k, min(forward, reverse)


def _is_simple_two_base_window(window_codes: list[int], index: int, k: int) -> bool:
    """Reject only exact homopolymer/dinucleotide patterns, not all two-base biology."""
    ordered = [window_codes[(index + 1 + offset) % k] for offset in range(k)]
    if len(set(ordered)) <= 1:
        return True
    return k >= 4 and all(ordered[offset] == ordered[offset % 2] for offset in range(2, k))
