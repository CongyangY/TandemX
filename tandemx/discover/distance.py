"""Thresholded, sequence-based comparisons for operational monomer clusters."""
from __future__ import annotations

from collections import Counter, defaultdict
from math import floor


def circular_words(sequence: str, k: int) -> Counter[str]:
    extended = sequence + sequence[:k - 1]
    return Counter(extended[i:i + k] for i in range(len(sequence)))


def bounded_edit_distance(a: str, b: str, limit: int, backend: str = "python") -> int:
    """Exact global unit-cost distance if <=limit, otherwise limit+1; N mismatches."""
    if limit < 0 or not (a + b).isascii():
        raise ValueError("Need nonnegative edit limit and ASCII sequences")
    if max(1, len(a)) * (2 * limit + 1) > 32_000_000:
        raise ValueError("Cluster comparison exceeds 32 million cells")
    a, b = a.upper(), b.upper()
    if backend == "rust":
        from tandemx.discover.rust_backend import RustBackendUnavailable
        try:
            from tandemx import _rust_core
        except ImportError as exc:
            raise RustBackendUnavailable("Rebuild the Rust extension for sequence clustering") from exc
        if not hasattr(_rust_core, "bounded_edit_distance"):
            raise RustBackendUnavailable("Rust extension lacks sequence clustering; reinstall TandemX")
        return _rust_core.bounded_edit_distance(a, b, limit)
    if backend != "python":
        raise ValueError("Distance backend must be python or rust")
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    start = 0
    while start < min(len(a), len(b)) and a[start] == b[start] and a[start] in "ACGT":
        start += 1
    a, b = a[start:], b[start:]
    while a and b and a[-1] == b[-1] and a[-1] in "ACGT":
        a, b = a[:-1], b[:-1]
    if not a or not b:
        return min(limit + 1, max(len(a), len(b)))
    previous = {j: j for j in range(min(len(b), limit) + 1)}
    for i, x in enumerate(a, 1):
        current = {}
        for j in range(max(0, i - limit), min(len(b), i + limit) + 1):
            current[j] = (i if j == 0 else min(previous.get(j, limit + 1) + 1,
                         current.get(j - 1, limit + 1) + 1,
                         previous.get(j - 1, limit + 1) + int(x != b[j - 1] or x not in "ACGT")))
        if min(current.values(), default=limit + 1) > limit:
            return limit + 1
        previous = current
    return min(limit + 1, previous.get(len(b), limit + 1))


def cyclic_merge_evidence(reference: str, query: str, minimum_identity: float,
                          backend: str = "python") -> tuple[int, float] | None:
    """Return a witnessed distance upper bound, never a false exact-minimum claim.

    Shared q-grams provide only a necessary-condition rejection. Seed-supported
    rotations are tried first; all remaining rotations/strands are tried before
    rejecting a pair that passed the q-gram bound. An accepted alignment may not
    be the minimum-distance alignment, so its similarity is a lower bound.
    """
    if not reference or not query or not 0 < minimum_identity <= 1:
        raise ValueError("Need nonempty monomers and cluster identity in (0,1]")
    reference, query = reference.upper(), query.upper()
    if set(reference + query) - set("ACGTN"):
        raise ValueError("Monomers must contain only ACGTN")
    size = max(len(reference), len(query))
    limit = floor((1 - minimum_identity) * size + 1e-9)
    if abs(len(reference) - len(query)) > limit:
        return None
    k = min(9, len(reference), len(query))
    words_a = circular_words(reference, k)
    reverse = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for strand in (query, reverse):
        words_b = circular_words(strand, k)
        # At most k circular q-grams per unit edit are destroyed/created.
        if sum((words_a & words_b).values()) < max(0, size - k * limit):
            continue
        positions: dict[str, list[int]] = defaultdict(list)
        extended = strand + strand[:k - 1]
        for j in range(len(strand)):
            word = extended[j:j + k]
            if len(positions[word]) < 32:
                positions[word].append(j)
        votes: Counter[int] = Counter()
        extended_a = reference + reference[:k - 1]
        for i in range(0, len(reference), max(1, len(reference) // 32)):
            for j in positions.get(extended_a[i:i + k], []):
                votes[(j - i) % len(strand)] += 1
        preferred = [offset for offset, _ in sorted(votes.items(), key=lambda item: (-item[1], item[0]))]
        tried = set(preferred)
        for offset in [*preferred, *(i for i in range(len(strand)) if i not in tried)]:
            rotated = strand[offset:] + strand[:offset]
            distance = bounded_edit_distance(reference, rotated, limit, backend)
            if distance <= limit:
                return distance, 1 - distance / size
    return None
