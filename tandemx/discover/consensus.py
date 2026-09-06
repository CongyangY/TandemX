"""Deterministic majority consensus from alignment-derived repeat units."""
from __future__ import annotations

from bisect import bisect_left
from collections import Counter
from statistics import median_low

from tandemx.discover.alignment import AlignmentHit, global_align_ops


def extract_aligned_units(sequence: str, hit: AlignmentHit, max_units: int = 32) -> list[str]:
    """Follow the adjacent-copy correspondence to cut homologous unit boundaries.

    Up to 32 evenly spaced complete units are retained for consensus; read-local
    alignment still uses the entire region. Terminal partial units are omitted.
    """
    if max_units < 2:
        raise ValueError("max_units must be at least two")
    mapping = dict(hit.pairs)
    keys = sorted(mapping)
    if not keys:
        return []
    boundaries = [keys[0]]
    while boundaries[-1] < hit.end:
        start = boundaries[-1]
        index = bisect_left(keys, start)
        if index == len(keys) or keys[index] - start > 5:
            break
        anchor = keys[index]
        end = mapping[anchor] - (anchor - start)
        if not 0.7 * hit.period <= end - start <= 1.3 * hit.period or end > hit.end:
            break
        boundaries.append(end)
    if 0.9 * hit.period <= hit.end - boundaries[-1] <= 1.1 * hit.period:
        boundaries.append(hit.end)
    units = [sequence[left:right] for left, right in zip(boundaries, boundaries[1:])]
    if len(units) > max_units:
        indices = [i * (len(units) - 1) // (max_units - 1) for i in range(max_units)]
        units = [units[i] for i in indices]
    return units


def _vote(reference: str, units: list[str], backend: str) -> str:
    columns = [Counter() for _ in reference]
    insertions = [Counter() for _ in range(len(reference) + 1)]
    for unit in units:
        ops = global_align_ops(reference, unit, backend=backend)
        r = q = 0
        inserted = ["" for _ in insertions]
        for op in ops:
            if op == "M":
                columns[r][unit[q]] += 1; r += 1; q += 1
            elif op == "D":
                columns[r]["-"] += 1; r += 1
            else:
                inserted[r] += unit[q]; q += 1
        for position, word in enumerate(inserted):
            insertions[position][word] += 1
    # Ties retain a reference base and prefer no insertion. This conservatism is
    # explicitly a majority consensus, not an inferred ancestral monomer.
    result = []
    for position in range(len(reference) + 1):
        insertion, support = min(insertions[position].items(), key=lambda item: (-item[1], len(item[0]), item[0]))
        if support * 2 > len(units):
            result.append(insertion)
        if position < len(reference):
            base = min(columns[position], key=lambda b: (-columns[position][b], b != reference[position], b))
            if base != "-":
                result.append(base)
    return "".join(result)


def aligned_unit_consensus(sequence: str, hit: AlignmentHit, *, backend: str = "python") -> tuple[str, int]:
    units = extract_aligned_units(sequence, hit)
    if len(units) < 2:
        return "", len(units)
    length = median_low([len(unit) for unit in units])
    # A terminal unit may include unrelated flank sequence. Choose a central
    # length with high cross-unit seed support instead of privileging the first
    # unit. Sets give each observed unit one vote per seed, not abundance votes.
    k = min(7, length)
    sketches = [set(unit[i:i + k] for i in range(len(unit) - k + 1)) for unit in units]
    support = Counter(word for sketch in sketches for word in sketch)
    reference = min(enumerate(units), key=lambda item: (
        abs(len(item[1]) - length),
        -sum(support[word] - 1 for word in sketches[item[0]]) / max(1, len(sketches[item[0]])),
        item[0],
    ))[1]
    for _ in range(2):
        updated = _vote(reference, units, backend)
        if not updated or updated == reference:
            break
        reference = updated
    return reference, len(units)
