"""Development-only periodic-context gate; not a production abundance estimator.

Exact seeds vote for read/template phase diagonals. A bounded-gap diagonal
must span two units and cover enough distinct template start phases. Indels
split diagonals: sensitivity loss is an explicit experimental endpoint.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from tandemx.utils.kmers import reverse_complement


@dataclass(frozen=True)
class ContextConfig:
    k: int = 11
    minimum_units: int = 2
    phase_fraction: float = 0.7
    maximum_seed_gap: int = 22


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("Invalid half-open interval")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class PeriodicContextPrototype:
    """Frozen-catalogue exact-phase experiment with per-read working memory."""

    def __init__(self, catalogue: dict[str, str], config: ContextConfig = ContextConfig()):
        if config.k < 1 or config.minimum_units < 2:
            raise ValueError("Positive k and at least two units required")
        if not 0 < config.phase_fraction <= 1 or config.maximum_seed_gap < config.k:
            raise ValueError("Invalid phase fraction or seed gap")
        self.config = config
        self.periods: dict[str, int] = {}
        self.index: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        if not catalogue:
            raise ValueError("Catalogue is empty")
        for family, sequence in catalogue.items():
            sequence = sequence.upper()
            if not family or len(sequence) < config.k or set(sequence) - set("ACGT"):
                raise ValueError("Catalogue requires named ACGT units of length >= k")
            period = len(sequence)
            self.periods[family] = period
            for strand, unit in enumerate((sequence, reverse_complement(sequence))):
                circular = unit + unit[:config.k - 1]
                for phase in range(period):
                    self.index[circular[phase:phase + config.k]].append((family, strand, phase))

    def intervals(self, sequence: str) -> dict[str, list[tuple[int, int]]]:
        """Return context-qualified spans, retaining cross-family overlaps."""
        sequence = sequence.upper()
        k = self.config.k
        diagonals: dict[tuple[str, int, int], list[tuple[int, int]]] = defaultdict(list)
        for pos in range(len(sequence) - k + 1):
            for family, strand, phase in self.index.get(sequence[pos:pos + k], ()):
                diagonal = (pos - phase) % self.periods[family]
                diagonals[family, strand, diagonal].append((pos, phase))
        accepted: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for (family, _, _), hits in diagonals.items():
            period = self.periods[family]
            start = last = hits[0][0]
            phases = {hits[0][1]}
            for pos, phase in hits[1:] + [(len(sequence) + self.config.maximum_seed_gap + 1, -1)]:
                if pos - last > self.config.maximum_seed_gap:
                    end = last + k
                    if (end - start >= self.config.minimum_units * period
                            and len(phases) / period >= self.config.phase_fraction):
                        accepted[family].append((start, end))
                    start = pos
                    phases = set()
                phases.add(phase)
                last = pos
        return {family: merge_intervals(rows) for family, rows in accepted.items()}


def partition_evidence(
    intervals: dict[str, list[tuple[int, int]]], read_length: int,
) -> dict[str, object]:
    """Count unique assignments and overlapping ambiguity without duplicating bp.

    Eligible bp is a per-family potential count, not an additive total or a
    statistical upper confidence bound. Unassigned includes undetected repeats.
    """
    events: dict[int, list[tuple[str, int]]] = defaultdict(list)
    unique = {family: 0 for family in intervals}
    eligible = {family: 0 for family in intervals}
    for family, rows in intervals.items():
        for start, end in merge_intervals(rows):
            if end > read_length:
                raise ValueError("Interval exceeds read length")
            events[start].append((family, 1))
            events[end].append((family, -1))
    active: set[str] = set()
    last = assigned = ambiguous = 0
    for pos, changes in sorted(events.items()):
        width = pos - last
        for family in active:
            eligible[family] += width
        if len(active) == 1:
            unique[next(iter(active))] += width
            assigned += width
        elif len(active) > 1:
            ambiguous += width
        for family, delta in changes:
            if delta == 1:
                active.add(family)
            else:
                active.remove(family)
        last = pos
    return {"unique_bp": unique, "eligible_bp": eligible,
            "ambiguous_bp": ambiguous, "unassigned_bp": read_length - assigned - ambiguous}
