"""Round-one bounded local re-phasing prototype; independent of public CLI."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from benchmarks.challenge.context_prototype import PeriodicContextPrototype, merge_intervals


@dataclass(frozen=True)
class LocalPhaseConfig:
    maximum_seed_gap: int = 45
    maximum_drift: int = 20
    reset_penalty: float = 1.0
    drift_penalty: float = 0.05
    predecessor_cap: int = 64
    phase_fraction: float = 0.7
    minimum_units: int = 2
    specificity_weights: bool = True
    bridge_accepted_short_gaps: bool = False


class LocalPhasePrototype(PeriodicContextPrototype):
    """Sparse chaining with a bounded change in local read/template offset.

    A phase jump has a fixed reset cost plus a per-base drift cost. Full phase
    coverage and template-coordinate recurrence are checked after traceback.
    This is a development scoring rule, not a calibrated probability model.
    """

    def __init__(self, catalogue: dict[str, str], config: LocalPhaseConfig = LocalPhaseConfig()):
        super().__init__(catalogue)
        if config.maximum_seed_gap < self.config.k or config.maximum_drift < 0:
            raise ValueError("Invalid gap/drift bound")
        if config.predecessor_cap < 1 or config.reset_penalty < 0 or config.drift_penalty < 0:
            raise ValueError("Invalid chaining costs")
        if not 0 <= config.phase_fraction <= 1 or config.minimum_units < 0:
            raise ValueError("Invalid evidence gates")
        self.local_config = config
        self.family_frequency = {
            word: len({family for family, _, _ in hits}) for word, hits in self.index.items()
        }

    def intervals(self, sequence: str) -> dict[str, list[tuple[int, int]]]:
        k = self.config.k
        cfg = self.local_config
        groups: dict[tuple[str, int], list[tuple[int, int, float]]] = defaultdict(list)
        for pos in range(len(sequence) - k + 1):
            word = sequence[pos:pos + k].upper()
            for family, strand, phase in self.index.get(word, ()):
                weight = 1 / self.family_frequency[word] if cfg.specificity_weights else 1.0
                groups[family, strand].append((pos, phase, weight))
        accepted: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for (family, _), hits in groups.items():
            period = self.periods[family]
            scores = [hit[2] for hit in hits]
            previous = [-1] * len(hits)
            progress = [0] * len(hits)
            for i, (pos, phase, weight) in enumerate(hits):
                for j in range(i - 1, max(-1, i - cfg.predecessor_cap - 1), -1):
                    distance = pos - hits[j][0]
                    if distance > cfg.maximum_seed_gap:
                        break
                    if distance <= 0:
                        continue
                    # Nearest circular displacement, then require positive
                    # template advance. This admits small insertions/deletions
                    # while rejecting arbitrary backwards phase ordering.
                    delta = (phase - hits[j][1]) % period
                    advance = delta + round((distance - delta) / period) * period
                    drift = abs(distance - advance)
                    if advance <= 0 or drift > cfg.maximum_drift:
                        continue
                    penalty = cfg.reset_penalty + cfg.drift_penalty * drift if drift else 0.0
                    candidate = scores[j] + weight - penalty
                    if candidate > scores[i]:
                        scores[i] = candidate
                        previous[i] = j
                        progress[i] = progress[j] + advance
            # Trace terminal paths only, preserving separate interrupted blocks.
            parents = {p for p in previous if p >= 0}
            for tip in range(len(hits)):
                if tip in parents or progress[tip] + k < cfg.minimum_units * period:
                    continue
                path = []
                cursor = tip
                while cursor >= 0:
                    path.append(cursor)
                    cursor = previous[cursor]
                if len({hits[i][1] for i in path}) / period < cfg.phase_fraction:
                    continue
                # Do not fill seed-free inserted sequence as repeat bp.
                pieces = [(hits[i][0], hits[i][0] + k) for i in path]
                if cfg.bridge_accepted_short_gaps:
                    pieces = [(min(start for start, _ in pieces), max(end for _, end in pieces))]
                accepted[family].extend(merge_intervals(pieces))
        return {family: merge_intervals(rows) for family, rows in accepted.items()}
