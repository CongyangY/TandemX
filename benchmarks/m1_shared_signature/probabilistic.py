"""Development-only read likelihood mixture with explicit unknown rejection.

This model sees catalogue sequences and reads only. It estimates each read's
substitution error from its best circular alignment, fits catalogue-group and
unknown mixing weights by EM, and assigns only posterior-supported reads.
Exact cyclic-equivalent units remain an unresolved group.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .model import DNA, rotations


@dataclass(frozen=True)
class Settings:
    posterior_min: float = 0.95
    max_mismatch: int = 18
    error_floor: float = 0.02
    error_ceiling: float = 0.30
    max_iterations: int = 100
    tolerance: float = 1e-9


def canonical_rotation(sequence: str) -> str:
    return min(sequence[i:] + sequence[:i] for i in range(len(sequence)))


def distance_matrix(reads: list[str], catalogue: dict[str, str]
                    ) -> tuple[np.ndarray, list[list[str]]]:
    if not catalogue:
        raise ValueError("Need nonempty catalogue")
    lengths = {len(seq) for seq in catalogue.values()}
    if len(lengths) != 1 or next(iter(lengths)) == 0:
        raise ValueError("This prototype requires nonempty equal-length units")
    groups_by_rotation: dict[str, list[str]] = {}
    for name, sequence in catalogue.items():
        if set(sequence) - set(DNA):
            raise ValueError("Catalogue must contain only A/C/G/T")
        groups_by_rotation.setdefault(canonical_rotation(sequence), []).append(name)
    groups = list(groups_by_rotation.values())
    references = [rotations(catalogue[names[0]]) for names in groups]
    distances = np.empty((len(reads), len(groups)), dtype=float)
    for row, read in enumerate(reads):
        if len(read) not in lengths or set(read) - set(DNA):
            raise ValueError("Reads must be equal-length A/C/G/T sequences")
        query = np.fromiter((DNA.index(base) for base in read), dtype=np.uint8)
        distances[row] = [np.min(np.count_nonzero(reference != query, axis=1))
                          for reference in references]
    return distances, groups


def infer(reads: list[str], catalogue: dict[str, str], settings: Settings = Settings()) -> dict:
    if not 0.5 < settings.posterior_min < 1 or settings.max_mismatch < 0 or not (
        0 < settings.error_floor < settings.error_ceiling < 0.75
    ) or settings.max_iterations < 1:
        raise ValueError("Invalid probabilistic development settings")
    distances, groups = distance_matrix(reads, catalogue)
    names = ["+".join(sorted(group)) for group in groups]
    family = {name: 0.0 if len(group) == 1 else None
              for group in groups for name in group}
    group_estimates = {name: 0.0 for name, group in zip(names, groups) if len(group) > 1}
    if not reads:
        return dict(family=family, groups=group_estimates, rejected_unknown=0,
                    rejected_ambiguous=0, rejected_gate=0, assigned=0,
                    read_assignment=[],
                    posterior_soft={name: 0.0 for name in names},
                    mixture_weights={name: 0.0 for name in names}, iterations=0)

    length = len(reads[0])
    minimum = np.min(distances, axis=1)
    estimated_error = np.clip(minimum / length, settings.error_floor, settings.error_ceiling)
    log_likelihood = ((length - distances) * np.log1p(-estimated_error[:, None])
                      + distances * np.log(estimated_error[:, None] / 3))
    # The unknown model is uniform DNA; it has no catalogue-specific exposure.
    unknown_log_likelihood = np.full((len(reads), 1), -length * math.log(4))
    log_likelihood = np.hstack((log_likelihood, unknown_log_likelihood))
    weights = np.full(len(groups) + 1, 1.0 / (len(groups) + 1))
    iterations = 0
    for iterations in range(1, settings.max_iterations + 1):
        joint = log_likelihood + np.log(np.maximum(weights, 1e-12))
        joint -= np.max(joint, axis=1, keepdims=True)
        posterior = np.exp(joint)
        posterior /= np.sum(posterior, axis=1, keepdims=True)
        updated = np.maximum(np.mean(posterior, axis=0), 1e-12)
        updated /= np.sum(updated)
        if np.max(np.abs(updated - weights)) < settings.tolerance:
            weights = updated
            break
        weights = updated
    joint = log_likelihood + np.log(np.maximum(weights, 1e-12))
    joint -= np.max(joint, axis=1, keepdims=True)
    posterior = np.exp(joint)
    posterior /= np.sum(posterior, axis=1, keepdims=True)

    winner = np.argmax(posterior, axis=1)
    confidence = np.max(posterior, axis=1)
    beyond_gate = minimum > settings.max_mismatch
    unknown = (~beyond_gate) & (winner == len(groups))
    ambiguous = (~beyond_gate) & (~unknown) & (confidence < settings.posterior_min)
    assignable = (~beyond_gate) & (~unknown) & (~ambiguous)
    for index, (name, group) in enumerate(zip(names, groups)):
        assigned = float(np.sum(assignable & (winner == index)))
        if len(group) == 1:
            family[group[0]] = assigned
        else:
            group_estimates[name] = assigned
    read_assignment = ["gate" if beyond_gate[i] else
                       "unknown" if unknown[i] else
                       "ambiguous" if ambiguous[i] else names[winner[i]]
                       for i in range(len(reads))]
    return dict(family=family, groups=group_estimates,
                rejected_unknown=int(np.sum(unknown)),
                rejected_ambiguous=int(np.sum(ambiguous)),
                rejected_gate=int(np.sum(beyond_gate)),
                assigned=int(np.sum(assignable)),
                read_assignment=read_assignment,
                posterior_soft={name: float(np.sum(posterior[:, index]))
                                for index, name in enumerate(names)},
                mixture_weights={name: float(weight) for name, weight in
                                 zip([*names, "unknown"], weights)},
                iterations=iterations)
