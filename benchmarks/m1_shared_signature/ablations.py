"""Fixed development-only A/C/D abundance candidates on a shared catalogue.

A: exclusive 5-mer depth with a fixed error-survival correction.
C: constrained Poisson count regression on an error-aware 5-mer design.
D: per-read composite-multinomial mixture EM on the same design.

No candidate uses simulator source labels. The common ordinary-mapping gate is
applied before calling these functions. This is bounded to <=6 family groups.
"""
from __future__ import annotations

from collections import Counter
from itertools import product

import numpy as np

from .model import DNA, circular_kmers

K = 5
ASSUMED_SUBSTITUTION = 0.08
MAX_ITERATIONS = 200
POSTERIOR_MIN = 0.95
MIN_SINGULAR_RATIO = 0.05


def design(catalogue: dict[str, str], k: int = K, error: float = ASSUMED_SUBSTITUTION
           ) -> tuple[list[list[str]], list[str], np.ndarray, list[Counter[str]]]:
    if not catalogue or len(catalogue) > 6 or not 0 < error < 0.75:
        raise ValueError("Need 1-6 catalogue entries and a valid error probability")
    lengths = {len(sequence) for sequence in catalogue.values()}
    if len(lengths) != 1 or next(iter(lengths)) < k:
        raise ValueError("This toy requires equal-length units at least k bases long")
    signatures: dict[tuple[tuple[str, int], ...], list[str]] = {}
    for name, sequence in catalogue.items():
        if set(sequence) - set(DNA):
            raise ValueError("Catalogue must contain only A/C/G/T")
        signature = circular_kmers(sequence, k)
        signatures.setdefault(tuple(sorted(signature.items())), []).append(name)
    groups = list(signatures.values())
    words = ["".join(parts) for parts in product(DNA, repeat=k)]
    profiles = [circular_kmers(catalogue[group[0]], k) for group in groups]
    expected = np.zeros((len(words), len(groups)), dtype=float)
    for column, profile in enumerate(profiles):
        for template, multiplicity in profile.items():
            mismatch = np.fromiter((sum(a != b for a, b in zip(template, word))
                                    for word in words), dtype=np.int16)
            expected[:, column] += multiplicity * (
                (1 - error) ** (k - mismatch) * (error / 3) ** mismatch)
    return groups, words, expected, profiles


def observations(reads: list[str], words: list[str], k: int = K) -> np.ndarray:
    location = {word: i for i, word in enumerate(words)}
    matrix = np.zeros((len(reads), len(words)), dtype=float)
    for row, read in enumerate(reads):
        if len(read) < k or set(read) - set(DNA):
            raise ValueError("Reads must be valid A/C/G/T units")
        for word, count in circular_kmers(read, k).items():
            matrix[row, location[word]] = count
    return matrix


def expand_estimates(groups: list[list[str]], counts: np.ndarray
                     ) -> tuple[dict[str, float | None], dict[str, float]]:
    families: dict[str, float | None] = {}
    unresolved: dict[str, float] = {}
    for group, count in zip(groups, counts):
        if len(group) == 1:
            families[group[0]] = float(count)
        else:
            for name in group:
                families[name] = None
            unresolved["+".join(sorted(group))] = float(count)
    return families, unresolved


def design_refusal(expected: np.ndarray, catalogue: dict[str, str], n_reads: int
                   ) -> dict | None:
    """Refuse all individual counts when the input design is nearly rank deficient."""
    singular = np.linalg.svd(expected, compute_uv=False)
    if singular[-1] / singular[0] >= MIN_SINGULAR_RATIO:
        return None
    names = sorted(catalogue)
    return dict(family={name: None for name in names},
                groups={"+".join(names): float(n_reads)},
                design_refusal="singular_ratio_below_0.05")


def discriminative(reads: list[str], catalogue: dict[str, str]
                   ) -> dict:
    groups, words, expected, profiles = design(catalogue)
    refused = design_refusal(expected, catalogue, len(reads))
    if refused is not None:
        return dict(**refused, exclusive_kmers_per_read=[], total_estimated_reads=float(len(reads)))
    matrix = observations(reads, words)
    total = np.sum(matrix, axis=0)
    location = {word: i for i, word in enumerate(words)}
    estimates = []
    exclusive_counts = []
    for index, profile in enumerate(profiles):
        other_words = set().union(*(set(candidate) for j, candidate in enumerate(profiles)
                                    if j != index))
        exclusive = [word for word in profile if word not in other_words]
        expected_per_read = sum(profile[word] for word in exclusive)
        exclusive_counts.append(expected_per_read)
        if expected_per_read == 0:
            estimates.append(0.0)
        else:
            observed = sum(total[location[word]] for word in exclusive)
            estimates.append(observed / (expected_per_read * (1 - ASSUMED_SUBSTITUTION) ** K))
    families, unresolved = expand_estimates(groups, np.array(estimates))
    # A group with no exclusive signature cannot support an individual number.
    for group, count in zip(groups, exclusive_counts):
        if count == 0 and len(group) == 1:
            families[group[0]] = None
    return dict(family=families, groups=unresolved, exclusive_kmers_per_read=exclusive_counts,
                total_estimated_reads=float(sum(estimates)))


def poisson_counts(reads: list[str], catalogue: dict[str, str]) -> dict:
    groups, words, expected, _ = design(catalogue)
    refused = design_refusal(expected, catalogue, len(reads))
    if refused is not None:
        return dict(**refused, iterations=0, total_estimated_reads=float(len(reads)))
    matrix = observations(reads, words)
    n = len(reads)
    weights = np.full(len(groups), n / len(groups), dtype=float)
    iterations = 0
    if n:
        observed = np.sum(matrix, axis=0)
        for iterations in range(1, MAX_ITERATIONS + 1):
            predicted = np.maximum(expected @ weights, 1e-12)
            update = weights * (expected.T @ (observed / predicted)) / np.sum(expected, axis=0)
            update *= n / np.sum(update)
            if np.max(np.abs(update - weights)) < 1e-8:
                weights = update
                break
            weights = update
    families, unresolved = expand_estimates(groups, weights)
    return dict(family=families, groups=unresolved, iterations=iterations,
                total_estimated_reads=float(sum(weights)))


def em_mixture(reads: list[str], catalogue: dict[str, str]) -> dict:
    groups, words, expected, _ = design(catalogue)
    refused = design_refusal(expected, catalogue, len(reads))
    if refused is not None:
        names = next(iter(refused["groups"]))
        return dict(**refused, iterations=0, rejected_ambiguous=0,
                    assigned=len(reads), read_assignment=[names] * len(reads))
    matrix = observations(reads, words)
    n = len(reads)
    group_names = ["+".join(sorted(group)) for group in groups]
    assignments: list[str] = []
    counts = np.zeros(len(groups), dtype=float)
    iterations = 0
    if n:
        likelihood = matrix @ np.log(expected / np.sum(expected, axis=0))
        weights = np.full(len(groups), 1 / len(groups))
        for iterations in range(1, MAX_ITERATIONS + 1):
            joint = likelihood + np.log(np.maximum(weights, 1e-12))
            joint -= np.max(joint, axis=1, keepdims=True)
            posterior = np.exp(joint)
            posterior /= np.sum(posterior, axis=1, keepdims=True)
            updated = np.maximum(np.mean(posterior, axis=0), 1e-12)
            updated /= np.sum(updated)
            if np.max(np.abs(updated - weights)) < 1e-9:
                weights = updated
                break
            weights = updated
        winner = np.argmax(posterior, axis=1)
        confidence = np.max(posterior, axis=1)
        for index, certainty in zip(winner, confidence):
            assignments.append(group_names[index] if certainty >= POSTERIOR_MIN else "ambiguous")
        counts = np.array([assignments.count(name) for name in group_names], dtype=float)
    families, unresolved = expand_estimates(groups, counts)
    return dict(family=families, groups=unresolved, iterations=iterations,
                rejected_ambiguous=assignments.count("ambiguous"),
                assigned=int(sum(counts)), read_assignment=assignments)
