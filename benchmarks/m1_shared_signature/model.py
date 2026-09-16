"""Bounded competitive mapping and shared-k-mer constrained regression.

Both methods receive the same complete, predeclared catalogue and reads. A
common alignment gate rejects reads that do not match any catalogue unit.
The regression fits aggregate circular k-mer counts, with nonnegative family
weights constrained to sum to the number of accepted reads. Equal signature
columns are reported as an unresolved group, never divided by a prior.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations

import numpy as np

DNA = "ACGT"
MIN_SINGULAR_RATIO = 0.05  # design-only safeguard; development choice, not calibrated


def rotations(sequence: str) -> np.ndarray:
    if not sequence or set(sequence) - set(DNA):
        raise ValueError("Units must be nonempty A/C/G/T sequences")
    codes = np.fromiter((DNA.index(base) for base in sequence), dtype=np.uint8)
    return np.stack([np.roll(codes, offset) for offset in range(len(codes))])


def circular_kmers(sequence: str, k: int) -> Counter[str]:
    if k < 1 or k > len(sequence):
        raise ValueError("k must be between 1 and the unit length")
    extended = sequence + sequence[: k - 1]
    return Counter(extended[i : i + k] for i in range(len(sequence)))


def align_and_gate(reads: list[str], catalogue: dict[str, str], max_mismatch: int
                   ) -> tuple[list[str], dict[str, float], int, int]:
    """Ordinary best cyclic-Hamming mapping; tied best hits are unassigned."""
    if not catalogue or max_mismatch < 0:
        raise ValueError("Need nonempty catalogue and nonnegative mismatch gate")
    lengths = {len(s) for s in catalogue.values()}
    if len(lengths) != 1:
        raise ValueError("This bounded prototype requires equal-length units")
    references = {name: rotations(seq) for name, seq in catalogue.items()}
    mapped = {name: 0.0 for name in catalogue}
    accepted: list[str] = []
    rejected = tied = 0
    for read in reads:
        if len(read) not in lengths or set(read) - set(DNA):
            raise ValueError("Read is not an equal-length A/C/G/T sequence")
        query = np.fromiter((DNA.index(base) for base in read), dtype=np.uint8)
        distances = {name: int(np.min(np.count_nonzero(ref != query, axis=1)))
                     for name, ref in references.items()}
        best = min(distances.values())
        if best > max_mismatch:
            rejected += 1
            continue
        accepted.append(read)
        winners = [name for name, distance in distances.items() if distance == best]
        if len(winners) == 1:
            mapped[winners[0]] += 1.0
        else:
            tied += 1
    return accepted, mapped, rejected, tied


def constrained_fit(columns: np.ndarray, observed: np.ndarray, total: int) -> np.ndarray:
    """Exact active-set enumeration for at most six catalogue groups."""
    count = columns.shape[1]
    if count > 6 or observed.shape != (columns.shape[0],) or total < 0:
        raise ValueError("Invalid bounded constrained regression input")
    if total == 0:
        return np.zeros(count)
    best_loss = float("inf")
    best = np.zeros(count)
    for size in range(1, count + 1):
        for indices in combinations(range(count), size):
            sub = columns[:, indices]
            # KKT equations for least squares with sum(weights) == total.
            gram = sub.T @ sub
            rhs = sub.T @ observed
            system = np.block([[gram, np.ones((size, 1))],
                               [np.ones((1, size)), np.zeros((1, 1))]])
            solution = np.linalg.lstsq(system, np.r_[rhs, total], rcond=None)[0][:size]
            if np.min(solution) < -1e-7:
                continue
            candidate = np.zeros(count)
            candidate[list(indices)] = np.maximum(solution, 0)
            loss = float(np.sum((columns @ candidate - observed) ** 2))
            if loss < best_loss:
                best_loss, best = loss, candidate
    return best


def shared_signature_fit(reads: list[str], catalogue: dict[str, str], k: int
                         ) -> tuple[dict[str, float | None], dict[str, float], list[list[str]]]:
    """Return identifiable family counts, unresolved group totals and groups."""
    if not catalogue:
        raise ValueError("Need nonempty catalogue")
    signatures = {name: circular_kmers(seq, k) for name, seq in catalogue.items()}
    groups: dict[tuple[tuple[str, int], ...], list[str]] = {}
    for name, signature in signatures.items():
        groups.setdefault(tuple(sorted(signature.items())), []).append(name)
    grouped = list(groups.values())
    vocabulary = sorted(set().union(*(signature.keys() for signature in signatures.values())))
    columns = np.array([[signatures[names[0]][word] for names in grouped]
                        for word in vocabulary], dtype=float)
    singular = np.linalg.svd(columns, compute_uv=False)
    if singular[-1] / singular[0] < MIN_SINGULAR_RATIO:
        names = sorted(catalogue)
        return ({name: None for name in names},
                {"+".join(names): float(len(reads))}, [names])
    observed_counts = Counter()
    for read in reads:
        observed_counts.update(circular_kmers(read, k))
    observed = np.array([observed_counts[word] for word in vocabulary], dtype=float)
    fitted = constrained_fit(columns, observed, len(reads))
    family: dict[str, float | None] = {}
    unresolved: dict[str, float] = {}
    for names, abundance in zip(grouped, fitted):
        if len(names) == 1:
            family[names[0]] = float(abundance)
        else:
            for name in names:
                family[name] = None
            unresolved["+".join(sorted(names))] = float(abundance)
    return family, unresolved, grouped
