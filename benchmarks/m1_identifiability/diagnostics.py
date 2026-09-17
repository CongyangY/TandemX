"""Input-only structural identifiability; no abundance estimator or fitted threshold."""

from __future__ import annotations

from collections import Counter
from math import log2
from typing import Mapping

import numpy as np

from benchmarks.m1_shared_signature.model import circular_kmers, rotations


def catalogue_diagnostics(catalogue: Mapping[str, str], k: int) -> dict:
    if not catalogue:
        raise ValueError("A nonempty catalogue is required")
    if any(not unit or set(unit) - set("ACGT") for unit in catalogue.values()):
        raise ValueError("Catalogue units must contain only nonempty A/C/G/T sequences")
    names = sorted(catalogue)
    if len(set(names)) != len(names):
        raise ValueError("Duplicate family identifiers")
    signatures = {name: circular_kmers(catalogue[name], k) for name in names}
    vocabulary = sorted(set().union(*(set(signature) for signature in signatures.values())))
    matrix = np.array(
        [[signatures[name][word] for name in names] for word in vocabulary], dtype=float
    )
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    _, _, vt = np.linalg.svd(matrix, full_matrices=True)
    rank = int(np.linalg.matrix_rank(matrix))
    nullspace = vt[rank:, :] if rank < len(names) else np.empty((0, len(names)))
    # The nullspace indicates exactly which coefficients cannot be recovered
    # from noiseless aggregate catalogue k-mer counts. It is not a noise model.
    null_involvement = np.sum(nullspace**2, axis=0)
    condition = (
        float(singular_values[0] / singular_values[-1])
        if rank == len(names) else None
    )
    family: dict[str, dict] = {}
    for index, name in enumerate(names):
        signature = signatures[name]
        unique = [word for word in signature if all(
            signatures[other][word] == 0 for other in names if other != name
        )]
        unique_bp_fraction = sum(signature[word] for word in unique) / len(catalogue[name])
        if null_involvement[index] > 1e-8:
            status = "NON_IDENTIFIABLE"
            reason = "aggregate_kmer_matrix_nullspace_involves_family"
        elif not unique:
            status = "PARTIALLY_IDENTIFIABLE"
            reason = "full_rank_but_no_exclusive_catalogue_kmer"
        else:
            status = "IDENTIFIABLE"
            reason = "noiseless_catalogue_signature_separable"
        family[name] = {
            "status": status,
            "reason": reason,
            "exclusive_kmer_count": len(unique),
            "exclusive_kmer_position_fraction": unique_bp_fraction,
            "nullspace_involvement": float(null_involvement[index]),
        }
    pairwise = []
    for left_index, left in enumerate(names):
        for right in names[left_index + 1:]:
            a, b = signatures[left], signatures[right]
            words = set(a) | set(b)
            shared = sum(min(a[word], b[word]) for word in words)
            union = sum(max(a[word], b[word]) for word in words)
            pairwise.append({
                "left": left, "right": right,
                "weighted_shared_kmer_fraction": shared / union if union else 1.0,
                "identical_kmer_signature": a == b,
            })
    return {
        "k": k,
        "family_order": names,
        "n_kmers": len(vocabulary),
        "signature_rank": rank,
        "family_count": len(names),
        "condition_number_full_rank_only": condition,
        "singular_values": [float(value) for value in singular_values],
        "effective_identifiable_family_dimensions": rank,
        "family": family,
        "pairwise": pairwise,
        "status_scope": "noiseless_given_catalogue_kmer_features_only",
    }


def read_assignment_diagnostics(
    catalogue: Mapping[str, str], reads: list[str], source_labels: list[str],
    max_mismatch: int,
) -> dict:
    """Observed cyclic-Hamming attribution on known-source development reads."""
    if len(reads) != len(source_labels):
        raise ValueError("Read and source-label counts differ")
    names = sorted(catalogue)
    lengths = {len(unit) for unit in catalogue.values()}
    if len(lengths) != 1:
        raise ValueError("Equal-length units are required for this diagnostic")
    references = {name: rotations(catalogue[name]) for name in names}
    dna = "ACGT"
    counts = {name: Counter() for name in sorted(set(source_labels))}
    tie_entropy = {name: 0.0 for name in counts}
    for read, source in zip(reads, source_labels):
        if len(read) not in lengths or set(read) - set(dna):
            raise ValueError("Reads must be equal-length A/C/G/T strings")
        query = np.fromiter((dna.index(base) for base in read), dtype=np.uint8)
        distances = {name: int(np.min(np.count_nonzero(ref != query, axis=1)))
                     for name, ref in references.items()}
        best = min(distances.values())
        counts[source]["total"] += 1
        if best > max_mismatch:
            counts[source]["rejected"] += 1
            continue
        winners = [name for name, distance in distances.items() if distance == best]
        counts[source]["accepted"] += 1
        tie_entropy[source] += log2(len(winners))
        if len(winners) > 1:
            counts[source]["tied_best"] += 1
        elif winners[0] == source:
            counts[source]["correct_unique"] += 1
        else:
            counts[source]["wrong_unique"] += 1
    return {
        name: {
            **{key: int(counts[name][key]) for key in (
                "total", "accepted", "rejected", "tied_best",
                "correct_unique", "wrong_unique",
            )},
            "uniquely_correct_fraction_all_source_reads": (
                counts[name]["correct_unique"] / counts[name]["total"]
                if counts[name]["total"] else None
            ),
            "mean_best_hit_tie_entropy_bits_per_accepted": (
                tie_entropy[name] / counts[name]["accepted"]
                if counts[name]["accepted"] else None
            ),
        }
        for name in counts
    }
