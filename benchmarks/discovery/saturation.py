"""Sequence-set and decision metrics for discovery saturation curves."""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
import math

from benchmarks.challenge.evaluate import canonical_monomer, maximum_matching
from benchmarks.challenge.sequence_metrics import cyclic_reaches_threshold


def unique_canonical(sequences: Iterable[str]) -> list[str]:
    """Return sorted unique monomers under rotation/strand canonicalization."""
    return sorted({canonical_monomer(sequence) for sequence in sequences if sequence})


def catalogue_matching(
    lower: Sequence[str], higher: Sequence[str], threshold: float = 0.90
) -> dict[int, int]:
    """One-to-one cyclic sequence matching from lower to higher catalogue."""
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    left = unique_canonical(lower)
    right = unique_canonical(higher)
    edges = [
        [
            j
            for j, sequence in enumerate(right)
            if cyclic_reaches_threshold(query, sequence, threshold)
        ]
        for query in left
    ]
    return maximum_matching(edges)


def transition_metrics(
    lower: Sequence[str],
    higher: Sequence[str],
    lower_depth: float,
    higher_depth: float,
    threshold: float = 0.90,
) -> dict[str, float | int | bool]:
    """Calculate adjacent-depth catalogue stability and discovery yield."""
    if not (math.isfinite(lower_depth) and math.isfinite(higher_depth)) or not (
        0 < lower_depth < higher_depth
    ):
        raise ValueError("depths must be finite, positive and increasing")
    left = unique_canonical(lower)
    right = unique_canonical(higher)
    matched = len(catalogue_matching(left, right, threshold))
    new = len(right) - matched
    lost = len(left) - matched
    union = len(left) + len(right) - matched
    jaccard = matched / union if union else 1.0
    new_fraction = new / len(right) if right else 0.0
    increment = higher_depth - lower_depth
    return {
        "lower_family_count": len(left),
        "higher_family_count": len(right),
        "matched_family_count": matched,
        "new_family_count": new,
        "lost_family_count": lost,
        "new_family_fraction": new_fraction,
        "new_families_per_added_x": new / increment,
        "family_jaccard": jaccard,
        "transition_pass": new_fraction < 0.05 and jaccard > 0.95,
    }


def saturation_decision(
    rows: Sequence[dict], seeds: Sequence[int], coverages: Sequence[float]
) -> dict:
    """Apply the all-seed, two-consecutive-transition saturation rule."""
    if len(set(seeds)) != len(seeds) or len(coverages) < 3:
        raise ValueError("require unique seeds and at least three coverages")
    ordered = sorted(float(value) for value in coverages)
    if ordered != list(coverages) or len(set(ordered)) != len(ordered):
        raise ValueError("coverages must be unique and increasing")
    by_higher: dict[float, list[dict]] = defaultdict(list)
    for row in rows:
        by_higher[float(row["higher_depth"])].append(row)
    transition_summary = []
    for index in range(1, len(ordered)):
        low, high = ordered[index - 1], ordered[index]
        group = by_higher.get(high, [])
        valid = {
            int(row["seed"]): row
            for row in group
            if float(row["lower_depth"]) == low
        }
        all_seed_pass = set(valid) == set(seeds) and all(
            bool(valid[seed]["transition_pass"]) for seed in seeds
        )
        transition_summary.append(
            {
                "lower_depth": low,
                "higher_depth": high,
                "all_seed_pass": all_seed_pass,
                "passing_seeds": sum(
                    bool(valid.get(seed, {}).get("transition_pass")) for seed in seeds
                ),
            }
        )
    saturation_depth = None
    for previous, current in zip(transition_summary, transition_summary[1:]):
        if previous["all_seed_pass"] and current["all_seed_pass"]:
            saturation_depth = current["higher_depth"]
            break
    return {
        "saturation_reached": saturation_depth is not None,
        "saturation_depth": saturation_depth,
        "rule": "two_consecutive_transitions_all_seeds_new_fraction_lt_0.05_and_jaccard_gt_0.95",
        "transitions": transition_summary,
    }
