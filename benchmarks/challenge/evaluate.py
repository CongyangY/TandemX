"""Independent truth scoring; no TandemX detector or matcher is used here."""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from functools import lru_cache

from .schema import ArrayRecord


def interval_iou(a: ArrayRecord, b: ArrayRecord) -> float:
    if a.read_id != b.read_id:
        return 0.0
    intersection = max(0, min(a.end, b.end) - max(a.start, b.start))
    return intersection / (a.end - a.start + b.end - b.start - intersection)


def maximum_matching(edges: list[list[int]]) -> dict[int, int]:
    """Deterministic maximum-cardinality bipartite matching, left -> right."""
    owners: dict[int, int] = {}

    def augment(left: int, visited: set[int]) -> bool:
        for right in edges[left]:
            if right in visited:
                continue
            visited.add(right)
            if right not in owners or augment(owners[right], visited):
                owners[right] = left
                return True
        return False

    for left in range(len(edges)):
        augment(left, set())
    return {left: right for right, left in owners.items()}


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if not 0 <= successes <= total:
        raise ValueError("Invalid binomial counts")
    if total == 0:
        return math.nan, math.nan
    p = successes / total
    divisor = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / divisor
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / divisor
    return (0.0 if successes == 0 else max(0.0, centre - half),
            1.0 if successes == total else min(1.0, centre + half))


def canonical_monomer(sequence: str) -> str:
    """Normalize rotations and reverse complements independently of TandemX."""
    if not sequence:
        return ""
    sequence = sequence.upper()
    reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(strand[i:] + strand[:i] for strand in (sequence, reverse) for i in range(len(sequence)))


@lru_cache(maxsize=4096)
def circular_identity(a: str, b: str) -> float:
    """Exact maximum ungapped identity across both strands and all rotations.

    Only equal-length monomers are compared. This stricter sequence endpoint is
    reported separately from tolerance-based period/array recall.
    """
    if not a or len(a) != len(b):
        return 0.0
    a, b = a.upper(), b.upper()
    reverse = b.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    if a in b + b or a in reverse + reverse:
        return 1.0
    best = 0
    for strand in (b, reverse):
        doubled = strand + strand
        for offset in range(len(a)):
            best = max(best, sum(x == y and x in "ACGT" for x, y in zip(a, doubled[offset:offset + len(a)])))
    return best / len(a)


def score_arrays(predicted: list[ArrayRecord], truth: list[ArrayRecord],
                 read_lengths: dict[str, int], min_iou: float = 0.5) -> tuple[dict, list[dict]]:
    if not 0 < min_iou <= 1 or not read_lengths:
        raise ValueError("Need reads and IoU in (0,1]")
    for record in predicted + truth:
        if record.read_id not in read_lengths or record.end > read_lengths[record.read_id]:
            raise ValueError(f"Unknown read or out-of-bounds interval: {record}")
    by_read: dict[str, list[int]] = defaultdict(list)
    for j, record in enumerate(truth):
        by_read[record.read_id].append(j)
    edges = [[j for j in by_read[p.read_id]
              if interval_iou(p, truth[j]) >= min_iou
              and abs(p.period - truth[j].period) <= max(2, round(0.02 * truth[j].period))]
             for p in predicted]
    matched = maximum_matching(edges)
    positive_reads = {r.read_id for r in truth}
    predicted_reads = {r.read_id for r in predicted}
    negative_reads = set(read_lengths) - positive_reads
    read_tp = len(predicted_reads & positive_reads)
    read_fp = len(predicted_reads & negative_reads)
    tp = len(matched)
    metrics = {
        "truth_array_count": len(truth), "predicted_array_count": len(predicted), "matched_array_count": tp,
        "array_recall": ratio(tp, len(truth)), "array_precision": ratio(tp, len(predicted)),
        "array_f1": ratio(2 * tp, len(predicted) + len(truth)),
        "read_detection_recall": ratio(read_tp, len(positive_reads)),
        "read_detection_precision": ratio(read_tp, len(predicted_reads)),
        "negative_read_count": len(negative_reads), "false_positive_read_count": read_fp,
        "negative_read_call_rate": ratio(read_fp, len(negative_reads)),
        "matched_period_mae_bp": statistics.fmean(abs(predicted[i].period - truth[j].period)
                                                 for i, j in matched.items()) if matched else math.nan,
        "matched_boundary_mae_bp": statistics.fmean((abs(predicted[i].start - truth[j].start)
                                                      + abs(predicted[i].end - truth[j].end)) / 2
                                                     for i, j in matched.items()) if matched else math.nan,
        "warning": "planted_truth;conditional_errors_on_matches;wilson_intervals_are_descriptive",
    }
    for name, numerator, denominator in (("array_recall", tp, len(truth)),
                                         ("array_precision", tp, len(predicted)),
                                         ("negative_read_call_rate", read_fp, len(negative_reads))):
        metrics[name + "_ci_low"], metrics[name + "_ci_high"] = wilson_interval(numerator, denominator)
    details = []
    for i, p in enumerate(predicted):
        target = truth[matched[i]] if i in matched else None
        details.append({"prediction_index": i, "read_id": p.read_id, "start": p.start, "end": p.end,
                        "period": p.period, "matched_truth_index": matched.get(i, "NA"),
                        "truth_family_id": target.family_id if target else "NA",
                        "iou": interval_iou(p, target) if target else 0.0,
                        "status": "matched" if target else "unmatched"})
    return metrics, details


def score_families(predicted_sequences: list[str], truth_families: dict[str, str],
                   identity_threshold: float = 0.9) -> tuple[dict, list[dict]]:
    """Sequence-supported recovery with one-to-one truth/prediction matching.

    Exact duplicate sequences are deduplicated; predictions are not clustered
    into families by length. This endpoint is recovery, not family precision.
    """
    if not 0 <= identity_threshold <= 1:
        raise ValueError("Identity threshold must be in [0,1]")
    sequences = sorted(set(canonical_monomer(s) for s in predicted_sequences if s))
    families = sorted(truth_families)
    scores = [[circular_identity(truth_families[f], seq) for seq in sequences] for f in families]
    edges = [[j for j, value in enumerate(row) if value >= identity_threshold] for row in scores]
    matches = maximum_matching(edges)
    details = [{"family_id": family, "best_identity": max(scores[i], default=0.0),
                "recovered": int(i in matches),
                "assigned_sequence_index": matches.get(i, "NA"),
                "criterion": f"equal_length_circular_ungapped_identity_ge_{identity_threshold:g}"}
               for i, family in enumerate(families)]
    return {"truth_family_count": len(families), "recovered_family_count": len(matches),
            "sequence_family_recall": ratio(len(matches), len(families)),
            "family_criterion": f"equal_length_circular_ungapped_identity_ge_{identity_threshold:g}"}, details
