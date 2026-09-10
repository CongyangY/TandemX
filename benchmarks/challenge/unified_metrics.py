"""Independent unified base-union metrics for native comparator outputs.

``interval_metrics`` is a toy/development scorer. Predictions and truth remain
separate inputs. Native family IDs are mapped to truth family IDs only through
the supplied correspondence table; unmatched or ambiguous native outputs stay
in the global union denominator and are reported as unassigned native bases.
All intervals are 0-based, half-open, and checked against ``read_lengths``.
"""

from __future__ import annotations

from collections import defaultdict

from .schema import ArrayRecord


def _union(rows: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(rows):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def interval_metrics(
    predictions: list[ArrayRecord],
    truths: list[ArrayRecord],
    read_lengths: dict[str, int],
    native_to_truth: dict[str, str | None],
) -> dict[str, object]:
    """Return global and mapped-family base-union endpoints.

    Same-family duplicate predictions are unioned. Bases covered by two or more
    mapped truth families are ambiguous and excluded from family TP/FP rather
    than forced to one family. ``negative_read_predicted_bp`` counts the global
    prediction union on reads with no truth array. ``unassigned_native_bp`` is
    the union of predictions whose native family has no unique mapping.
    """
    if any(length < 0 for length in read_lengths.values()):
        raise ValueError("Read lengths must be non-negative")
    for record in predictions + truths:
        if record.read_id not in read_lengths:
            raise ValueError(f"Unknown read ID: {record.read_id}")
        if record.end > read_lengths[record.read_id]:
            raise ValueError(f"Interval exceeds read length: {record}")
    truth_by_read: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for record in truths:
        truth_by_read[record.read_id][record.family_id].append((record.start, record.end))
    for by_family in truth_by_read.values():
        for family in by_family:
            by_family[family] = _union(by_family[family])

    pred_by_read: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    mapped_by_read: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    unmatched_by_read: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for record in predictions:
        pred_by_read[record.read_id][record.family_id].append((record.start, record.end))
        target = native_to_truth.get(record.family_id)
        if target:
            mapped_by_read[record.read_id][target].append((record.start, record.end))
        else:
            unmatched_by_read[record.read_id][record.family_id].append((record.start, record.end))
    for by_read in (pred_by_read, mapped_by_read, unmatched_by_read):
        for by_family in by_read.values():
            for family in by_family:
                by_family[family] = _union(by_family[family])

    all_truth_families = {record.family_id for record in truths}
    mapped_truth_families = {target for target in native_to_truth.values() if target}
    families = sorted(all_truth_families | mapped_truth_families)
    per = {
        family: {"truth_bp": 0, "predicted_bp": 0, "true_positive_bp": 0,
                 "false_positive_bp": 0, "ambiguous_bp": 0}
        for family in families
    }
    truth_union_by_read: dict[str, list[tuple[int, int]]] = {}
    prediction_union_by_read: dict[str, list[tuple[int, int]]] = {}
    for read_id, length in read_lengths.items():
        truth_union_by_read[read_id] = _union(
            [interval for rows in truth_by_read.get(read_id, {}).values() for interval in rows]
        )
        prediction_union_by_read[read_id] = _union(
            [interval for rows in pred_by_read.get(read_id, {}).values() for interval in rows]
        )
        points = {0, length}
        for rows in (truth_by_read.get(read_id, {}).values(), pred_by_read.get(read_id, {}).values()):
            points.update(point for family_rows in rows for start, end in family_rows for point in (start, end))
        for pos, nxt in zip(sorted(points), sorted(points)[1:]):
            if nxt <= pos:
                continue
            truth_active = {family for family, rows in truth_by_read.get(read_id, {}).items()
                            if any(start <= pos and nxt <= end for start, end in rows)}
            pred_active = {family for family, rows in mapped_by_read.get(read_id, {}).items()
                           if any(start <= pos and nxt <= end for start, end in rows)}
            if len(pred_active) > 1:
                for family in pred_active:
                    per[family]["ambiguous_bp"] += nxt - pos
            elif len(pred_active) == 1:
                family = next(iter(pred_active))
                per[family]["predicted_bp"] += nxt - pos
                if truth_active == {family}:
                    per[family]["true_positive_bp"] += nxt - pos
                else:
                    per[family]["false_positive_bp"] += nxt - pos
        for family, rows in truth_by_read.get(read_id, {}).items():
            per[family]["truth_bp"] += sum(end - start for start, end in rows)

    global_truth_bp = sum(end - start for rows in truth_union_by_read.values() for start, end in rows)
    global_predicted_bp = sum(end - start for rows in prediction_union_by_read.values() for start, end in rows)
    intersection_bp = 0
    for read_id in read_lengths:
        for start, end in prediction_union_by_read[read_id]:
            intersection_bp += sum(max(0, min(end, right) - max(start, left))
                                   for left, right in truth_union_by_read[read_id])
    negative_reads = set(read_lengths) - set(truth_by_read)
    negative_predicted_bp = sum(sum(end - start for start, end in prediction_union_by_read[read_id])
                                for read_id in negative_reads)
    # Native families without a unique correspondence remain in the global
    # prediction union.  Their diagnostic count is a read-level union so two
    # unmatched native reports covering the same bases are not double counted.
    unassigned_native_bp = sum(
        sum(end - start for start, end in _union(
            [interval for rows in by_family.values() for interval in rows]
        ))
        for by_family in unmatched_by_read.values()
    )
    return {
        "global": {
            "truth_union_bp": global_truth_bp,
            "predicted_union_bp": global_predicted_bp,
            "true_positive_bp": intersection_bp,
            "base_union_recall": intersection_bp / global_truth_bp if global_truth_bp else None,
            "base_union_precision": intersection_bp / global_predicted_bp if global_predicted_bp else None,
            "negative_read_predicted_bp": negative_predicted_bp,
            "unassigned_native_bp": unassigned_native_bp,
        },
        "per_family": per,
    }
