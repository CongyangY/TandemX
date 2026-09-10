"""Toy-only, mask-based metrics for the phase-gate fixture.

This scorer consumes predictions separately from fixture truth. It does not call
a model, select thresholds, or write truth into observable inputs.
"""

from __future__ import annotations

from collections import defaultdict

from .context_prototype import merge_intervals


def score_predictions(
    predictions: dict[str, dict[str, list[tuple[int, int]]]],
    truth_intervals: list[tuple[str, int, int, str]],
    read_lengths: dict[str, int],
    supplied_families: set[str] | list[str] | tuple[str, ...],
) -> dict[str, object]:
    """Score predictions against generated read-coordinate truth masks.

    ``supplied_families`` comes from the independent fixture catalogue, so a
    catalogue family with zero sampled truth remains a valid output row.
    Per-base masks are appropriate only for the small development fixture. A
    predicted base covered by more than one family is always ambiguous and is
    excluded from family TP/FP counts. ``unassigned_bp`` is the number of read
    bases with no prediction; unique, ambiguous and unassigned bases conserve
    total read length.
    """
    known_families = set(supplied_families)
    if len(known_families) != len(list(supplied_families)):
        raise ValueError("Catalogue family names must be unique")
    if any(length < 0 for length in read_lengths.values()):
        raise ValueError("Read lengths must be non-negative")
    truth: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for read_id, start, end, family in truth_intervals:
        if read_id not in read_lengths:
            raise ValueError("Unknown truth read ID")
        if family not in known_families or start < 0 or end <= start or end > read_lengths[read_id]:
            raise ValueError("Invalid truth interval")
        truth[read_id][family].append((start, end))
    for read_id in truth:
        for family in truth[read_id]:
            truth[read_id][family] = merge_intervals(truth[read_id][family])

    for read_id, families in predictions.items():
        if read_id not in read_lengths:
            raise ValueError("Unknown prediction read ID")
        for family, rows in families.items():
            if family not in known_families:
                raise ValueError("Unknown prediction family")
            for start, end in rows:
                if start < 0 or end <= start or end > read_lengths[read_id]:
                    raise ValueError("Prediction interval exceeds read length")

    families = sorted(known_families)
    per = {
        family: {
            "read_truth_bp": 0,
            "predicted_unique_bp": 0,
            "ambiguous_eligible_bp": 0,
            "true_positive_bp": 0,
            "false_positive_bp": 0,
            "background_false_bp": 0,
            "recall": None,
            "precision": None,
            "absolute_relative_error": None,
        }
        for family in families
    }
    ambiguous_total = 0
    unassigned_total = 0
    unique_total = 0
    for read_id, length in read_lengths.items():
        pred_rows: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for family, rows in predictions.get(read_id, {}).items():
            pred_rows[family].extend(rows)
        for family in pred_rows:
            pred_rows[family] = merge_intervals(pred_rows[family])
        events: dict[int, list[tuple[str, int]]] = defaultdict(list)
        for family, rows in pred_rows.items():
            for start, end in rows:
                events[start].append((family, 1))
                events[end].append((family, -1))
        active: set[str] = set()
        last = 0
        truth_points = {
            point
            for rows in truth.get(read_id, {}).values()
            for start, end in rows
            for point in (start, end)
        }
        for pos in sorted(set(events) | truth_points | {length}):
            width = pos - last
            truth_at = {
                family for family, rows in truth.get(read_id, {}).items()
                if any(start <= last and pos <= end for start, end in rows)
            }
            if len(active) > 1:
                ambiguous_total += width
                for family in active:
                    per[family]["ambiguous_eligible_bp"] += width
            elif len(active) == 1:
                family = next(iter(active))
                per[family]["predicted_unique_bp"] += width
                unique_total += width
                if truth_at == {family}:
                    per[family]["true_positive_bp"] += width
                else:
                    per[family]["false_positive_bp"] += width
                    if not truth_at:
                        per[family]["background_false_bp"] += width
            else:
                unassigned_total += width
            for family, delta in events.get(pos, []):
                if delta == 1:
                    active.add(family)
                else:
                    active.remove(family)
            last = pos

    for read_id, family_rows in truth.items():
        for family, rows in family_rows.items():
            per[family]["read_truth_bp"] += sum(end - start for start, end in rows)
    for row in per.values():
        truth_bp = row["read_truth_bp"]
        predicted_bp = row["predicted_unique_bp"]
        tp = row["true_positive_bp"]
        row["recall"] = tp / truth_bp if truth_bp else None
        row["precision"] = tp / predicted_bp if predicted_bp else None
        row["absolute_relative_error"] = (
            abs(predicted_bp - truth_bp) / truth_bp if truth_bp else None
        )
    total_read_bases = sum(read_lengths.values())
    if unique_total + ambiguous_total + unassigned_total != total_read_bases:
        raise AssertionError("Prediction mask mass is not conserved")
    return {
        "per_family": per,
        "overall": {
            "ambiguous_bp": ambiguous_total,
            "unassigned_bp": unassigned_total,
            "predicted_unique_bp": unique_total,
            "read_bases": total_read_bases,
            "mass_conserved": True,
        },
    }
