"""A deliberately simple, two-flank array-length discordance baseline.

This method cannot classify order-only edits or reads lacking both unique flanks.
It is not an assembly-error oracle: it compares observed read and assembly spans.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Mapping
from typing import Any

import edlib


def _unique_flank_hit(flank: str, sequence: str) -> tuple[int, int] | None:
    if not flank or not sequence:
        return None
    tolerance = max(2, math.ceil(len(flank) * 0.10))
    hit = edlib.align(flank, sequence, mode="HW", task="locations", k=tolerance)
    locations = hit.get("locations", [])
    if hit["editDistance"] < 0 or len(locations) != 1:
        return None
    start, end_inclusive = locations[0]
    return start, end_inclusive + 1


def _array_span(sequence: str, left: str, right: str) -> tuple[int, int] | None:
    left_hit = _unique_flank_hit(left, sequence)
    right_hit = _unique_flank_hit(right, sequence)
    if left_hit is None or right_hit is None or left_hit[1] > right_hit[0]:
        return None
    return left_hit[1], right_hit[0]


def predict_case(row: Mapping[str, Any]) -> dict[str, Any]:
    """Return schema-v2 event decision from paired-flank span evidence only."""
    case_id = str(row["case_id"])
    base: dict[str, Any] = {
        "case_id": case_id,
        "route": "two_flank_array_length_v1",
        "status": "insufficient_read_support",
        "audit_state": "INSUFFICIENT_READ_SUPPORT",
        "event_score": None,
        "event_type": None,
        "predicted_edited_interval_bp": None,
        "predicted_edited_label_path": None,
        "predicted_edited_orientation_path": None,
        "predicted_copy_spans_bp": None,
        "predicted_signed_bp_delta": None,
        "reason": "missing_unique_two_flank_span",
    }
    left = str(row["left_flank_sequence"])
    right = str(row["right_flank_sequence"])
    assembly = _array_span(str(row["assembly_sequence"]), left, right)
    if assembly is None:
        return base
    reads = row["raw_read_sequences"]
    if not isinstance(reads, dict):
        return base
    spans = [_array_span(str(sequence), left, right) for sequence in reads.values()]
    read_lengths = [end - start for span in spans if span is not None for start, end in [span]]
    if len(read_lengths) < 3:
        return base
    assembly_length = assembly[1] - assembly[0]
    median_read_length = statistics.median(read_lengths)
    delta = median_read_length - assembly_length
    threshold = max(2, math.ceil(max(assembly_length, median_read_length) * 0.05))
    if abs(delta) > threshold:
        base.update(
            status="ok",
            audit_state="DISCORDANT",
            event_score=1.0,
            event_type="array_length_discordance",
            predicted_edited_interval_bp=list(assembly),
            predicted_signed_bp_delta=-delta,
            reason="median_two_flank_read_span_differs_from_assembly",
        )
    else:
        base.update(
            status="ok",
            audit_state="SUPPORTED",
            event_score=0.0,
            predicted_edited_interval_bp=list(assembly),
            predicted_signed_bp_delta=-delta,
            reason="median_two_flank_read_span_within_fixed_tolerance",
        )
    return base
