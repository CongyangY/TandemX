"""Audit assembly-localized TandemX arrays against author annotations.

This is an orthogonal interpretation audit. It does not alter family discovery,
localization, source eligibility, or collapse-classification thresholds.
Coordinates are normalized to 0-based half-open intervals before overlap.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable, Sequence

from benchmarks.challenge.schema import write_table
from tandemx.compare.mvp import read_arrays_bed, union_interval_length


DEFAULT_REPEAT_CLASSES = ("centromere", "5S_rDNA", "45S_rDNA", "telomere")

CLASS_FIELDS = [
    "annotation_class",
    "annotation_bp",
    "tandemx_overlap_bp",
    "annotation_bp_recall",
    "overlapping_family_count",
]
PAIR_FIELDS = ["family_id", "annotation_class", "overlap_bp"]
FAMILY_FIELDS = [
    "family_id",
    "tandemx_array_bp",
    "selected_annotation_overlap_bp",
    "selected_annotation_fraction_of_array",
    "best_annotation_class",
    "best_annotation_overlap_bp",
    "overlap_class_count",
    "warning",
]


@dataclass(frozen=True)
class AnnotationInterval:
    chrom: str
    start: int
    end: int
    annotation_class: str


def read_gff_intervals(path: Path) -> list[AnnotationInterval]:
    """Read strict GFF3 rows and convert 1-based inclusive to BED coordinates."""
    intervals: list[AnnotationInterval] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path} line {line_number} must contain 9 GFF3 fields")
            try:
                start_1based = int(fields[3])
                end_1based = int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{path} line {line_number} start/end must be integers") from exc
            if start_1based < 1 or end_1based < start_1based:
                raise ValueError(f"{path} line {line_number} has invalid 1-based inclusive coordinates")
            if not fields[0] or not fields[2] or fields[2] == ".":
                raise ValueError(f"{path} line {line_number} requires chromosome and feature type")
            intervals.append(
                AnnotationInterval(
                    chrom=fields[0],
                    start=start_1based - 1,
                    end=end_1based,
                    annotation_class=fields[2],
                )
            )
    if not intervals:
        raise ValueError(f"No annotation intervals found: {path}")
    return intervals


def _union_by_chrom(intervals: Iterable[tuple[str, int, int]]) -> dict[str, list[tuple[int, int]]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        grouped[chrom].append((start, end))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, values in grouped.items():
        ordered = sorted(values)
        current: list[list[int]] = []
        for start, end in ordered:
            if not current or start > current[-1][1]:
                current.append([start, end])
            else:
                current[-1][1] = max(current[-1][1], end)
        merged[chrom] = [(start, end) for start, end in current]
    return merged


def _intersection_length(
    left: dict[str, list[tuple[int, int]]],
    right: dict[str, list[tuple[int, int]]],
) -> int:
    total = 0
    for chrom in left.keys() & right.keys():
        left_intervals = left[chrom]
        right_intervals = right[chrom]
        left_index = right_index = 0
        while left_index < len(left_intervals) and right_index < len(right_intervals):
            left_start, left_end = left_intervals[left_index]
            right_start, right_end = right_intervals[right_index]
            total += max(0, min(left_end, right_end) - max(left_start, right_start))
            if left_end <= right_end:
                left_index += 1
            else:
                right_index += 1
    return total


def audit_annotation_overlap(
    arrays_path: Path,
    annotation_path: Path,
    selected_classes: Sequence[str] = DEFAULT_REPEAT_CLASSES,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Return class, nonzero family-class, family, and global audit records."""
    selected = tuple(dict.fromkeys(selected_classes))
    if not selected or any(not value for value in selected):
        raise ValueError("selected_classes must contain at least one non-empty class")
    selected_set = set(selected)
    arrays = read_arrays_bed(arrays_path)
    if not arrays:
        raise ValueError(f"No TandemX arrays found: {arrays_path}")
    annotations = read_gff_intervals(annotation_path)
    available_classes = {interval.annotation_class for interval in annotations}
    missing_classes = selected_set - available_classes
    if missing_classes:
        raise ValueError(f"Selected annotation classes absent from GFF3: {sorted(missing_classes)}")

    array_by_family: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    arrays_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for array in arrays:
        array_by_family[array.family_id][array.chrom].append((array.start, array.end))
        arrays_by_chrom[array.chrom].append((array.start, array.end, array.family_id))

    annotation_by_class: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    selected_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for interval in annotations:
        if interval.annotation_class in selected_set:
            annotation_by_class[interval.annotation_class][interval.chrom].append((interval.start, interval.end))
            selected_by_chrom[interval.chrom].append((interval.start, interval.end, interval.annotation_class))

    pair_intersections: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    for chrom, chrom_arrays in arrays_by_chrom.items():
        chrom_annotations = sorted(selected_by_chrom.get(chrom, ()))
        for array_start, array_end, family_id in sorted(chrom_arrays):
            for annotation_start, annotation_end, annotation_class in chrom_annotations:
                if annotation_end <= array_start:
                    continue
                if annotation_start >= array_end:
                    break
                start = max(array_start, annotation_start)
                end = min(array_end, annotation_end)
                if end > start:
                    pair_intersections[(family_id, annotation_class, chrom)].append((start, end))

    pair_bp: dict[tuple[str, str], int] = defaultdict(int)
    for (family_id, annotation_class, _chrom), intervals in pair_intersections.items():
        pair_bp[(family_id, annotation_class)] += union_interval_length(intervals)
    pair_rows = [
        {"family_id": family_id, "annotation_class": annotation_class, "overlap_bp": overlap_bp}
        for (family_id, annotation_class), overlap_bp in sorted(pair_bp.items())
    ]

    family_rows: list[dict[str, object]] = []
    selected_union = _union_by_chrom(
        (chrom, start, end)
        for values in annotation_by_class.values()
        for chrom, intervals in values.items()
        for start, end in intervals
    )
    for family_id in sorted(array_by_family):
        family_union = _union_by_chrom(
            (chrom, start, end)
            for chrom, intervals in array_by_family[family_id].items()
            for start, end in intervals
        )
        array_bp = sum(union_interval_length(intervals) for intervals in family_union.values())
        selected_overlap_bp = _intersection_length(family_union, selected_union)
        class_overlaps = {
            annotation_class: pair_bp[(family_id, annotation_class)]
            for annotation_class in selected
            if pair_bp.get((family_id, annotation_class), 0) > 0
        }
        best_class = max(class_overlaps, key=lambda value: (class_overlaps[value], value)) if class_overlaps else None
        family_rows.append(
            {
                "family_id": family_id,
                "tandemx_array_bp": array_bp,
                "selected_annotation_overlap_bp": selected_overlap_bp,
                "selected_annotation_fraction_of_array": selected_overlap_bp / array_bp,
                "best_annotation_class": best_class,
                "best_annotation_overlap_bp": class_overlaps.get(best_class, 0) if best_class else 0,
                "overlap_class_count": len(class_overlaps),
                "warning": "" if class_overlaps else "no_selected_author_annotation_overlap",
            }
        )

    all_array_union = _union_by_chrom(
        (array.chrom, array.start, array.end) for array in arrays
    )
    class_rows: list[dict[str, object]] = []
    for annotation_class in selected:
        class_union = _union_by_chrom(
            (chrom, start, end)
            for chrom, intervals in annotation_by_class[annotation_class].items()
            for start, end in intervals
        )
        annotation_bp = sum(union_interval_length(intervals) for intervals in class_union.values())
        overlap_bp = _intersection_length(all_array_union, class_union)
        class_rows.append(
            {
                "annotation_class": annotation_class,
                "annotation_bp": annotation_bp,
                "tandemx_overlap_bp": overlap_bp,
                "annotation_bp_recall": overlap_bp / annotation_bp if annotation_bp else None,
                "overlapping_family_count": len(
                    {family_id for family_id, value in pair_bp if value == annotation_class}
                ),
            }
        )

    total_array_bp = sum(union_interval_length(intervals) for intervals in all_array_union.values())
    total_annotation_bp = sum(union_interval_length(intervals) for intervals in selected_union.values())
    total_overlap_bp = _intersection_length(all_array_union, selected_union)
    summary = {
        "schema_version": 1,
        "interpretation": "orthogonal_author_annotation_overlap_not_used_for_threshold_selection",
        "coordinate_normalization": "GFF3_1_based_inclusive_to_0_based_half_open",
        "selected_annotation_classes": list(selected),
        "tandemx_family_count": len(family_rows),
        "tandemx_array_union_bp": total_array_bp,
        "selected_annotation_union_bp": total_annotation_bp,
        "overlap_union_bp": total_overlap_bp,
        "selected_annotation_bp_recall": total_overlap_bp / total_annotation_bp if total_annotation_bp else None,
        "tandemx_array_fraction_in_selected_annotation": total_overlap_bp / total_array_bp if total_array_bp else None,
        "families_with_selected_annotation_overlap": sum(
            int(row["selected_annotation_overlap_bp"]) > 0 for row in family_rows
        ),
        "families_without_selected_annotation_overlap": sum(
            int(row["selected_annotation_overlap_bp"]) == 0 for row in family_rows
        ),
    }
    return class_rows, pair_rows, family_rows, summary


def run_annotation_audit(
    arrays_path: Path,
    annotation_path: Path,
    outdir: Path,
    selected_classes: Sequence[str] = DEFAULT_REPEAT_CLASSES,
) -> dict[str, object]:
    class_rows, pair_rows, family_rows, summary = audit_annotation_overlap(
        arrays_path, annotation_path, selected_classes
    )
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir / "annotation_class_summary.tsv", class_rows, CLASS_FIELDS)
    write_table(outdir / "family_annotation_class_overlap.tsv", pair_rows, PAIR_FIELDS)
    write_table(outdir / "family_annotation_summary.tsv", family_rows, FAMILY_FIELDS)
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary
