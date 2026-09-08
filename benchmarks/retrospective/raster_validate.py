"""Independent byte-raster validation for the interval-sweep annotation audit."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path
from typing import Sequence


def validate_rasterized_overlap(
    arrays_path: Path,
    annotation_path: Path,
    audit_dir: Path,
    selected_classes: Sequence[str],
) -> dict[str, object]:
    """Recompute union and overlap bp without using the sweep implementation."""
    selected = tuple(dict.fromkeys(selected_classes))
    if not selected:
        raise ValueError("At least one annotation class is required")
    arrays: dict[str, list[tuple[int, int]]] = defaultdict(list)
    annotations: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    max_end: dict[str, int] = defaultdict(int)
    with arrays_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed BED line {line_number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED coordinates at line {line_number}")
            arrays[fields[0]].append((start, end))
            max_end[fields[0]] = max(max_end[fields[0]], end)
    with annotation_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF3 line {line_number}")
            if fields[2] not in selected:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid GFF3 coordinates at line {line_number}")
            annotations[fields[0]].append((start, end, fields[2]))
            max_end[fields[0]] = max(max_end[fields[0]], end)

    class_annotation_bp = {value: 0 for value in selected}
    class_overlap_bp = {value: 0 for value in selected}
    array_bp = annotation_bp = overlap_bp = 0
    for chrom, length in max_end.items():
        array_mask = bytearray(length)
        annotation_mask = bytearray(length)
        for start, end in arrays.get(chrom, ()):
            array_mask[start:end] = b"\1" * (end - start)
        for start, end, _annotation_class in annotations.get(chrom, ()):
            annotation_mask[start:end] = b"\1" * (end - start)
        array_bp += array_mask.count(1)
        annotation_bp += annotation_mask.count(1)
        overlap_bp += sum(left & right for left, right in zip(array_mask, annotation_mask))
        for annotation_class in selected:
            class_mask = bytearray(length)
            for start, end, value in annotations.get(chrom, ()):
                if value == annotation_class:
                    class_mask[start:end] = b"\1" * (end - start)
            class_annotation_bp[annotation_class] += class_mask.count(1)
            class_overlap_bp[annotation_class] += sum(
                left & right for left, right in zip(array_mask, class_mask)
            )

    expected_summary = json.loads((audit_dir / "summary.json").read_text(encoding="utf-8"))
    with (audit_dir / "annotation_class_summary.tsv").open(encoding="utf-8", newline="") as handle:
        expected_classes = {row["annotation_class"]: row for row in csv.DictReader(handle, delimiter="\t")}
    checks = {
        "tandemx_array_union_bp": array_bp == int(expected_summary["tandemx_array_union_bp"]),
        "selected_annotation_union_bp": annotation_bp == int(expected_summary["selected_annotation_union_bp"]),
        "overlap_union_bp": overlap_bp == int(expected_summary["overlap_union_bp"]),
        "class_annotation_bp": all(
            value in expected_classes
            and class_annotation_bp[value] == int(expected_classes[value]["annotation_bp"])
            for value in selected
        ),
        "class_overlap_bp": all(
            value in expected_classes
            and class_overlap_bp[value] == int(expected_classes[value]["tandemx_overlap_bp"])
            for value in selected
        ),
    }
    receipt = {
        "schema_version": 1,
        "method": "independent_per_chromosome_byte_raster",
        "coordinate_normalization": "GFF3_1_based_inclusive_to_0_based_half_open",
        "selected_annotation_classes": list(selected),
        "observed": {
            "tandemx_array_union_bp": array_bp,
            "selected_annotation_union_bp": annotation_bp,
            "overlap_union_bp": overlap_bp,
            "class_annotation_bp": class_annotation_bp,
            "class_overlap_bp": class_overlap_bp,
        },
        "checks": checks,
        "complete": all(checks.values()),
    }
    if not receipt["complete"]:
        raise ValueError(f"Independent raster validation disagrees with interval audit: {checks}")
    return receipt


def write_validation_receipt(path: Path, receipt: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
