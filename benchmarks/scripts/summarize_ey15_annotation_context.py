#!/usr/bin/env python3
"""Post hoc biological-context audit for frozen Ey15 reference states.

This does not alter or validate the primary donor-matched classifier.  It asks
whether reference-collapse families defined by the frozen old/new assembly
rule overlap selected repeat annotations supplied by the assembly authors.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_bed_unions(path: Path) -> dict[str, int]:
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED row has fewer than four fields at {path}:{line_number}")
            chrom, family_id = fields[0], fields[3]
            start, end = int(fields[1]), int(fields[2])
            if not chrom or not family_id or start < 0 or end <= start:
                raise ValueError(f"invalid BED interval at {path}:{line_number}")
            intervals.setdefault((family_id, chrom), []).append((start, end))
    totals: dict[str, int] = {}
    for (family_id, _chrom), rows in intervals.items():
        merged = 0
        current_start: int | None = None
        current_end: int | None = None
        for start, end in sorted(rows):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= current_end:
                current_end = max(current_end, end)
            else:
                merged += current_end - current_start
                current_start, current_end = start, end
        if current_start is not None and current_end is not None:
            merged += current_end - current_start
        totals[family_id] = totals.get(family_id, 0) + merged
    return totals


def read_annotation_summary(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "family_id",
            "selected_annotation_overlap_bp",
            "selected_annotation_fraction_of_array",
            "best_annotation_class",
        }
        if not reader.fieldnames or required - set(reader.fieldnames):
            raise ValueError(f"missing annotation fields: {sorted(required - set(reader.fieldnames or []))}")
        result: dict[str, dict[str, str]] = {}
        for line_number, row in enumerate(reader, 2):
            family_id = row["family_id"]
            if not family_id or family_id in result:
                raise ValueError(f"empty or duplicate family_id at {path}:{line_number}")
            result[family_id] = row
    return result


def fisher_right_tail(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher probability for enrichment in the first row."""
    row_total = a + b
    column_total = a + c
    total = a + b + c + d
    denominator = math.comb(total, row_total)
    maximum = min(row_total, column_total)
    return sum(
        math.comb(column_total, value)
        * math.comb(total - column_total, row_total - value)
        / denominator
        for value in range(a, maximum + 1)
    )


def summarize(
    old_arrays: Path,
    new_arrays: Path,
    annotation_summary: Path,
    *,
    collapse_threshold: float,
    overexpansion_threshold: float,
    min_new_bp: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    old_bp = read_bed_unions(old_arrays)
    new_bp = read_bed_unions(new_arrays)
    annotations = read_annotation_summary(annotation_summary)
    rows: list[dict[str, Any]] = []
    for family_id in sorted(set(old_bp) | set(new_bp)):
        old_value = old_bp.get(family_id, 0)
        new_value = new_bp.get(family_id, 0)
        if new_value < min_new_bp:
            continue
        ratio = old_value / new_value
        state = (
            "reference_collapse"
            if ratio < collapse_threshold
            else "old_exceeds_new"
            if ratio > overexpansion_threshold
            else "reference_retained"
        )
        annotation = annotations.get(family_id)
        if annotation is None:
            raise ValueError(f"family absent from annotation summary: {family_id}")
        overlap = int(annotation["selected_annotation_overlap_bp"])
        rows.append(
            {
                "family_id": family_id,
                "old_assembly_bp": old_value,
                "new_assembly_bp": new_value,
                "old_new_ratio": ratio,
                "reference_state": state,
                "selected_annotation_overlap_bp": overlap,
                "selected_annotation_fraction_of_array": float(
                    annotation["selected_annotation_fraction_of_array"]
                ),
                "best_annotation_class": annotation["best_annotation_class"],
                "has_selected_annotation_overlap": overlap > 0,
            }
        )
    collapse_rows = [row for row in rows if row["reference_state"] == "reference_collapse"]
    other_rows = [row for row in rows if row["reference_state"] != "reference_collapse"]
    a = sum(row["has_selected_annotation_overlap"] for row in collapse_rows)
    b = len(collapse_rows) - a
    c = sum(row["has_selected_annotation_overlap"] for row in other_rows)
    d = len(other_rows) - c
    class_counts = Counter(
        row["best_annotation_class"]
        for row in collapse_rows
        if row["has_selected_annotation_overlap"]
    )
    summary = {
        "schema_version": 1,
        "analysis_status": "posthoc_descriptive_audit_after_reference_state_inspection",
        "primary_benchmark_changed": False,
        "config": {
            "collapse_threshold": collapse_threshold,
            "overexpansion_threshold": overexpansion_threshold,
            "min_new_bp": min_new_bp,
            "annotation_positive_rule": "selected_annotation_overlap_bp_greater_than_zero",
        },
        "inputs": {
            "old_arrays": {"path": str(old_arrays), "sha256": sha256(old_arrays)},
            "new_arrays": {"path": str(new_arrays), "sha256": sha256(new_arrays)},
            "annotation_summary": {
                "path": str(annotation_summary),
                "sha256": sha256(annotation_summary),
            },
        },
        "eligible_families": len(rows),
        "reference_collapse_families": len(collapse_rows),
        "other_reference_state_families": len(other_rows),
        "contingency": {
            "collapse_annotated": a,
            "collapse_unannotated": b,
            "other_annotated": c,
            "other_unannotated": d,
        },
        "fisher_exact_greater_p": fisher_right_tail(a, b, c, d),
        "odds_ratio": None if b * c == 0 else a * d / (b * c),
        "odds_ratio_warning": "undefined_due_to_zero_cell" if b * c == 0 else "",
        "collapse_best_annotation_class_counts": dict(sorted(class_counts.items())),
        "warning": (
            "posthoc_association_with_selected_author_annotations_not_independent_truth;"
            "new_assembly_shares_HiFi_evidence_with_read_estimator;"
            "localizer_has_low_author_centromere_annotation_recall"
        ),
    }
    return rows, summary


def write_outputs(rows: list[dict[str, Any]], summary: dict[str, Any], outdir: Path) -> None:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {outdir}")
    outdir.mkdir(parents=True)
    fields = list(rows[0]) if rows else ["family_id"]
    with (outdir / "reference_annotation_rows.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-arrays", required=True, type=Path)
    parser.add_argument("--new-arrays", required=True, type=Path)
    parser.add_argument("--annotation-summary", required=True, type=Path)
    parser.add_argument("--collapse-threshold", required=True, type=float)
    parser.add_argument("--overexpansion-threshold", required=True, type=float)
    parser.add_argument("--min-new-bp", required=True, type=int)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    rows, summary = summarize(
        args.old_arrays,
        args.new_arrays,
        args.annotation_summary,
        collapse_threshold=args.collapse_threshold,
        overexpansion_threshold=args.overexpansion_threshold,
        min_new_bp=args.min_new_bp,
    )
    write_outputs(rows, summary, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
