#!/usr/bin/env python3
"""Evaluate frozen Ey15 assembly-alignment coverage at localized families."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Iterable

from benchmarks.challenge.schema import digest_file, write_table


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
SCOPES = (
    "any_alignment",
    "primary_alignment",
    "same_chromosome_primary_alignment",
    "same_chromosome_primary_mapq20_alignment",
)


def merge(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1]:
            result.append((start, end))
        else:
            result[-1] = (result[-1][0], max(result[-1][1], end))
    return result


def interval_bp(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge(intervals))


def intersection_bp(
    left: Iterable[tuple[int, int]], right: Iterable[tuple[int, int]]
) -> int:
    first, second = merge(left), merge(right)
    total = i = j = 0
    while i < len(first) and j < len(second):
        start = max(first[i][0], second[j][0])
        end = min(first[i][1], second[j][1])
        if end > start:
            total += end - start
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def read_family_intervals(path: Path) -> dict[tuple[str, str], list[tuple[int, int]]]:
    result: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED row has fewer than four fields: {path}:{line_number}")
            chrom, family = fields[0], fields[3]
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start or not chrom or not family:
                raise ValueError(f"invalid BED row: {path}:{line_number}")
            result[(family, chrom)].append((start, end))
    return {key: merge(value) for key, value in result.items()}


def family_totals(
    intervals: dict[tuple[str, str], list[tuple[int, int]]]
) -> dict[str, int]:
    result: dict[str, int] = defaultdict(int)
    for (family, _chrom), rows in intervals.items():
        result[family] += interval_bp(rows)
    return dict(result)


def aligned_query_segments(
    qstart: int, qend: int, strand: str, cigar: str
) -> list[tuple[int, int]]:
    tokens = CIGAR.findall(cigar)
    if not tokens or "".join(f"{length}{op}" for length, op in tokens) != cigar:
        raise ValueError(f"invalid or unsupported CIGAR: {cigar}")
    segments: list[tuple[int, int]] = []
    position = qstart if strand == "+" else qend
    for length_text, operation in tokens:
        length = int(length_text)
        consumes_query = operation in {"M", "I", "S", "=", "X"}
        aligned = operation in {"M", "=", "X"}
        if strand == "+":
            if aligned:
                segments.append((position, position + length))
            if consumes_query:
                position += length
        elif strand == "-":
            if aligned:
                segments.append((position - length, position))
            if consumes_query:
                position -= length
        else:
            raise ValueError(f"invalid PAF strand: {strand}")
    expected = qend if strand == "+" else qstart
    if position != expected:
        raise ValueError(
            f"CIGAR query consumption disagrees with PAF span: {position} != {expected}"
        )
    return merge(segments)


def read_alignment_scopes(
    path: Path,
) -> tuple[dict[str, dict[str, list[tuple[int, int]]]], dict[str, int]]:
    scopes: dict[str, dict[str, list[tuple[int, int]]]] = {
        scope: defaultdict(list) for scope in SCOPES
    }
    counters = {"paf_rows": 0, "primary_rows": 0, "secondary_rows": 0}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF row has fewer than twelve fields: {path}:{line_number}")
            qname, qlength, qstart, qend, strand = (
                fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]
            )
            tname, mapq = fields[5], int(fields[11])
            if qstart < 0 or qend <= qstart or qend > qlength:
                raise ValueError(f"invalid PAF query coordinates: {path}:{line_number}")
            tags = {field.split(":", 2)[0]: field for field in fields[12:] if ":" in field}
            if "cg" not in tags or "tp" not in tags:
                raise ValueError(f"PAF row lacks cg/tp tag: {path}:{line_number}")
            cigar = tags["cg"].split(":", 2)[2]
            segments = aligned_query_segments(qstart, qend, strand, cigar)
            primary = tags["tp"] == "tp:A:P"
            counters["paf_rows"] += 1
            counters["primary_rows" if primary else "secondary_rows"] += 1
            scopes["any_alignment"][qname].extend(segments)
            if primary:
                scopes["primary_alignment"][qname].extend(segments)
                if qname == tname:
                    scopes["same_chromosome_primary_alignment"][qname].extend(segments)
                    if mapq >= 20:
                        scopes["same_chromosome_primary_mapq20_alignment"][qname].extend(
                            segments
                        )
    merged_scopes = {
        scope: {chrom: merge(rows) for chrom, rows in by_chrom.items()}
        for scope, by_chrom in scopes.items()
    }
    return merged_scopes, counters


def evaluate(
    paf: Path,
    old_arrays: Path,
    new_arrays: Path,
    *,
    min_new_bp: int,
    collapse_threshold: float,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    old_intervals = read_family_intervals(old_arrays)
    new_intervals = read_family_intervals(new_arrays)
    old_totals, new_totals = family_totals(old_intervals), family_totals(new_intervals)
    coverage, paf_counts = read_alignment_scopes(paf)
    rows: list[dict[str, object]] = []
    for family in sorted(set(old_totals) | set(new_totals)):
        new_bp = new_totals.get(family, 0)
        if new_bp < min_new_bp:
            continue
        old_bp = old_totals.get(family, 0)
        ratio = old_bp / new_bp
        state = "reference_collapse" if ratio < collapse_threshold else "other_reference_state"
        row: dict[str, object] = {
            "family_id": family,
            "old_assembly_bp": old_bp,
            "new_assembly_bp": new_bp,
            "old_new_ratio": ratio,
            "reference_state_binary": state,
        }
        family_rows = {
            chrom: intervals
            for (candidate, chrom), intervals in new_intervals.items()
            if candidate == family
        }
        for scope in SCOPES:
            aligned_bp = sum(
                intersection_bp(intervals, coverage[scope].get(chrom, []))
                for chrom, intervals in family_rows.items()
            )
            row[f"{scope}_bp"] = aligned_bp
            row[f"{scope}_fraction"] = aligned_bp / new_bp
        rows.append(row)
    groups: dict[str, object] = {}
    for state in ("reference_collapse", "other_reference_state"):
        selected = [row for row in rows if row["reference_state_binary"] == state]
        groups[state] = {
            "families": len(selected),
            **{
                f"{scope}_fraction_mean": mean(
                    float(row[f"{scope}_fraction"]) for row in selected
                )
                if selected
                else None
                for scope in SCOPES
            },
            **{
                f"{scope}_fraction_median": median(
                    float(row[f"{scope}_fraction"]) for row in selected
                )
                if selected
                else None
                for scope in SCOPES
            },
        }
    summary = {
        "schema_version": 1,
        "complete": True,
        "analysis_status": "frozen_post_primary_orthogonal_audit",
        "eligible_family_count": len(rows),
        "paf": paf_counts,
        "groups": groups,
        "inputs": {
            "paf": {"path": str(paf), "sha256": digest_file(paf)},
            "old_arrays": {"path": str(old_arrays), "sha256": digest_file(old_arrays)},
            "new_arrays": {"path": str(new_arrays), "sha256": digest_file(new_arrays)},
        },
        "warning": (
            "assembly_alignment_context_not_independent_copy_truth;"
            "insertions_excluded_from_aligned_query_coverage;"
            "repetitive_multimapping_and_mapper_heuristics_remain"
        ),
    }
    return rows, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--old-arrays", required=True, type=Path)
    parser.add_argument("--new-arrays", required=True, type=Path)
    parser.add_argument("--min-new-bp", required=True, type=int)
    parser.add_argument("--collapse-threshold", required=True, type=float)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {args.outdir}")
    rows, summary = evaluate(
        args.paf,
        args.old_arrays,
        args.new_arrays,
        min_new_bp=args.min_new_bp,
        collapse_threshold=args.collapse_threshold,
    )
    args.outdir.mkdir(parents=True)
    write_table(args.outdir / "family_alignment_context.tsv", rows, list(rows[0]))
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

