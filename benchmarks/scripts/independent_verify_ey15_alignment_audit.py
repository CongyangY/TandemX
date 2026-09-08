#!/usr/bin/env python3
"""Independently verify Ey15 assembly-alignment family coverage.

The verifier does not import the primary alignment evaluator.  It reparses BED
and PAF, manually decodes CIGAR, and marks coverage inside each eligible family
interval with byte arrays before comparing every reported field and summary.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Any, Iterable


SCOPES = (
    "any_alignment",
    "primary_alignment",
    "same_chromosome_primary_alignment",
    "same_chromosome_primary_mapq20_alignment",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def merged(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[list[int]] = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(start, end) for start, end in result]


def read_bed(path: Path) -> dict[tuple[str, str], list[tuple[int, int]]]:
    raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"invalid BED row: {path}:{line_number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start or not fields[0] or not fields[3]:
                raise ValueError(f"invalid BED interval: {path}:{line_number}")
            raw[(fields[3], fields[0])].append((start, end))
    return {key: merged(rows) for key, rows in raw.items()}


def totals(rows: dict[tuple[str, str], list[tuple[int, int]]]) -> dict[str, int]:
    result: dict[str, int] = defaultdict(int)
    for (family, _chrom), intervals in rows.items():
        result[family] += sum(end - start for start, end in intervals)
    return dict(result)


def cigar_tokens(value: str) -> list[tuple[int, str]]:
    tokens: list[tuple[int, str]] = []
    digits = ""
    for character in value:
        if character.isdigit():
            digits += character
        else:
            if not digits or character not in "MIDNSHP=X":
                raise ValueError(f"invalid CIGAR: {value}")
            tokens.append((int(digits), character))
            digits = ""
    if digits or not tokens:
        raise ValueError(f"invalid CIGAR: {value}")
    return tokens


def query_segments(qstart: int, qend: int, strand: str, cigar: str) -> list[tuple[int, int]]:
    position = qstart if strand == "+" else qend
    result: list[tuple[int, int]] = []
    for length, operation in cigar_tokens(cigar):
        query_consuming = operation in "MIS=X"
        aligned = operation in "M=X"
        if strand == "+":
            if aligned:
                result.append((position, position + length))
            if query_consuming:
                position += length
        elif strand == "-":
            if aligned:
                result.append((position - length, position))
            if query_consuming:
                position -= length
        else:
            raise ValueError(f"invalid strand: {strand}")
    if position != (qend if strand == "+" else qstart):
        raise ValueError("CIGAR query consumption differs from PAF query span")
    return result


def parse_paf(path: Path) -> tuple[list[dict[str, Any]], dict[str, int]]:
    alignments: list[dict[str, Any]] = []
    counts = {"paf_rows": 0, "primary_rows": 0, "secondary_rows": 0}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"invalid PAF row: {path}:{line_number}")
            tags = {field.partition(":")[0]: field for field in fields[12:]}
            if not tags.get("cg", "").startswith("cg:Z:") or not tags.get(
                "tp", ""
            ).startswith("tp:A:"):
                raise ValueError(f"missing cg/tp tag: {path}:{line_number}")
            qstart, qend = int(fields[2]), int(fields[3])
            primary = tags["tp"] == "tp:A:P"
            alignments.append(
                {
                    "query": fields[0],
                    "target": fields[5],
                    "mapq": int(fields[11]),
                    "primary": primary,
                    "segments": query_segments(
                        qstart, qend, fields[4], tags["cg"].split(":", 2)[2]
                    ),
                }
            )
            counts["paf_rows"] += 1
            counts["primary_rows" if primary else "secondary_rows"] += 1
    return alignments, counts


def mark_overlap(
    interval: tuple[int, int],
    segments: Iterable[tuple[int, int]],
    markers: bytearray,
) -> None:
    interval_start, interval_end = interval
    for segment_start, segment_end in segments:
        start = max(interval_start, segment_start)
        end = min(interval_end, segment_end)
        if end > start:
            markers[start - interval_start : end - interval_start] = b"\x01" * (end - start)


def expected(
    paf: Path,
    old_arrays: Path,
    new_arrays: Path,
    *,
    min_new_bp: int,
    collapse_threshold: float,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    old_rows, new_rows = read_bed(old_arrays), read_bed(new_arrays)
    old_totals, new_totals = totals(old_rows), totals(new_rows)
    families = sorted(
        family for family, value in new_totals.items() if value >= min_new_bp
    )
    alignments, counts = parse_paf(paf)
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for alignment in alignments:
        by_chrom[alignment["query"]].append(alignment)
    rows: dict[str, dict[str, Any]] = {}
    for family in families:
        old_bp, new_bp = old_totals.get(family, 0), new_totals[family]
        row: dict[str, Any] = {
            "family_id": family,
            "old_assembly_bp": old_bp,
            "new_assembly_bp": new_bp,
            "old_new_ratio": old_bp / new_bp,
            "reference_state_binary": (
                "reference_collapse"
                if old_bp / new_bp < collapse_threshold
                else "other_reference_state"
            ),
        }
        covered = {scope: 0 for scope in SCOPES}
        for (candidate, chrom), intervals in new_rows.items():
            if candidate != family:
                continue
            for interval in intervals:
                markers = {
                    scope: bytearray(interval[1] - interval[0]) for scope in SCOPES
                }
                for alignment in by_chrom.get(chrom, []):
                    scopes = ["any_alignment"]
                    if alignment["primary"]:
                        scopes.append("primary_alignment")
                        if alignment["target"] == chrom:
                            scopes.append("same_chromosome_primary_alignment")
                            if alignment["mapq"] >= 20:
                                scopes.append("same_chromosome_primary_mapq20_alignment")
                    for scope in scopes:
                        mark_overlap(interval, alignment["segments"], markers[scope])
                for scope in SCOPES:
                    covered[scope] += sum(markers[scope])
        for scope in SCOPES:
            row[f"{scope}_bp"] = covered[scope]
            row[f"{scope}_fraction"] = covered[scope] / new_bp
        rows[family] = row
    groups: dict[str, Any] = {}
    for state in ("reference_collapse", "other_reference_state"):
        selected = [row for row in rows.values() if row["reference_state_binary"] == state]
        groups[state] = {
            "families": len(selected),
            **{
                f"{scope}_fraction_mean": (
                    mean(row[f"{scope}_fraction"] for row in selected)
                    if selected
                    else None
                )
                for scope in SCOPES
            },
            **{
                f"{scope}_fraction_median": (
                    median(row[f"{scope}_fraction"] for row in selected)
                    if selected
                    else None
                )
                for scope in SCOPES
            },
        }
    return rows, {"eligible_family_count": len(rows), "paf": counts, "groups": groups}


def close(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return left is None and right is None
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(float(left), float(right), rel_tol=1e-9, abs_tol=1e-9)
    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(close(left[key], right[key]) for key in left)
    return left == right


def parse_value(field: str, value: str) -> Any:
    if field in {"family_id", "reference_state_binary"}:
        return value
    if field.endswith("_bp"):
        return int(value)
    return float(value)


def verify(
    paf: Path,
    old_arrays: Path,
    new_arrays: Path,
    family_context: Path,
    summary_path: Path,
    *,
    min_new_bp: int,
    collapse_threshold: float,
) -> dict[str, Any]:
    expected_rows, expected_summary = expected(
        paf,
        old_arrays,
        new_arrays,
        min_new_bp=min_new_bp,
        collapse_threshold=collapse_threshold,
    )
    with family_context.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed_rows = {
            row["family_id"]: {key: parse_value(key, value) for key, value in row.items()}
            for row in reader
        }
    failures: list[str] = []
    if set(expected_rows) != set(observed_rows):
        failures.append("family_universe_mismatch")
    for family in sorted(set(expected_rows) & set(observed_rows)):
        for field, expected_value in expected_rows[family].items():
            if field not in observed_rows[family] or not close(
                expected_value, observed_rows[family][field]
            ):
                failures.append(f"family_field_mismatch:{family}:{field}")
    observed_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    for field, expected_value in expected_summary.items():
        if field not in observed_summary or not close(expected_value, observed_summary[field]):
            failures.append(f"summary_field_mismatch:{field}")
    paths = {
        "paf": paf,
        "old_arrays": old_arrays,
        "new_arrays": new_arrays,
        "family_context": family_context,
        "summary": summary_path,
    }
    return {
        "schema_version": 1,
        "verification_method": "independent_manual_CIGAR_parser_and_per_interval_byte_marking",
        "verification_passed": not failures,
        "eligible_families_recomputed": len(expected_rows),
        "failures": failures,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--old-arrays", required=True, type=Path)
    parser.add_argument("--new-arrays", required=True, type=Path)
    parser.add_argument("--family-context", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--min-new-bp", required=True, type=int)
    parser.add_argument("--collapse-threshold", required=True, type=float)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(f"refusing to overwrite receipt: {args.receipt}")
    result = verify(
        args.paf,
        args.old_arrays,
        args.new_arrays,
        args.family_context,
        args.summary,
        min_new_bp=args.min_new_bp,
        collapse_threshold=args.collapse_threshold,
    )
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0 if result["verification_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

