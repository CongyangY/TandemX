"""Summarize pre-registered TRF and TideHunter array calls on one locus.

This is a result parser only.  It does not run either detector or infer
sequence novelty.  Both input formats use one-based, inclusive coordinates
relative to the extracted locus.  Coverage is the union of clipped intervals;
the dominant period is the period with the largest per-period interval union.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Sequence

from tandemx.importers.tidehunter import iter_tidehunter_f2


TARGET_PERIOD_BP = 785
DEFAULT_REGION_LENGTH_BP = 124_029
DEFAULT_COVERAGE_THRESHOLD = 0.90
DEFAULT_PERIOD_TOLERANCE = 0.05

OUTPUT_FIELDS = [
    "tool",
    "tool_version",
    "region_length_bp",
    "record_count",
    "covered_bp",
    "coverage_fraction",
    "dominant_period",
    "period_relation_to_785",
    "passes_coverage",
    "passes_period",
    "verdict",
]


@dataclass(frozen=True)
class ArrayHit:
    """One detector interval in one-based, inclusive locus coordinates."""

    start: int
    end: int
    period: int

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < self.start or self.period < 1:
            raise ValueError(f"Invalid array hit: {self}")


def _parse_int(value: str, path: Path, line_number: int, field: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"Invalid TRF {field} at {path}:{line_number}: {value!r}") from error


def iter_trf_hits(path: Path) -> Iterable[ArrayHit]:
    """Yield TRF ``-d`` data rows while ignoring its human-readable headers.

    TRF data rows have at least the standard 15 columns.  The first three
    columns are start, end and period.  Rows beginning with a nonnumeric token
    are headers; a numeric row with malformed coordinates is an input error.
    """

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith(("@", "#")):
                continue
            fields = stripped.split()
            if not fields or not fields[0].lstrip("+-").isdigit():
                continue
            if len(fields) < 15:
                raise ValueError(
                    f"TRF numeric row has fewer than 15 fields at {path}:{line_number}"
                )
            start = _parse_int(fields[0], path, line_number, "start")
            end = _parse_int(fields[1], path, line_number, "end")
            period = _parse_int(fields[2], path, line_number, "period")
            if start < 1 or end < start or period < 1:
                raise ValueError(f"Invalid TRF coordinates at {path}:{line_number}")
            yield ArrayHit(start, end, period)


def iter_tidehunter_hits(path: Path) -> Iterable[ArrayHit]:
    """Yield TideHunter ``-f 2`` intervals using the repository's strict parser."""

    for record in iter_tidehunter_f2(path):
        yield ArrayHit(record.start, record.end, record.period)


def _clip_interval(hit: ArrayHit, region_length_bp: int) -> tuple[int, int] | None:
    start = max(1, hit.start)
    end = min(region_length_bp, hit.end)
    return (start, end) if start <= end else None


def _union_length(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    covered = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    return covered + current_end - current_start + 1


def period_relation_to_target(
    period: int | None,
    *,
    target_period_bp: int = TARGET_PERIOD_BP,
    tolerance: float = DEFAULT_PERIOD_TOLERANCE,
) -> str:
    """Classify a period as target, integer multiple/divisor, or unresolved."""

    if period is None:
        return "unresolved_no_period"
    if target_period_bp < 1 or not 0 <= tolerance < 1:
        raise ValueError("target period must be positive and tolerance must be in [0,1)")
    if period == target_period_bp:
        return f"exact_{target_period_bp}"
    if abs(period - target_period_bp) / target_period_bp <= tolerance:
        return f"near_{target_period_bp}"

    max_multiple = max(2, math.ceil(max(period, target_period_bp) / min(period, target_period_bp)) + 1)
    for multiple in range(2, max_multiple + 1):
        expected_multiple = target_period_bp * multiple
        if abs(period - expected_multiple) / expected_multiple <= tolerance:
            return f"{multiple}x_{target_period_bp}"
        expected_divisor = target_period_bp / multiple
        if abs(period - expected_divisor) / expected_divisor <= tolerance:
            return f"1/{multiple}_of_{target_period_bp}"
    return "unresolved"


def summarize_hits(
    tool: str,
    tool_version: str,
    hits: Sequence[ArrayHit],
    *,
    region_length_bp: int = DEFAULT_REGION_LENGTH_BP,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
    period_tolerance: float = DEFAULT_PERIOD_TOLERANCE,
    target_period_bp: int = TARGET_PERIOD_BP,
) -> dict[str, object]:
    """Summarize detector hits into the cross-tool result schema."""

    if not tool.strip() or not tool_version.strip():
        raise ValueError("tool and tool_version are required")
    if region_length_bp < 1:
        raise ValueError("region_length_bp must be positive")
    if not 0 < coverage_threshold <= 1:
        raise ValueError("coverage_threshold must be in (0,1]")
    if not 0 <= period_tolerance < 1:
        raise ValueError("period_tolerance must be in [0,1)")

    clipped: list[tuple[ArrayHit, tuple[int, int]]] = []
    for hit in hits:
        interval = _clip_interval(hit, region_length_bp)
        if interval is not None:
            clipped.append((hit, interval))

    all_intervals = [interval for _, interval in clipped]
    covered_bp = _union_length(all_intervals)
    coverage_fraction = covered_bp / region_length_bp
    coverage_by_period: dict[int, list[tuple[int, int]]] = {}
    for hit, interval in clipped:
        coverage_by_period.setdefault(hit.period, []).append(interval)
    dominant_period = (
        max(
            ((period, _union_length(intervals)) for period, intervals in coverage_by_period.items()),
            key=lambda item: (item[1], -item[0]),
        )[0]
        if coverage_by_period
        else None
    )
    relation = period_relation_to_target(
        dominant_period, target_period_bp=target_period_bp, tolerance=period_tolerance
    )
    passes_coverage = coverage_fraction >= coverage_threshold
    passes_period = relation not in {"unresolved", "unresolved_no_period"}
    if not clipped:
        verdict = "no_records"
    elif passes_coverage and passes_period:
        verdict = "confirmed_tandem_structure"
    elif not passes_coverage and not passes_period:
        verdict = "unresolved_coverage_and_period"
    elif not passes_coverage:
        verdict = "unresolved_coverage"
    else:
        verdict = "unresolved_period"
    return {
        "tool": tool,
        "tool_version": tool_version,
        "region_length_bp": region_length_bp,
        "record_count": len(hits),
        "covered_bp": covered_bp,
        "coverage_fraction": f"{coverage_fraction:.6f}",
        "dominant_period": dominant_period if dominant_period is not None else "NA",
        "period_relation_to_785": relation,
        "passes_coverage": str(passes_coverage).lower(),
        "passes_period": str(passes_period).lower(),
        "verdict": verdict,
    }


def parse_and_summarize(
    trf_path: Path,
    tidehunter_path: Path,
    *,
    region_length_bp: int = DEFAULT_REGION_LENGTH_BP,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
    period_tolerance: float = DEFAULT_PERIOD_TOLERANCE,
) -> list[dict[str, object]]:
    """Parse both raw outputs without invoking either detector."""

    trf_hits = list(iter_trf_hits(trf_path))
    tidehunter_hits = list(iter_tidehunter_hits(tidehunter_path))
    return [
        summarize_hits(
            "TRF",
            "4.10.0-rc.2",
            trf_hits,
            region_length_bp=region_length_bp,
            coverage_threshold=coverage_threshold,
            period_tolerance=period_tolerance,
        ),
        summarize_hits(
            "TideHunter",
            "1.5.5",
            tidehunter_hits,
            region_length_bp=region_length_bp,
            coverage_threshold=coverage_threshold,
            period_tolerance=period_tolerance,
        ),
    ]


def write_summary(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trf", required=True, type=Path, help="TRF 4.10.0-rc.2 -d raw output")
    parser.add_argument("--tidehunter", required=True, type=Path, help="TideHunter 1.5.5 -f 2 raw output")
    parser.add_argument("--region-length", type=int, default=DEFAULT_REGION_LENGTH_BP)
    parser.add_argument("--coverage-threshold", type=float, default=DEFAULT_COVERAGE_THRESHOLD)
    parser.add_argument("--period-tolerance", type=float, default=DEFAULT_PERIOD_TOLERANCE)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    for input_path in (args.trf, args.tidehunter):
        if not input_path.is_file():
            parser.error(f"Missing input: {input_path}")
    try:
        rows = parse_and_summarize(
            args.trf,
            args.tidehunter,
            region_length_bp=args.region_length,
            coverage_threshold=args.coverage_threshold,
            period_tolerance=args.period_tolerance,
        )
        write_summary(args.output, rows)
    except ValueError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
