import csv
from pathlib import Path

import pytest

from benchmarks.scripts.confirm_array_tools import (
    ArrayHit,
    iter_trf_hits,
    parse_and_summarize,
    period_relation_to_target,
    summarize_hits,
    write_summary,
)


def _trf_row(start: int, end: int, period: int) -> str:
    return (
        f"{start} {end} {period} 10.0 {period} 99 0 1000 25 25 25 25 2.00 "
        "ACGT ACGT"
    )


def _tidehunter_row(start: int, end: int, period: int, read_length: int) -> str:
    consensus = "A" * period
    return f"read1 repeat1 10.0 {read_length} {start} {end} {period} 99.0 1 {start},{min(end, start + period)} {consensus}"


def test_period_relation_accepts_target_multiples_and_divisors() -> None:
    assert period_relation_to_target(785) == "exact_785"
    assert period_relation_to_target(800) == "near_785"
    assert period_relation_to_target(1570) == "2x_785"
    assert period_relation_to_target(392) == "1/2_of_785"
    assert period_relation_to_target(500) == "unresolved"
    assert period_relation_to_target(None) == "unresolved_no_period"


def test_summarize_hits_unions_intervals_and_uses_period_coverage() -> None:
    row = summarize_hits(
        "TRF",
        "4.10.0-rc.2",
        [
            ArrayHit(1, 100, 785),
            ArrayHit(90, 200, 785),
            ArrayHit(200, 400, 1570),
        ],
        region_length_bp=400,
    )
    assert row["record_count"] == 3
    assert row["covered_bp"] == 400
    assert row["coverage_fraction"] == "1.000000"
    assert row["dominant_period"] == 1570
    assert row["period_relation_to_785"] == "2x_785"
    assert row["passes_coverage"] == "true"
    assert row["passes_period"] == "true"
    assert row["verdict"] == "confirmed_tandem_structure"


def test_parse_raw_trf_and_tidehunter_outputs_and_write_summary(tmp_path: Path) -> None:
    trf = tmp_path / "trf.txt"
    trf.write_text(
        "Tandem Repeats Finder\n"
        "Sequence: TXF000708_array\n"
        f"{_trf_row(1, 50, 785)}\n"
        f"{_trf_row(40, 124029, 785)}\n"
    )
    tidehunter = tmp_path / "tidehunter.tsv"
    tidehunter.write_text(_tidehunter_row(1, 124029, 785, 124029) + "\n")
    rows = parse_and_summarize(trf, tidehunter)
    assert [row["tool"] for row in rows] == ["TRF", "TideHunter"]
    assert all(row["covered_bp"] == 124029 for row in rows)
    assert all(row["dominant_period"] == 785 for row in rows)
    assert all(row["passes_coverage"] == "true" for row in rows)
    assert all(row["passes_period"] == "true" for row in rows)

    output = tmp_path / "summary.tsv"
    write_summary(output, rows)
    parsed = list(csv.DictReader(output.open(), delimiter="\t"))
    assert parsed[0]["tool_version"] == "4.10.0-rc.2"
    assert parsed[1]["tool_version"] == "1.5.5"
    assert parsed[0]["verdict"] == "confirmed_tandem_structure"


def test_trf_numeric_short_row_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "bad_trf.txt"
    path.write_text("1 10 785\n")
    with pytest.raises(ValueError, match="fewer than 15"):
        list(iter_trf_hits(path))


def test_empty_detector_output_is_explicit(tmp_path: Path) -> None:
    trf = tmp_path / "empty_trf.txt"
    trf.write_text("Tandem Repeats Finder\n")
    tidehunter = tmp_path / "empty_tidehunter.tsv"
    tidehunter.write_text("# no calls\n")
    rows = parse_and_summarize(trf, tidehunter)
    assert all(row["verdict"] == "no_records" for row in rows)

    row = summarize_hits("TideHunter", "1.5.5", [], region_length_bp=100)
    assert row["covered_bp"] == 0
    assert row["dominant_period"] == "NA"
    assert row["passes_coverage"] == "false"
    assert row["passes_period"] == "false"
    assert row["verdict"] == "no_records"
