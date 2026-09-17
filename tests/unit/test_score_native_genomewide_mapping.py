"""Reject malformed or ambiguous full-reference PAF records."""

from __future__ import annotations

import pytest

from benchmarks.scripts.score_native_genomewide_mapping import parse_paf_line


def test_primary_paf_coordinates_and_identity() -> None:
    line = "read1\t1000\t10\t990\t-\tChr3\t2000\t500\t1480\t970\t980\t42\ttp:A:P\tcg:Z:980M\n"
    row = parse_paf_line(line)
    assert row["primary"] is True
    assert row["target_start0"] == 500
    assert row["target_end0"] == 1480
    assert row["identity"] == 970 / 980


@pytest.mark.parametrize("line", [
    "read1\t1000\t10\t1001\t+\tChr3\t2000\t500\t1480\t970\t980\t42\ttp:A:P",
    "read1\t1000\t10\t990\t+\tChr3\t2000\t500\t1480\t981\t980\t42\ttp:A:P",
    "read1\t1000\t10\t990\t+\tChr3\t2000\t500\t1480\t970\t980\t42",
    "read1\t1000\t10\t990\t+\tChr3\t2000\t500\t1480\t970\t980\t42\ttp:A:P\ttp:A:S",
])
def test_invalid_paf_is_rejected(line: str) -> None:
    with pytest.raises(ValueError):
        parse_paf_line(line)
