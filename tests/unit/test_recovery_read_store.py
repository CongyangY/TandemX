from __future__ import annotations

import csv
from pathlib import Path

import pytest

from tandemx.recovery.read_store import collect_recruitment


def _paf(query: str, start: int, end: int, target: str, target_length: int, target_start: int,
         target_end: int, *, matches: int | None = None) -> str:
    block = end - start
    return "\t".join(map(str, (
        query, 10_000, start, end, "+", target, target_length, target_start, target_end,
        block if matches is None else matches, block, 60,
    ))) + "\n"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_streamed_collector_matches_small_oracle(tmp_path: Path) -> None:
    paf = tmp_path / "reads.paf"
    paf.write_text("".join((
        _paf("r1", 100, 700, "repeat_F1", 30_000, 0, 600),
        _paf("r1", 1000, 3000, "L1_left", 2_000, 0, 2_000),
        _paf("r2", 5000, 7000, "L1_right", 2_000, 0, 2_000),
        _paf("r2", 100, 550, "repeat_F1", 30_000, 0, 450),  # below repeat span
        _paf("r3", 1000, 2700, "L1_left", 2_000, 0, 1700),  # below anchor coverage
    )), encoding="utf-8")
    anchors = {
        "L1_left": {"family_id": "F1", "locus_id": "L1"},
        "L1_right": {"family_id": "F1", "locus_id": "L1"},
    }

    repeat_reads, anchor_hits, lengths, row_count = collect_recruitment(
        paf, anchors, tmp_path / "recruited.tsv", {"F1"}, 10, 10
    )

    assert repeat_reads == {"F1": {"r1"}}
    assert set(anchor_hits["L1"]) == {"r1", "r2"}
    assert [hit.target for hit in anchor_hits["L1"]["r1"]] == ["L1_left"]
    assert lengths == {"r1": 10_000, "r2": 10_000}
    assert row_count == 3
    assert [(row["read_id"], row["evidence_type"]) for row in _rows(tmp_path / "recruited.tsv")] == [
        ("r1", "repeat_support_only"), ("r1", "flank_anchored"), ("r2", "flank_anchored"),
    ]


@pytest.mark.parametrize("maximum_rows,maximum_reads,error", [
    (1, 10, "row limit"),
    (10, 1, "read-ID limit"),
])
def test_limits_keep_partial_without_final(tmp_path: Path, maximum_rows: int, maximum_reads: int, error: str) -> None:
    paf = tmp_path / "reads.paf"
    paf.write_text("".join((
        _paf("r1", 100, 700, "repeat_F1", 30_000, 0, 600),
        _paf("r2", 100, 700, "repeat_F1", 30_000, 0, 600),
    )), encoding="utf-8")
    output = tmp_path / "recruited.tsv"

    with pytest.raises(ValueError, match=error):
        collect_recruitment(paf, {}, output, {"F1"}, maximum_reads, maximum_rows)
    assert not output.exists()
    partial = output.with_suffix(".tsv.partial")
    assert partial.is_file()
    assert len(_rows(partial)) == 1


def test_unknown_repeat_family_and_target_are_rejected(tmp_path: Path) -> None:
    paf = tmp_path / "reads.paf"
    output = tmp_path / "recruited.tsv"
    paf.write_text(_paf("r1", 100, 700, "repeat_F2", 30_000, 0, 600), encoding="utf-8")
    with pytest.raises(ValueError, match="unselected repeat family: F2"):
        collect_recruitment(paf, {}, output, {"F1"}, 10, 10)
    assert not output.exists()
    assert output.with_suffix(".tsv.partial").is_file()
