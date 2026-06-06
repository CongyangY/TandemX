from __future__ import annotations

from pathlib import Path

from tandemx.discover.mvp import FastaRecord
from tandemx.locate.mvp import (
    ArrayHit,
    classify_assembly_read_ratio,
    compare_assembly_to_reads,
    covered_bp,
    merge_intervals,
    window_density,
)


def test_merge_intervals_with_gap() -> None:
    assert merge_intervals([(0, 10), (12, 20), (40, 50)], max_gap=2) == [
        (0, 20, 2),
        (40, 50, 1),
    ]


def test_window_density() -> None:
    record = FastaRecord(read_id="chr1", description="chr1", sequence="A" * 100)
    windows = window_density(record, [(10, 30), (50, 70)], window_size=50, step_size=50)
    assert windows[0].score == 0.4
    assert windows[1].score == 0.4
    assert covered_bp([(10, 30)], 0, 20) == 10


def test_classify_assembly_read_ratio() -> None:
    assert classify_assembly_read_ratio(1000, 300)[0] == "possible_collapse"
    assert classify_assembly_read_ratio(1000, 1800)[0] == "possible_overexpansion"
    assert classify_assembly_read_ratio(1000, 900)[0] == "consistent"
    assert classify_assembly_read_ratio(0, 0)[0] == "low_confidence"


def test_compare_assembly_to_reads(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    copy_number.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
        "haploid_depth\testimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF000001\t100\t10\t10\t1\t10\t1000\thigh\t\n",
        encoding="utf-8",
    )
    comparisons = compare_assembly_to_reads(
        [
            ArrayHit("chr1", 0, 300, "TXF000001", 900, ".", "high", ""),
        ],
        copy_number,
    )
    assert comparisons[0].status == "possible_collapse"
    assert comparisons[0].assembly_estimated_bp == 300
