import json
from pathlib import Path

from benchmarks.scripts.evaluate_ey15_old_new_alignment_audit import (
    aligned_query_segments,
    evaluate,
    intersection_bp,
)


ROOT = Path(__file__).resolve().parents[2]


def test_alignment_audit_is_frozen_before_output_inspection() -> None:
    config = json.loads(
        (ROOT / "benchmarks/configs/ey15_old_new_alignment_audit_v1.json").read_text()
    )
    assert config["status"] == "frozen_before_alignment_output_inspection"
    assert config["preregistration_parent_commit"] == (
        "5608a9176fad039f54685e69225ef238eef31a48"
    )
    assert config["arguments"] == [
        "-x",
        "asm5",
        "-t",
        "4",
        "-c",
        "--eqx",
        "--secondary=yes",
        "-N",
        "100",
    ]
    assert config["evaluation"]["excluded_query_insertion_operation"] == "I"


def test_cigar_query_segments_exclude_insertions_on_both_strands() -> None:
    assert aligned_query_segments(10, 27, "+", "5=3I4X2D5=") == [(10, 15), (18, 27)]
    assert aligned_query_segments(10, 27, "-", "5=3I4X2D5=") == [(10, 19), (22, 27)]
    assert intersection_bp([(0, 10), (8, 20)], [(5, 15)]) == 10


def test_evaluate_reports_all_eligible_families_and_scopes(tmp_path: Path) -> None:
    old_bed = tmp_path / "old.bed"
    new_bed = tmp_path / "new.bed"
    paf = tmp_path / "old_new.paf"
    old_bed.write_text("chr1\t0\t20\tcollapse\nchr1\t0\t100\tretained\n")
    new_bed.write_text("chr1\t0\t100\tcollapse\nchr1\t0\t100\tretained\n")
    paf.write_text(
        "chr1\t100\t0\t100\t+\tchr1\t100\t0\t80\t80\t100\t60\t"
        "tp:A:P\tcg:Z:40=20I40=\n"
    )
    rows, summary = evaluate(
        paf,
        old_bed,
        new_bed,
        min_new_bp=50,
        collapse_threshold=0.6,
    )
    assert len(rows) == 2
    assert summary["eligible_family_count"] == 2
    assert summary["paf"]["primary_rows"] == 1
    for row in rows:
        assert row["same_chromosome_primary_alignment_bp"] == 80
        assert row["same_chromosome_primary_alignment_fraction"] == 0.8
