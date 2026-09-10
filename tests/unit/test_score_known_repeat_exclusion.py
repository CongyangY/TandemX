from __future__ import annotations

from benchmarks.scripts.score_known_repeat_exclusion import (
    circular_glocal_identity,
    read_masked_known_fasta,
    score_pair,
)


def test_circular_glocal_identity_handles_rotation_and_reverse_complement() -> None:
    identity, distance, orientation = circular_glocal_identity("AAGTC", "TTGAC")
    assert identity == 1.0
    assert distance == 0
    assert orientation == "reverse"


def test_score_pair_marks_an_integer_multiple_known_match() -> None:
    row = score_pair(
        "ACGT" * 20,
        "GTAC" * 10,
        possible_threshold=0.80,
        strong_threshold=0.90,
        minimum_aligned_bp=20,
    )
    assert row["glocal_edit_identity"] == 1.0
    assert row["integer_multiple"] == 2
    assert row["exclusion_state"] == "strong_known_clone_match"


def test_read_masked_known_fasta_records_iupac_masking(tmp_path) -> None:
    fasta = tmp_path / "known.fa"
    fasta.write_text(">rDNA\nACGTWYNR\n", encoding="utf-8")
    records, masked_count = read_masked_known_fasta(fasta)
    assert masked_count == 4
    assert records[0].identifier == "rDNA"
    assert records[0].sequence == "ACGTNNNN"
