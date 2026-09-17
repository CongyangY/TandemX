"""Gates for the scoped natural-flank native-read pilot."""

import json
from pathlib import Path

import pytest

from benchmarks.m2_routes.native_read_pilot import (
    primary_paf_array_interval,
    reverse_complement,
    run,
    trim_between_flanks,
)


def test_natural_flank_trim_and_repeated_anchor_abstention() -> None:
    left = "AACCGGTTAACCGGTT"
    right = "TTGGAACCTTGGAACC"
    sequence = "GGG" + left + "ACGTACGT" + right + "TTT"
    result = trim_between_flanks(sequence, left, right)
    assert result["status"] == "ELIGIBLE"
    assert result["array_sequence"] == "ACGTACGT"
    assert result["trim_interval_0based"] == [19, 27]
    assert trim_between_flanks(sequence + left, left, right)["status"] == "ABSTAIN"
    assert trim_between_flanks("", left, right)["status"] == "ABSTAIN"
    assert trim_between_flanks(right + "ACGT" + left, left, right)["reason"] == \
        "flank_order_or_overlap_invalid"


def test_reverse_complement_restores_reference_orientation() -> None:
    sequence = "AAACCGGTTATCGA"
    assert reverse_complement(reverse_complement(sequence)) == sequence


def test_primary_paf_projects_minus_strand_and_rejects_missing(tmp_path: Path) -> None:
    paf = tmp_path / "reads.paf"
    paf.write_text("r1\t110\t5\t105\t-\tC3|target\t100\t0\t100\t100\t100\t60\t"
                   "tp:A:P\tcg:Z:100M\n")
    assert primary_paf_array_interval(paf, "r1", "C3", 20, 80) == [25, 85]
    with pytest.raises(ValueError, match="exactly one"):
        primary_paf_array_interval(paf, "absent", "C3", 20, 80)
    paf.write_text(paf.read_text() * 2)
    with pytest.raises(ValueError, match="exactly one"):
        primary_paf_array_interval(paf, "r1", "C3", 20, 80)


def test_existing_output_fails_closed(tmp_path: Path) -> None:
    destination = tmp_path / "existing"
    destination.mkdir()
    (destination / "sentinel").write_text("retain")
    with pytest.raises(FileExistsError, match="pilot output exists"):
        run(destination)
    assert (destination / "sentinel").read_text() == "retain"


def test_archived_denominator_and_negative_are_explicit() -> None:
    evidence = Path(__file__).resolve().parents[2] / \
        "benchmarks/m2_routes/native_read_evidence_20260917/per_case.jsonl"
    rows = [json.loads(line) for line in evidence.read_text().splitlines()]
    assert len(rows) == 10
    edited = [row for row in rows if row["case_id"] != "C3_source_original"]
    assert len(edited) == 9
    negative = [row for row in edited if row["injected_deleted_bp"] == 0]
    positive = [row for row in edited if row["injected_deleted_bp"] > 0]
    assert [row["case_id"] for row in negative] == ["C3_terminal_000"]
    assert len(positive) == 8
    for route in ("m2_technical", "length_baseline"):
        assert negative[0][route]["state"] == "SUPPORTED"
        assert all(row[route]["state"] == "DISCORDANT" for row in positive)
    assert all(row["biological_audit_state"] == "NOT_EVALUATED_PAIRING_UNVERIFIED"
               for row in rows)
