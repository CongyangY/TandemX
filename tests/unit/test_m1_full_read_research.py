from tandemx.io.sequences import SequenceRecord

from benchmarks.m1_shared_signature.full_read_research import (
    FullReadConfig, classify_read, summarize_records,
)
from benchmarks.m1_shared_signature.run_full_read_dev import score_segments, verify_input


def test_mixed_read_conserves_mass_and_refuses_identical_families() -> None:
    a = "ACGTTGCA" * 10
    b = "TCCGATGA" * 10
    c = "GGGATCCA" * 10
    catalogue = {"a": a, "b": b, "c": c, "c_duplicate": c}
    config = FullReadConfig()
    sequence = a * 4 + b * 4 + c * 4
    segments = classify_read(sequence, catalogue, config)
    assert segments[0]["start"] == 0
    assert segments[-1]["end"] == len(sequence)
    assert all(left["end"] == right["start"] for left, right in zip(segments, segments[1:]))
    assert not any(row["identity"] in ("c", "c_duplicate") and row["status"] == "assigned"
                   for row in segments)
    summary = summarize_records([SequenceRecord(id="one", sequence=sequence)], catalogue, config)
    assert summary["total_read_bp"] == len(sequence)
    assert summary["assigned_by_family"]["a"] > 0
    assert summary["assigned_by_family"]["b"] > 0
    assert summary["ambiguity_groups"]["c+c_duplicate"] > 0
    assert sum(summary[key] for key in ("assigned_read_bp", "ambiguous_read_bp", "unknown_read_bp")) == len(sequence)


def test_empty_and_read_cap_fail_closed() -> None:
    catalogue = {"a": "ACGT" * 20}
    config = FullReadConfig(max_read_bp=100)
    for sequence in ("", "A" * 101):
        try:
            classify_read(sequence, catalogue, config)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid full read accepted")


def test_frozen_input_and_boundary_confusion_are_explicit() -> None:
    manifest = verify_input()
    assert manifest["read_count"] == 18
    truth = [{"read_id": "r", "length_bp": 100, "spans": [
        {"start": 0, "end": 50, "label": "f1"},
        {"start": 50, "end": 100, "label": "background"},
    ]}]
    segments = {"r": [{"start": 0, "end": 80, "status": "assigned", "identity": "f1"},
                       {"start": 80, "end": 100, "status": "unknown", "identity": None}]}
    score = score_segments(segments, truth, ["f1", "decoy_zero"])
    assert score["wrong_family_assigned_bp"] == 30
    assert score["background_false_assigned_bp"] == 30
    assert score["unknown_background_bp"] == 20
    assert score["unknown_positive_bp"] == 0
    assert score["mass_by_status"] == {"assigned": 80, "unknown": 20}
