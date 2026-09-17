from tandemx.io.sequences import SequenceRecord

from benchmarks.m1_shared_signature.full_read_research import (
    FullReadConfig, classify_read, summarize_records,
)


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
