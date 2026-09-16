"""Experimental occupancy behavior and bounded streaming tests."""
from __future__ import annotations

import random
import json
import tracemalloc

import pytest

from tandemx.io.sequences import SequenceRecord
from benchmarks.m1_shared_signature.occupancy_research import (
    CompetitiveConfig, classify_window, periodic_edit_distance,
    reverse_complement, run_competitive_occupancy, summarize_records,
)
from benchmarks.m1_shared_signature import occupancy_research


def units() -> tuple[str, str]:
    rng = random.Random(817)
    return ("".join(rng.choices("ACGT", k=80)),
            "".join(rng.choices("ACGT", k=80)))


def config(tmp_path, **overrides) -> CompetitiveConfig:
    values = dict(reads=tmp_path / "reads.fa", monomers=tmp_path / "catalogue.fa",
                  genome_size=8000, outdir=tmp_path / "output", haploid_depth=2.0)
    values.update(overrides)
    return CompetitiveConfig(**values)


def test_periodic_edit_supports_phase_rc_and_indel(tmp_path) -> None:
    a, b = units()
    cfg = config(tmp_path)
    assert periodic_edit_distance(a[17:] + a[:17], a) == 0
    assert periodic_edit_distance(a[:40] + "G" + a[40:], a) <= 1
    assert classify_window(reverse_complement(a), {"a": a, "b": b}, cfg) == ("assigned", "a")
    query = a[:40] + "G" + a[40:-1]
    assert classify_window(query, {"a": a, "b": b}, cfg) == ("assigned", "a")


def test_near_family_ambiguity_and_unknown(tmp_path) -> None:
    a, _ = units()
    b = a[:12] + ("A" if a[12] != "A" else "C") + a[13:]
    cfg = config(tmp_path)
    assert classify_window(a, {"a": a, "b": b}, cfg) == ("ambiguous", "a+b")
    assert classify_window(b, {"a": a, "b": b}, cfg) == ("ambiguous", "a+b")
    assert classify_window("N" * 80, {"a": a, "b": b}, cfg) == ("unknown", None)
    assert classify_window(a[:20], {"a": a}, cfg) == ("unknown", None)


def test_python_fallback_retains_edit_objective(tmp_path, monkeypatch) -> None:
    a, _ = units()
    monkeypatch.setattr(occupancy_research, "edlib", None)
    assert periodic_edit_distance(a[:40] + "G" + a[40:], a) <= 1
    assert classify_window(reverse_complement(a), {"a": a}, config(tmp_path)) == ("assigned", "a")


def test_read_bp_conservation_depth_and_long_read(tmp_path) -> None:
    a, b = units()
    cfg = config(tmp_path)
    records = (SequenceRecord(id=f"r{i}", sequence=sequence)
               for i, sequence in enumerate((a * 125, reverse_complement(b), "N" * 80)))
    summary = summarize_records(records, {"a": a, "b": b}, cfg)
    assert summary["read_count"] == 3
    assert summary["assigned_read_bp"] == 10080
    assert summary["unknown_read_bp"] == 80
    assert summary["ambiguous_read_bp"] == 0
    assert sum(summary["assigned_by_family"].values()) == summary["assigned_read_bp"]
    assert summary["haploid_depth"] == 2.0


def test_fixed_window_tail_is_counted_unknown(tmp_path) -> None:
    a, _ = units()
    summary = summarize_records(
        iter([SequenceRecord(id="tail", sequence=a + a[:20])]),
        {"a": a}, config(tmp_path),
    )
    assert summary["assigned_read_bp"] == 80
    assert summary["unknown_read_bp"] == 20
    assert summary["total_read_bp"] == 100


def test_empty_and_read_cap_fail_without_silent_clipping(tmp_path) -> None:
    a, _ = units()
    cfg = config(tmp_path, max_read_bp=80)
    with pytest.raises(ValueError, match="No reads"):
        summarize_records(iter(()), {"a": a}, cfg)
    with pytest.raises(ValueError, match="read-length cap"):
        summarize_records(iter([SequenceRecord(id="long", sequence=a * 2)]), {"a": a}, cfg)


def test_research_output_receipt_and_failure_cleanup(tmp_path) -> None:
    a, b = units()
    cfg = config(tmp_path)
    cfg.monomers.write_text(f">family_id=a\n{a}\n>family_id=b\n{b}\n")
    cfg.reads.write_text(f">r1\n{a}\n>r2\n{'N' * 80}\n")
    rows = run_competitive_occupancy(cfg)
    assert rows[0]["assigned_read_bp"] == 80
    receipt = json.loads((cfg.outdir / "competitive_occupancy_summary.json").read_text())
    assert receipt["complete"] and receipt["assigned_read_bp"] == 80
    assert receipt["unknown_read_bp"] == 80
    assert (cfg.outdir / "competitive_occupancy.tsv").read_text().startswith("family_id\t")
    empty_cfg = config(tmp_path, reads=tmp_path / "empty.fa", outdir=tmp_path / "empty-output")
    empty_cfg.reads.write_text("")
    with pytest.raises(ValueError):
        run_competitive_occupancy(empty_cfg)
    assert not (empty_cfg.outdir / "competitive_occupancy_summary.json").exists()


def test_catalogue_cap_is_explicit(tmp_path) -> None:
    a, _ = units()
    cfg = config(tmp_path)
    cfg.monomers.write_text("".join(f">family_id=f{i}\n{a}\n" for i in range(9)))
    cfg.reads.write_text(f">r1\n{a}\n")
    with pytest.raises(ValueError, match="family cap"):
        run_competitive_occupancy(cfg)
    assert not cfg.outdir.exists()


def test_streaming_peak_does_not_scale_with_read_count(tmp_path) -> None:
    a, _ = units()
    cfg = config(tmp_path)

    def peak(count: int) -> int:
        tracemalloc.start()
        try:
            records = (SequenceRecord(id=str(i), sequence=a) for i in range(count))
            summary = summarize_records(records, {"a": a}, cfg)
            assert summary["read_count"] == count
            return tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()

    small, larger = peak(5), peak(50)
    assert larger < small + 250_000
