"""Independent checks for the frozen engineered-edit truth and scorer gates."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.scripts.build_b1_structure_cases import build, load_protocol, reverse_complement, sha256_file, write_jsonl
from benchmarks.scripts.score_b1_structure_cases import read_jsonl, score


ROOT = Path(__file__).resolve().parents[2]
V2 = ROOT / "benchmarks/controlled_collapse/v2"


def _bundle(tmp_path: Path) -> Path:
    out = tmp_path / "bundle"
    build(V2 / "protocol.json", ROOT, out)
    return out


def test_frozen_source_hash_and_reproducible_bundle(tmp_path: Path) -> None:
    protocol, motifs = load_protocol(V2 / "protocol.json", ROOT)
    assert len(motifs) == 4 and {len(x) for x in motifs.values()} == {24}
    assert sha256_file(ROOT / protocol["source_template_path"]) == protocol["source_template_sha256"]
    out = _bundle(tmp_path)
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["inputs_sha256"] == sha256_file(V2 / "development_bundle/inputs.jsonl")
    assert receipt["truth_sha256"] == sha256_file(V2 / "development_bundle/truth.jsonl")
    assert receipt["case_count"] == 13


def test_exact_edited_sequences_and_truth_spans(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path)
    inputs = read_jsonl(bundle / "inputs.jsonl")
    truths = read_jsonl(bundle / "truth.jsonl")
    by_id = {row["case_id"]: row for row in inputs}
    assert len(by_id) == len(truths) == 13
    # Repeated identical HOR copies make some distinct edit coordinates
    # observationally identical from sequence alone.
    assert len({row["assembly_sequence"] for row in inputs}) == 10
    assert {row["event_type"] for row in truths} == {
        "intact", "copy_loss", "complete_deletion", "internal_deletion",
        "boundary_truncation", "HOR_copy_deletion", "monomer_deletion",
        "duplication", "rearrangement", "inversion", "compression"}
    assert sum(row["event_status"] == "intact_negative" for row in truths) == 2
    for truth in truths:
        public = by_id[truth["case_id"]]
        assert not set(public).intersection({"event_status", "event_type", "signed_bp_delta", "edited_label_path"})
        left, right = public["left_flank_sequence"], public["right_flank_sequence"]
        assembly = public["assembly_sequence"]
        assert assembly.startswith(left) and assembly.endswith(right)
        assert assembly.count(left) == assembly.count(right) == 1
        assert truth["edited_array_interval_bp"] == [len(left), len(assembly) - len(right)]
        assert truth["signed_bp_delta"] == len(assembly) - len(next(iter(public["raw_read_sequences"].values())))
        assert sum(b - a for a, b in truth["source_deleted_spans_bp"]) == truth["removed_bp"]
        assert sum(b - a for a, b in truth["edited_inserted_spans_bp"]) == truth["inserted_bp"]
        # Reconstruct each edited segment from the immutable source read and the
        # declared operation; this catches coordinate or orientation mistakes.
        source = next(iter(public["raw_read_sequences"].values()))
        rebuilt = [None] * len(assembly)
        for segment in truth["coordinate_chain"]:
            a, b = segment["edited_start0"], segment["edited_end0"]
            if segment["operation"] == "delete":
                assert a == b
                continue
            original = source[segment["source_start0"]:segment["source_end0"]]
            if segment["operation"] == "reverse_complement":
                original = reverse_complement(original)
            elif segment["operation"] == "compress_tail":
                original = original[:b - a]
            assert len(original) == b - a
            assert all(x is None for x in rebuilt[a:b])
            rebuilt[a:b] = original
        assert "".join(rebuilt) == assembly
    complete = next(row for row in truths if row["event_type"] == "complete_deletion")
    assert complete["edited_copy_count"] == 0 and complete["signed_bp_delta"] == -384
    inversion = next(row for row in truths if row["event_type"] == "inversion")
    assert inversion["signed_bp_delta"] == 0 and inversion["edited_orientation_path"].count("-") == 1
    compression = next(row for row in truths if row["event_type"] == "compression")
    assert compression["signed_bp_delta"] == -32 and compression["edited_copy_count"] == 16


def _perfect_predictions(bundle: Path, path: Path) -> None:
    rows = []
    for truth in read_jsonl(bundle / "truth.jsonl"):
        rows.append(dict(case_id=truth["case_id"], status="ok",
                         event_score=0.0 if truth["event_status"] == "intact_negative" else 1.0,
                         event_type=truth["event_type"],
                         predicted_edited_label_path=truth["edited_label_path"],
                         predicted_edited_orientation_path=truth["edited_orientation_path"],
                         predicted_signed_bp_delta=truth["signed_bp_delta"],
                         predicted_edited_interval_bp=truth["edited_array_interval_bp"]))
    write_jsonl(path, rows)


def test_scorer_perfect_is_tautological_smoke(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)  # derived from truth: scorer test, never route performance
    summary = score(bundle, V2 / "metrics_frozen.json", pred, tmp_path / "score")
    assert summary["primary_status"] == "complete"
    assert summary["primary"] == {"tp": 11, "fn": 0, "fp": 0, "tn": 2,
                                  "sensitivity": 1.0, "fpr": 0.0, "ppv": 1.0}
    assert summary["secondary"]["signed_bp_delta_mae"] == 0
    assert summary["secondary"]["edited_label_path_mean_distance"] == 0
    assert summary["secondary"]["edited_array_interval_mean_iou"] == 1
    assert summary["observable_equivalence_class_count"] == 10
    assert summary["event_type_nonidentifiable_case_count"] == 3
    assert summary["secondary"]["event_type_accuracy"] is None
    ambiguous = [row for row in read_jsonl(tmp_path / "score/per_case.jsonl")
                 if row["event_type_identifiability_status"] == "not_identifiable_from_available_sequences"]
    assert {row["truth_event_type"] for row in ambiguous} == {
        "copy_loss", "HOR_copy_deletion", "boundary_truncation"}
    assert all(row["event_type_match"] is None for row in ambiguous)


def test_abstain_and_missing_fields_block_aggregate(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)
    rows = read_jsonl(pred)
    rows[0].update(status="abstain", event_score=None)
    rows[1]["predicted_signed_bp_delta"] = None
    write_jsonl(pred, rows)
    summary = score(bundle, V2 / "metrics_frozen.json", pred, tmp_path / "score")
    assert summary["denominator"] == 13 and summary["status_counts"]["abstain"] == 1
    assert summary["primary_status"] == "blocked_incomplete_predictions"
    assert all(value is None for value in summary["primary"].values())
    assert all(value is None for value in summary["secondary"].values())


@pytest.mark.parametrize("bad", [None, float("nan"), float("inf"), -0.1, 1.1])
def test_ok_score_must_be_finite(tmp_path: Path, bad: float | None) -> None:
    bundle = _bundle(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)
    rows = read_jsonl(pred)
    rows[0]["event_score"] = bad
    write_jsonl(pred, rows)
    with pytest.raises(ValueError, match="finite event_score"):
        score(bundle, V2 / "metrics_frozen.json", pred, tmp_path / "score")


def test_hash_tamper_fails_closed(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)
    (bundle / "truth.jsonl").write_text((bundle / "truth.jsonl").read_text() + "\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        score(bundle, V2 / "metrics_frozen.json", pred, tmp_path / "score")
