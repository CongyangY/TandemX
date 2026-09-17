"""Synthetic source, injected edits, noisy reads and frozen B2 scorer checks."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.scripts.build_b1_structure_cases import reverse_complement, sha256_file, write_jsonl
from benchmarks.scripts.build_b2_cases import build, load_source
from benchmarks.scripts.generate_b2_source import make_source
from benchmarks.scripts.score_b1_structure_cases import read_jsonl
from benchmarks.scripts.score_b2_cases import score
from benchmarks.scripts.score_b1_colcen_cases import edited_breakpoints


ROOT = Path(__file__).resolve().parents[2]
B2 = ROOT / "benchmarks/controlled_collapse/b2"


def _build(tmp_path: Path) -> Path:
    out = tmp_path / "bundle"
    build(B2 / "protocol.json", B2 / "source_manifest.json", out)
    return out


def _perfect_predictions(bundle: Path, path: Path) -> None:
    rows = []
    for truth in read_jsonl(bundle / "truth.jsonl"):
        rows.append(dict(case_id=truth["case_id"], status="ok",
                         event_score=1.0 if truth["event_status"] == "positive" else 0.0,
                         event_type=truth["event_type"],
                         predicted_signed_bp_delta=truth["signed_bp_delta"],
                         predicted_edited_label_path=truth["edited_label_path"],
                         predicted_HOR_period_monomers=truth["source_HOR_period_monomers"],
                         predicted_breakpoints_bp=edited_breakpoints(truth)))
    write_jsonl(path, rows)


def test_independent_source_manifest_replays(tmp_path: Path) -> None:
    manifest = make_source(B2 / "protocol.json", tmp_path / "source.json")
    assert sha256_file(tmp_path / "source.json") == sha256_file(B2 / "source_manifest.json")
    assert manifest["split"] == "synthetic_held_out"
    assert [row["monomer_length_bp"] for row in manifest["lineages"]] == [100, 200, 300, 400]
    assert all(row["source_label_path"][10:14] == list("BCBC") for row in manifest["lineages"])
    assert all(len(row["source_copies"]) == 24 for row in manifest["lineages"])
    protocol, lineages = load_source(B2 / "protocol.json", B2 / "source_manifest.json")
    assert len(lineages) == 4 and protocol["source_seed"] != 17862026


def test_case_bundle_replays_and_edit_ledger_reconstructs(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    receipt = json.loads((bundle / "receipt.json").read_text())
    frozen = json.loads((B2 / "held_out_bundle/receipt.json").read_text())
    assert receipt["inputs_sha256"] == frozen["inputs_sha256"] == sha256_file(bundle / "inputs.jsonl")
    assert receipt["truth_sha256"] == frozen["truth_sha256"] == sha256_file(bundle / "truth.jsonl")
    inputs = {row["case_id"]: row for row in read_jsonl(bundle / "inputs.jsonl")}
    truths = read_jsonl(bundle / "truth.jsonl")
    assert len(inputs) == len(truths) == 52
    assert len({row["assembly_sequence"] for row in inputs.values()}) == 48
    assert sum(row["event_status"] == "positive" for row in truths) == 44
    assert {row["monomer_length_bp"] for row in truths} == {100, 200, 300, 400}
    assert {row["synthetic_read_support"] for row in inputs.values()} == {1, 3}
    for truth in truths:
        public = inputs[truth["case_id"]]
        assert public["split"] == "synthetic_held_out"
        assert not set(public).intersection({"event_status", "event_type", "signed_bp_delta", "edited_label_path"})
        assert truth["unit_truth_status"] == "exact_synthetic"
        assert truth["source_HOR_period_monomers"] == 5
        assert truth["source_label_path"][10:14] == list("BCBC")
        source = next(row for row in json.loads((B2 / "source_manifest.json").read_text())["lineages"]
                      if row["lineage_id"] == truth["material_id"])
        source_sequence = public["left_flank_sequence"] + "".join(x["sequence"] for x in source["source_copies"]) + public["right_flank_sequence"]
        edited = public["assembly_sequence"]
        assert truth["signed_bp_delta"] == len(edited) - len(source_sequence)
        rebuilt = [None] * len(edited)
        for segment in truth["coordinate_chain"]:
            a, b = segment["edited_start0"], segment["edited_end0"]
            if segment["operation"] == "delete":
                assert a == b
                continue
            original = source_sequence[segment["source_start0"]:segment["source_end0"]]
            if segment["operation"] == "reverse_complement":
                original = reverse_complement(original)
            elif segment["operation"] == "compress_tail":
                original = original[:b - a]
            assert len(original) == b - a
            assert all(x is None for x in rebuilt[a:b])
            rebuilt[a:b] = original
        assert "".join(rebuilt) == edited


def test_independent_read_error_ledger_and_support(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    inputs = {row["case_id"]: row for row in read_jsonl(bundle / "inputs.jsonl")}
    truths = read_jsonl(bundle / "truth.jsonl")
    manifest = {row["lineage_id"]: row for row in json.loads((B2 / "source_manifest.json").read_text())["lineages"]}
    noisy_cases = 0
    for truth in truths:
        public = inputs[truth["case_id"]]
        source_array = "".join(x["sequence"] for x in manifest[truth["material_id"]]["source_copies"])
        left, right = public["left_flank_sequence"], public["right_flank_sequence"]
        assert len(public["raw_read_sequences"]) == public["synthetic_read_support"]
        for read_id, observed in public["raw_read_sequences"].items():
            assert observed.startswith(left) and observed.endswith(right)
            by_position = {}
            for event in truth["read_error_events"][read_id]:
                by_position.setdefault(event["source_position0"] - len(left), []).append(event)
            rebuilt = []
            for position, base in enumerate(source_array):
                current = base
                deleted = False
                for event in by_position.get(position, []):
                    if event["operation"] == "insert_before":
                        rebuilt.append(event["inserted_base"])
                    elif event["operation"] == "delete":
                        deleted = True
                    elif event["operation"] == "substitute":
                        current = event["replacement_base"]
                if not deleted:
                    rebuilt.append(current)
            assert "".join(rebuilt) == observed[len(left):-len(right)]
        if public["error_profile"]["substitution_rate"] > 0:
            noisy_cases += 1
        else:
            assert all(not events for events in truth["read_error_events"].values())
    assert noisy_cases == 44


def test_tautological_scorer_and_abstention_cost(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)
    full = score(bundle, B2 / "metrics_frozen.json", pred, tmp_path / "full")
    assert full["primary_intent_to_diagnose"]["tp"] == 44
    assert full["primary_intent_to_diagnose"]["tn"] == 8
    assert full["signed_bp_delta_mae"] == 0
    assert full["breakpoint_localization_error_bp"] == 0
    assert full["HOR_period_error_monomers"] == 0
    assert full["repeat_order_accuracy"] == 1
    rows = read_jsonl(pred)
    truth = {row["case_id"]: row for row in read_jsonl(bundle / "truth.jsonl")}
    positive = next(row for row in rows if truth[row["case_id"]]["event_status"] == "positive")
    intact = next(row for row in rows if truth[row["case_id"]]["event_status"] == "intact_negative")
    for row in (positive, intact):
        row.update(status="abstain", event_score=None)
    write_jsonl(pred, rows)
    partial = score(bundle, B2 / "metrics_frozen.json", pred, tmp_path / "partial")
    assert partial["primary_intent_to_diagnose"] == {
        "tp": 43, "fn": 1, "fp": 0, "tn": 7, "unresolved_negative": 1,
        "sensitivity": 43 / 44, "fpr": 0, "negative_failure_rate": 1 / 8}
    assert partial["ambiguity_rejected_fraction"] == 2 / 52
    assert partial["signed_bp_delta_mae"] is None
    assert partial["strata"]["monomer_length_bp"]


def test_invalid_score_or_tampered_truth_fails_closed(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _perfect_predictions(bundle, pred)
    rows = read_jsonl(pred)
    rows[0]["event_score"] = float("nan")
    write_jsonl(pred, rows)
    with pytest.raises(ValueError, match="finite event_score"):
        score(bundle, B2 / "metrics_frozen.json", pred, tmp_path / "bad")
    (bundle / "truth.jsonl").write_text((bundle / "truth.jsonl").read_text() + "\n")
    with pytest.raises(ValueError, match="truth.jsonl hash mismatch"):
        score(bundle, B2 / "metrics_frozen.json", pred, tmp_path / "tamper")
