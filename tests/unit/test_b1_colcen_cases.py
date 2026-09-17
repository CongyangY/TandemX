"""Source, edit-ledger and abstention tests for the Col-CEN development bundle."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.scripts.build_b1_colcen_cases import build, load_source
from benchmarks.scripts.build_b1_structure_cases import reverse_complement, sha256_file, write_jsonl
from benchmarks.scripts.score_b1_colcen_cases import score
from benchmarks.scripts.score_b1_structure_cases import read_jsonl


ROOT = Path(__file__).resolve().parents[2]
V3 = ROOT / "benchmarks/controlled_collapse/v3"


def _build(tmp_path: Path) -> Path:
    out = tmp_path / "bundle"
    build(V3 / "protocol.json", V3 / "source_arrays.fa", V3 / "source_catalogue.json", out)
    return out


def _predictions(bundle: Path, path: Path) -> None:
    # Deliberately derives predictions from truth solely to test the scorer.
    write_jsonl(path, [dict(case_id=t["case_id"], status="ok",
                            event_score=1.0 if t["event_status"] == "positive" else 0.0,
                            event_type=t["event_type"],
                            predicted_signed_bp_delta=t["signed_bp_delta"],
                            predicted_breakpoints_bp=sorted({s["edited_start0"] for s in t["coordinate_chain"]
                                if s["operation"] == "delete"} |
                                {p for s in t["coordinate_chain"] if s["operation"] in
                                 {"duplicate", "reorder", "reverse_complement"}
                                 for p in (s["edited_start0"], s["edited_end0"])}))
                       for t in read_jsonl(bundle / "truth.jsonl")])


def test_real_source_catalogue_and_bundle_replay(tmp_path: Path) -> None:
    protocol, arrays, monomers = load_source(V3 / "protocol.json", V3 / "source_arrays.fa",
                                             V3 / "source_catalogue.json")
    assert len(arrays) == len(monomers) == 3
    assert {len(sequence) for sequence in arrays.values()} == {3560}
    assert {len(monomer) for monomer in monomers.values()} == {178}
    assert protocol["source_assembly_sha256"] == json.loads((V3 / "source_catalogue.json").read_text())["source_assembly_sha256"]
    bundle = _build(tmp_path)
    receipt = json.loads((bundle / "receipt.json").read_text())
    frozen = json.loads((V3 / "development_bundle/receipt.json").read_text())
    assert receipt["inputs_sha256"] == frozen["inputs_sha256"] == sha256_file(bundle / "inputs.jsonl")
    assert receipt["truth_sha256"] == frozen["truth_sha256"] == sha256_file(bundle / "truth.jsonl")
    assert receipt["case_count"] == 39


def test_exact_native_interval_edits_and_truth_scope(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    inputs = {row["case_id"]: row for row in read_jsonl(bundle / "inputs.jsonl")}
    truths = read_jsonl(bundle / "truth.jsonl")
    assert len(inputs) == len(truths) == 39
    assert len({row["assembly_sequence"] for row in inputs.values()}) == 36
    assert {row["array_id"] for row in truths} == {"C1", "C2", "C3"}
    assert sum(row["event_status"] == "positive" for row in truths) == 33
    for t in truths:
        public = inputs[t["case_id"]]
        assert not set(public).intersection({"event_status", "event_type", "signed_bp_delta", "edited_tile_id_path"})
        source = next(iter(public["raw_read_sequences"].values()))
        assert len(set(public["raw_read_sequences"].values())) == 1  # no independent molecules
        assert public["read_pairing_status"] == "synthetic_simulated"
        assert t["truth_scope"] == "injected_edit_only"
        assert t["unit_truth_status"] == "operational_178bp_tiles_not_native_monomer_truth"
        assert t["hor_period_status"] == "unverified"
        assert t["signed_bp_delta"] == len(public["assembly_sequence"]) - len(source)
        rebuilt = [None] * len(public["assembly_sequence"])
        for segment in t["coordinate_chain"]:
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
        assert "".join(rebuilt) == public["assembly_sequence"]
    complete = [t for t in truths if t["event_type"] == "complete_deletion"]
    assert len(complete) == 3 and all(t["signed_bp_delta"] == -3560 for t in complete)
    compression = [t for t in truths if t["event_type"] == "compression"]
    assert len(compression) == 3 and all(t["signed_bp_delta"] == -40 for t in compression)


def test_intent_to_diagnose_and_na_gates(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _predictions(bundle, pred)
    rows = read_jsonl(pred)
    truth = {t["case_id"]: t for t in read_jsonl(bundle / "truth.jsonl")}
    positive = next(row for row in rows if truth[row["case_id"]]["event_status"] == "positive")
    intact = next(row for row in rows if truth[row["case_id"]]["event_status"] == "intact_negative")
    for row in (positive, intact):
        row["status"] = "abstain"
        row["event_score"] = None
        row["predicted_signed_bp_delta"] = None
    write_jsonl(pred, rows)
    summary = score(bundle, V3 / "metrics_frozen.json", pred, tmp_path / "score")
    assert summary["denominator"] == 39
    assert summary["primary_intent_to_diagnose"] == {
        "tp": 32, "fn": 1, "fp": 0, "tn": 5, "unresolved_negative": 1,
        "sensitivity": 32 / 33, "fpr": 0, "negative_failure_rate": 1 / 6}
    assert summary["eligible_only"]["denominator"] == 37
    assert summary["ambiguity_rejected_fraction"] == 2 / 39
    assert summary["signed_bp_delta_mae"] is None
    assert summary["hor_period_error_bp"] is None
    assert summary["repeat_order_accuracy"] is None
    assert summary["runtime_seconds"] is None and summary["peak_rss_bytes"] is None
    assert summary["biological_accuracy_status"].startswith("blocked_")


def test_tautological_full_predictions_only_check_scorer(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _predictions(bundle, pred)
    summary = score(bundle, V3 / "metrics_frozen.json", pred, tmp_path / "score")
    assert summary["primary_intent_to_diagnose"]["tp"] == 33
    assert summary["primary_intent_to_diagnose"]["tn"] == 6
    assert summary["signed_bp_delta_mae"] == 0
    assert summary["breakpoint_localization_error_bp"] == 0


def test_invalid_score_and_source_tamper_fail_closed(tmp_path: Path) -> None:
    bundle = _build(tmp_path)
    pred = tmp_path / "pred.jsonl"
    _predictions(bundle, pred)
    rows = read_jsonl(pred)
    rows[0]["event_score"] = float("nan")
    write_jsonl(pred, rows)
    with pytest.raises(ValueError, match="finite event_score"):
        score(bundle, V3 / "metrics_frozen.json", pred, tmp_path / "bad_score")
    source = tmp_path / "bad.fa"
    source.write_text((V3 / "source_arrays.fa").read_text() + "\n")
    with pytest.raises(ValueError, match="source catalogue/FASTA hash"):
        load_source(V3 / "protocol.json", source, V3 / "source_catalogue.json")
