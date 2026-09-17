"""Frozen C3 equal-length technical score and denominator gates."""

import json
from pathlib import Path

import pytest

from benchmarks.m2_routes.native_equal_length_score import (
    classify_length,
    classify_paths,
    run,
    summarize,
)


def test_path_decision_abstains_on_unresolved_or_mixed_reads() -> None:
    assembly = {"state": "RESOLVED", "labels": ["C3+", "C3-"]}
    reads = {str(i): {"state": "RESOLVED", "labels": ["C3+", "C3+"]}
             for i in range(3)}
    assert classify_paths(assembly, reads)["state"] == "DISCORDANT"
    assert classify_paths({"state": "AMBIGUOUS", "reason": "budget"}, reads) == \
        {"state": "ABSTAIN", "reason": "assembly_budget"}
    reads["2"]["labels"] = ["C3-", "C3+"]
    assert classify_paths(assembly, reads)["reason"] == "mixed_native_read_paths"
    reads["2"]["state"] = "AMBIGUOUS"
    assert classify_paths(assembly, reads)["reason"] == \
        "insufficient_resolved_native_read_paths"


def test_same_read_span_baseline_cannot_detect_equal_length_edit() -> None:
    equal = classify_length(3560, [3555, 3560, 3560])
    assert equal["state"] == "SUPPORTED"
    assert equal["signed_bp_delta_assembly_minus_read"] == 0
    assert classify_length(2670, [3555, 3560, 3560])["state"] == "DISCORDANT"
    assert classify_length(3560, [3560, 3560])["state"] == "ABSTAIN"


def test_existing_output_preserved(tmp_path: Path) -> None:
    destination = tmp_path / "existing"
    destination.mkdir()
    (destination / "sentinel").write_text("retain")
    with pytest.raises(FileExistsError, match="preserve prior evidence"):
        run(destination)
    assert (destination / "sentinel").read_text() == "retain"


def test_archived_six_case_denominator_and_operation_outcomes() -> None:
    evidence = Path(__file__).resolve().parents[2] / \
        "benchmarks/m2_routes/native_equal_length_evidence_20260917"
    rows = [json.loads(line) for line in (evidence / "per_case.jsonl").read_text().splitlines()]
    expected = {
        "intact": ("SUPPORTED", "SUPPORTED", "NOT_IDENTIFIABLE_BY_OPERATIONAL_LABEL_ORIENTATION_PATH"),
        "invert_tile_05": ("DISCORDANT", "SUPPORTED", "IDENTIFIABLE"),
        "invert_block_05_09": ("DISCORDANT", "SUPPORTED", "IDENTIFIABLE"),
        "swap_adjacent_05_06": ("SUPPORTED", "SUPPORTED", "NOT_IDENTIFIABLE_BY_OPERATIONAL_LABEL_ORIENTATION_PATH"),
        "swap_distant_05_15": ("SUPPORTED", "SUPPORTED", "NOT_IDENTIFIABLE_BY_OPERATIONAL_LABEL_ORIENTATION_PATH"),
        "replace_tile_05_with_06": ("SUPPORTED", "SUPPORTED", "NOT_IDENTIFIABLE_BY_OPERATIONAL_LABEL_ORIENTATION_PATH"),
    }
    assert len(rows) == 6 and {row["operation"] for row in rows} == set(expected)
    for row in rows:
        m2, span, identifiability = expected[row["operation"]]
        assert (row["m2_technical"]["state"], row["length_baseline"]["state"],
                row["path_identifiability"]) == (m2, span, identifiability)
        assert row["biological_audit_state"] == "NOT_EVALUATED_PAIRING_UNVERIFIED"
    assert summarize(rows) == json.loads((evidence / "summary.json").read_text())
    with pytest.raises(ValueError, match="six full-denominator"):
        summarize(rows[:-1])
