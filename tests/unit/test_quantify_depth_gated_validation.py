from __future__ import annotations

import json
from pathlib import Path

import yaml

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.evaluate_quantify_depth_gated_validation import (
    CANDIDATE,
    METHODS,
    build_candidate,
    evaluate_gates,
)


def _row(seed: int, coverage: int, family: str, method: str, error: float, control: float):
    return {
        "seed": seed,
        "condition_id": f"c{coverage}",
        "coverage": coverage,
        "error_model": "clean",
        "family_id": family,
        "method": method,
        "absolute_relative_error": error,
        "signed_relative_error": -error,
        "control_mean_depth": control,
    }


def test_depth_gate_is_condition_level_and_passes_predeclared_summary_gates() -> None:
    raw = []
    for seed in (1, 2):
        for coverage, baseline, controls, depth in ((1, 0.4, 0.8, 1.0), (5, 0.4, 0.1, 3.0)):
            for family in ("f1", "f2"):
                raw.extend(
                    [
                        _row(seed, coverage, family, METHODS[0], baseline, 0),
                        _row(seed, coverage, family, METHODS[1], controls, depth),
                    ]
                )
    candidate, paired = build_candidate(raw, 2.0)
    assert {row["method"] for row in candidate} == {CANDIDATE}
    assert {
        (row["coverage"], row["candidate_selected_method"]) for row in candidate
    } == {(1, METHODS[0]), (5, METHODS[1])}
    metrics = [*raw, *candidate]
    executions = [{"status": "ok"} for _ in range(8)]
    limits = {
        "failed_executions_max": 0,
        "aggregate_mare_reduction_min": 0.1,
        "candidate_mare_max": 0.3,
        "minimum_seed_mare_reduction_strictly_greater_than": 0,
        "minimum_coverage_mare_reduction_min": 0,
        "paired_improved_fraction_min": 0.4,
        "paired_nonworse_fraction_min": 1,
        "selected_condition_count_per_branch_min": 1,
    }
    result = evaluate_gates(metrics, paired, executions, limits)
    assert result["status"] == "passed"
    assert result["paired_improved"] == 4
    assert result["paired_equal"] == 4
    assert result["selected_condition_counts"] == {METHODS[0]: 2, METHODS[1]: 2}


def test_depth_gate_retains_a_failed_seed_gate() -> None:
    raw = [
        _row(1, 5, "f1", METHODS[0], 0.1, 0),
        _row(1, 5, "f1", METHODS[1], 0.2, 3),
    ]
    candidate, paired = build_candidate(raw, 2.0)
    limits = {
        "failed_executions_max": 0,
        "aggregate_mare_reduction_min": -1,
        "candidate_mare_max": 1,
        "minimum_seed_mare_reduction_strictly_greater_than": 0,
        "minimum_coverage_mare_reduction_min": -1,
        "paired_improved_fraction_min": 0,
        "paired_nonworse_fraction_min": 0,
        "selected_condition_count_per_branch_min": 0,
    }
    result = evaluate_gates([*raw, *candidate], paired, [{"status": "ok"}], limits)
    assert result["status"] == "failed"
    assert "minimum_seed_mare_reduction" in result["failed_gate_names"]


def test_frozen_validation_design_matches_development_except_for_seeds() -> None:
    root = Path(__file__).resolve().parents[2]
    development = json.loads(
        (root / "benchmarks/configs/factorial_scale_v1.json").read_text()
    )
    validation = json.loads(
        (
            root
            / "benchmarks/configs/factorial_scale_quantify_validation_v1.json"
        ).read_text()
    )
    assert {key: value for key, value in development.items() if key != "seeds"} == {
        key: value for key, value in validation.items() if key != "seeds"
    }
    assert validation["seeds"]["validation"] == [6401, 6402, 6403]
    assert not set(validation["seeds"]["validation"]) & set(
        development["seeds"]["development"]
    )
    gates = yaml.safe_load(
        (root / "benchmarks/configs/quantify_depth_gated_validation_v1.yaml").read_text()
    )
    assert gates["seeds"] == validation["seeds"]["validation"]
    assert gates["candidate"]["control_mean_depth_threshold"] == 2.0
    assert gates["acceptance_gates"]["aggregate_mare_reduction_min"] == 0.02
    assert gates["dataset_generation"]["config_sha256"] == digest_file(
        root / "benchmarks/configs/factorial_scale_quantify_validation_v1.json"
    )
    assert gates["dataset_generation"]["length_histogram_sha256"] == digest_file(
        root / "paper/evidence/Mo17_input_qc/qc/length_histogram.tsv"
    )
    for relative, expected in gates["dataset_generation"]["source_sha256"].items():
        assert digest_file(root / relative) == expected
