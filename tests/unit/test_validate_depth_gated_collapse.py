import json
from pathlib import Path

import pytest

from benchmarks.abundance.run import validate_frozen_classifier
from benchmarks.abundance.validate_depth_gated_collapse import select_depth_gated_rows
from benchmarks.challenge.schema import digest_file


MODEL = {
    "low_depth_cutoff": 2.0,
    "low_depth_strategy": {"decision_threshold": 0.6},
    "standard_depth_strategy": {"blend_alpha": 0.5, "decision_threshold": 0.5},
}


def _row(seed: int, depth: float, method: str, estimate: float, outcome: str) -> dict:
    return {
        "seed": str(seed),
        "unit_substitution_rate": "0.03",
        "array_fragments": "1",
        "coverage": "1" if depth < 2 else "5",
        "substitution_rate": "0.001",
        "assembly_fraction": "0.5",
        "family_id": "f1",
        "fragment_gap_bp": "500",
        "truth_assembly_read_ratio": "0.5",
        "truth_underrepresented": "True",
        "assembly_predicted_bp": "50",
        "period": "10",
        "estimated_haploid_depth": str(depth),
        "method": method,
        "read_estimated_copies": str(estimate),
        "read_estimated_bp": str(estimate * 10),
        "predicted_assembly_read_ratio": str(50 / (estimate * 10)),
        "native_status": "possible_collapse" if outcome == "TP" else "represented",
        "predicted_underrepresented": str(outcome == "TP"),
        "outcome": outcome,
        "fit_status": "single_k_baseline" if method == "single_k21" else "ok",
        "decision_threshold": "0.6",
    }


def test_select_depth_gated_rows_uses_depth_only() -> None:
    rows = []
    for seed, depth in ((1, 1.1), (2, 5.1)):
        rows.extend((
            _row(seed, depth, "single_k21", 8, "FN"),
            _row(seed, depth, "multik_loglinear", 18, "TP"),
        ))
    baseline, selected, high = select_depth_gated_rows(rows, MODEL, [1, 2])
    assert [row["outcome"] for row in baseline] == ["FN", "FN"]
    assert selected[0]["outcome"] == "FN"
    assert selected[0]["fit_status"].startswith("depth_lt_2")
    assert selected[1]["outcome"] == "TP"
    assert selected[1]["fit_status"].startswith("depth_ge_2")
    assert all(row["candidate_id"] == "blend_a0.5_t0.5" for row in high)


def test_select_depth_gated_rows_rejects_incomplete_pairs() -> None:
    rows = [_row(1, 5.1, "single_k21", 8, "FN")]
    with pytest.raises(ValueError, match="Incomplete"):
        select_depth_gated_rows(rows, MODEL, [1])


def test_select_depth_gated_rows_rejects_unexpected_method() -> None:
    row = _row(1, 5.1, "single_k21", 8, "FN")
    row["method"] = "posthoc"
    with pytest.raises(ValueError, match="Unexpected method"):
        select_depth_gated_rows([row], MODEL, [1])


def _frozen_config(tmp_path: Path) -> tuple[dict, Path]:
    development = tmp_path / "development"
    development.mkdir()
    seeds = [1, 2, 3, 4, 5, 6]
    heldout = [7, 8, 9]
    strategy = {
        "method": "depth_gated_log_space_blend",
        "estimated_depth_source": "single_k21_estimated_haploid_depth",
        "low_depth_cutoff": 2.0,
        "low_depth_strategy": {
            "method": "single_k21", "blend_alpha": 0.0, "decision_threshold": 0.6,
        },
        "standard_depth_strategy": {
            "method": "log_space_single_multik_blend",
            "blend_alpha": 0.5,
            "decision_threshold": 0.5,
        },
        "multik_k_values": [15, 21, 27, 31],
        "unavailable_or_nonpositive_rule": "fallback_single_k21",
    }
    development_config = {
        "development_sources": [
            {"seeds": seeds[:3]},
            {"seeds": seeds[3:]},
        ],
        "reserved_future_heldout_seeds": heldout,
        "classifier_rule": strategy,
    }
    (development / "run_config.json").write_text(json.dumps(development_config))
    (development / "comparison_metrics.tsv").write_text("seed\n1\n")
    metrics = {
        "full_sensitivity_delta": 0.06,
        "full_false_positive_rate_delta": 0.0,
        "full_precision_delta": 0.01,
        "minimum_cohort_sensitivity_delta": 0.05,
        "maximum_cohort_false_positive_rate_delta": 0.0,
        "minimum_cohort_precision_delta": 0.0,
        "minimum_seed_sensitivity_delta": 0.03,
        "maximum_seed_false_positive_rate_delta": 0.0,
        "minimum_seed_precision_delta": 0.0,
    }
    (development / "selection.tsv").write_text(
        "candidate_id\tmethod\tlow_depth_cutoff\tlow_depth_method\t"
        "low_depth_blend_alpha\tlow_depth_decision_threshold\tstandard_method\t"
        "standard_blend_alpha\tstandard_decision_threshold\tpassed\n"
        "depth_gated_blend_v3\tdepth_gated_log_space_blend\t2.0\tsingle_k21\t"
        "0.0\t0.6\tlog_space_single_multik_blend\t0.5\t0.5\tTrue\n"
    )
    (development / "environment.json").write_text(json.dumps({
        "consumed_development_seeds": seeds,
        "reserved_future_heldout_seeds": heldout,
        "config_sha256": digest_file(development / "run_config.json"),
    }))
    (development / "validation.json").write_text(json.dumps({
        "complete": True,
        "development_only": True,
        "post_failed_heldout_refinement": True,
        "selected_candidate": "depth_gated_blend_v3",
        "consumed_development_seeds": seeds,
        "reserved_future_heldout_seeds": heldout,
        "comparison_metrics_sha256": digest_file(development / "comparison_metrics.tsv"),
        "selection_sha256": digest_file(development / "selection.tsv"),
        "acceptance": {**metrics, "passed": True},
    }))
    model = {
        "method": strategy["method"],
        "selection_method": "post_failed_heldout_coverage_diagnostic",
        "estimated_depth_source": strategy["estimated_depth_source"],
        "low_depth_cutoff": strategy["low_depth_cutoff"],
        "low_depth_strategy": strategy["low_depth_strategy"],
        "standard_depth_strategy": strategy["standard_depth_strategy"],
        "k_values": strategy["multik_k_values"],
        "fallback_rule": strategy["unavailable_or_nonpositive_rule"],
        "development_seeds": seeds,
        "development_result": str(development),
        "selection_gates": metrics,
        "observed_development_metrics": metrics,
    }
    for field, filename in {
        "development_validation_sha256": "validation.json",
        "development_selection_sha256": "selection.tsv",
        "development_metrics_sha256": "comparison_metrics.tsv",
        "development_environment_sha256": "environment.json",
        "development_config_sha256": "run_config.json",
    }.items():
        model[field] = digest_file(development / filename)
    config = {
        "seeds": {"development": seeds, "heldout": heldout},
        "classifier_development_rule": {
            "method": "log_space_single_multik_blend",
            "alpha_grid": [0, 0.25, 0.5, 0.75, 1],
            "decision_threshold_grid": [0.45, 0.5, 0.55, 0.6],
            "multik_k_values": [15, 21, 27, 31],
            "unavailable_or_nonpositive_rule": "fallback_single_k21",
        },
        "classifier_model": model,
    }
    return config, development


def test_depth_gated_freeze_guard_verifies_all_evidence(tmp_path: Path) -> None:
    config, development = _frozen_config(tmp_path)
    assert validate_frozen_classifier(config, "heldout") == config[
        "classifier_model"
    ]["observed_development_metrics"]
    (development / "selection.tsv").write_text("tampered\n")
    with pytest.raises(ValueError, match="evidence differs"):
        validate_frozen_classifier(config, "heldout")


def test_depth_gated_freeze_guard_rejects_changed_rule(tmp_path: Path) -> None:
    config, _ = _frozen_config(tmp_path)
    config["classifier_model"]["low_depth_cutoff"] = 1.5
    with pytest.raises(ValueError, match="disjoint frozen"):
        validate_frozen_classifier(config, "heldout")
