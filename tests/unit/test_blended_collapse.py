import pytest

from benchmarks.abundance.evaluate_blended_collapse import blend_row, select_candidate


def _row(method: str, estimate: float, outcome: str) -> dict:
    return {
        "seed": "1", "unit_substitution_rate": "0.03", "array_fragments": "1",
        "fragment_gap_bp": "500", "coverage": "5", "substitution_rate": "0.001",
        "assembly_fraction": "0.5", "family_id": "f1",
        "truth_assembly_read_ratio": "0.5", "truth_underrepresented": "True",
        "assembly_predicted_bp": "50", "period": "10", "estimated_haploid_depth": "5",
        "method": method, "read_estimated_copies": str(estimate),
        "read_estimated_bp": str(estimate * 10), "predicted_assembly_read_ratio": "0.5",
        "native_status": "possible_collapse", "predicted_underrepresented": "True",
        "outcome": outcome, "fit_status": "test", "decision_threshold": "0.6",
    }


def test_log_space_blend_and_alpha_zero_anchor() -> None:
    pair = {
        "single_k21": _row("single_k21", 10, "TP"),
        "multik_loglinear": _row("multik_loglinear", 40, "TP"),
    }
    anchor = blend_row(pair, 0, .6)
    blended = blend_row(pair, .5, .6)
    assert anchor["read_estimated_copies"] == 10
    assert blended["read_estimated_copies"] == pytest.approx(20)
    assert blended["candidate_id"] == "blend_a0.5_t0.6"


def test_selection_maximizes_sensitivity_under_baseline_constraints() -> None:
    baseline = [
        {"outcome": outcome} for outcome in ("TP", "FN", "TN", "TN")
    ]
    candidates = {
        "anchor": [{"outcome": outcome, "blend_alpha": 0, "decision_threshold": .6}
                   for outcome in ("TP", "FN", "TN", "TN")],
        "better": [{"outcome": outcome, "blend_alpha": .25, "decision_threshold": .55}
                   for outcome in ("TP", "TP", "TN", "TN")],
        "more_fp": [{"outcome": outcome, "blend_alpha": .5, "decision_threshold": .6}
                    for outcome in ("TP", "TP", "FP", "TN")],
    }
    rule = {"selection_constraints": {
        "false_positive_rate_max_delta_vs_single_k21": 0,
        "precision_min_delta_vs_single_k21": 0,
    }}
    candidate_id, metrics = select_candidate(candidates, baseline, rule)
    assert candidate_id == "better"
    assert metrics["sensitivity"] == 1
    assert metrics["false_positive_rate"] == 0
