import pytest

from benchmarks.abundance.evaluate_depth_gated_classifier import (
    apply_depth_gate,
    evaluate_acceptance,
)


def _row(seed: int, depth: float, outcome: str, threshold: float) -> dict:
    return {
        "seed": str(seed), "unit_substitution_rate": ".03", "array_fragments": "1",
        "coverage": "1" if depth < 2 else "5", "substitution_rate": ".001",
        "assembly_fraction": ".5", "family_id": "f1", "estimated_haploid_depth": depth,
        "candidate_id": "single_k21" if threshold == .6 else "blend_a0.5_t0.5",
        "method": "single_k21" if threshold == .6 else "blend_a0.5_t0.5",
        "outcome": outcome, "decision_threshold": threshold, "blend_alpha": 0 if threshold == .6 else .5,
    }


RULE = {
    "low_depth_cutoff": 2.0,
    "low_depth_strategy": {"decision_threshold": .6},
    "standard_depth_strategy": {"decision_threshold": .5},
}


def test_depth_gate_retains_conservative_low_depth_and_high_depth_blend() -> None:
    low_single = _row(1, 1.1, "TN", .6)
    low_blend = _row(1, 1.1, "FP", .5)
    low = apply_depth_gate(low_single, low_blend, RULE)
    assert low["outcome"] == "TN" and low["decision_threshold"] == .6
    high_single = _row(1, 5.1, "FN", .6)
    high_blend = _row(1, 5.1, "TP", .5)
    high = apply_depth_gate(high_single, high_blend, RULE)
    assert high["outcome"] == "TP" and high["decision_threshold"] == .5


def test_depth_gate_rejects_unpaired_rows() -> None:
    single = _row(1, 1.1, "TN", .6)
    blend = _row(1, 1.1, "FP", .5)
    blend["family_id"] = "f2"
    with pytest.raises(ValueError, match="not paired"):
        apply_depth_gate(single, blend, RULE)


def test_depth_gated_acceptance_keeps_seed_guardrails() -> None:
    baseline = []
    selected = []
    for seed in (1, 2):
        for index, (old, new) in enumerate(
            zip(("TP", "FN", "FP", "TN"), ("TP", "TP", "FP", "TN")), 1
        ):
            base = _row(seed, 5, old, .6)
            base["family_id"] = f"f{index}"
            chosen = dict(base, outcome=new, method="depth_gated_blend_v3")
            baseline.append(base)
            selected.append(chosen)
    gates = {
        "full_sensitivity_min_delta_vs_single_k21": .1,
        "full_false_positive_rate_max_delta_vs_single_k21": 0,
        "full_precision_min_delta_vs_single_k21": 0,
        "minimum_cohort_sensitivity_min_delta_vs_single_k21": .1,
        "maximum_cohort_false_positive_rate_max_delta_vs_single_k21": 0,
        "minimum_cohort_precision_min_delta_vs_single_k21": 0,
        "minimum_seed_sensitivity_min_delta_vs_single_k21": .1,
        "maximum_seed_false_positive_rate_max_delta_vs_single_k21": 0,
        "minimum_seed_precision_min_delta_vs_single_k21": 0,
    }
    acceptance, cohorts, seeds = evaluate_acceptance(
        baseline, selected, {1: "a", 2: "b"}, gates
    )
    assert acceptance["passed"] is True
    assert set(cohorts) == {"a", "b"} and set(seeds) == {1, 2}
