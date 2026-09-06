from benchmarks.abundance.validate_blended_collapse import (
    comparison_output_rows,
    heldout_acceptance,
)


def _rows(seed_outcomes: dict[int, tuple[str, ...]]) -> list[dict]:
    return [
        {"seed": str(seed), "outcome": outcome}
        for seed, outcomes in seed_outcomes.items()
        for outcome in outcomes
    ]


def test_heldout_gate_requires_improvement_in_every_seed() -> None:
    baseline = _rows({seed: ("TP", "FN", "FP", "TN", "TN") for seed in (1, 2, 3)})
    selected = _rows({
        1: ("TP", "TP", "TN", "TN", "TN"),
        2: ("TP", "TP", "FP", "TN", "TN"),
        3: ("TP", "FN", "FP", "TN", "TN"),
    })
    gates = {
        "full_sensitivity_min_delta_vs_single_k21": .1,
        "full_false_positive_rate_max_delta_vs_single_k21": 0,
        "full_precision_min_delta_vs_single_k21": 0,
        "minimum_seed_sensitivity_min_delta_vs_single_k21": .1,
        "maximum_seed_false_positive_rate_max_delta_vs_single_k21": 0,
        "minimum_seed_precision_min_delta_vs_single_k21": 0,
    }
    acceptance, by_seed = heldout_acceptance(baseline, selected, [1, 2, 3], gates)
    assert acceptance["full_sensitivity_delta"] > 0
    assert acceptance["minimum_seed_sensitivity_delta"] == 0
    assert acceptance["passed"] is False
    assert by_seed[3]["sensitivity_delta"] == 0


def test_comparison_output_normalizes_baseline_and_selected_fields() -> None:
    baseline = [{"method": "single_k21", "outcome": "TN"}]
    selected = [{
        "method": "blend_a0.5_t0.5", "outcome": "TN",
        "candidate_id": "blend_a0.5_t0.5", "blend_alpha": .5,
    }]
    combined = comparison_output_rows(baseline, selected)
    assert set(combined[0]) == set(combined[1])
    assert combined[0]["candidate_id"] == "single_k21"
