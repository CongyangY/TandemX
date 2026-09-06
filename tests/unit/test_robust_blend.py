from benchmarks.abundance.select_robust_blend import select_robust_candidate


def _rows(seed_outcomes: dict[int, tuple[str, ...]], alpha: float, threshold: float) -> list[dict]:
    return [
        {
            "seed": str(seed), "outcome": outcome,
            "blend_alpha": str(alpha), "decision_threshold": str(threshold),
        }
        for seed, outcomes in seed_outcomes.items()
        for outcome in outcomes
    ]


def test_seed_robust_selection_rejects_pooled_gain_with_one_seed_fp() -> None:
    baseline = {seed: ("TP", "FN", "TN", "TN") for seed in (1, 2, 3)}
    moderate = {seed: ("TP", "TP", "TN", "TN") for seed in (1, 2, 3)}
    aggressive = {
        1: ("TP", "TP", "TN", "TN"),
        2: ("TP", "TP", "TN", "TN"),
        3: ("TP", "TP", "FP", "TN"),
    }
    candidates = {
        "blend_a0_t0.6": _rows(baseline, 0, .6),
        "blend_a0.5_t0.5": _rows(moderate, .5, .5),
        "blend_a0.75_t0.5": _rows(aggressive, .75, .5),
    }
    rule = {
        "baseline_candidate": "blend_a0_t0.6",
        "constraints": {
            "per_seed_false_positive_rate_max_delta": 0,
            "per_seed_precision_min_delta": 0,
        },
    }
    selected, summaries, _baseline = select_robust_candidate(candidates, [1, 2, 3], rule)
    assert selected == "blend_a0.5_t0.5"
    eligible = {row["candidate_id"]: row["eligible"] for row in summaries}
    assert eligible["blend_a0.5_t0.5"] is True
    assert eligible["blend_a0.75_t0.5"] is False
