from __future__ import annotations

from benchmarks.scripts.archive_quantify_calibration import build_decision


def _metric(seed: int, family: str, method: str, error: float, estimate: float, control: float) -> dict[str, object]:
    return {
        "seed": seed,
        "condition_id": "c1",
        "family_id": family,
        "method": method,
        "absolute_relative_error": error,
        "signed_relative_error": -error,
        "estimator_minus_sampling_oracle": -error,
        "estimated_copy_number": estimate,
        "control_mean_depth": control,
    }


def test_build_decision_retains_tradeoffs_and_selects_only_a_candidate() -> None:
    methods = (
        "baseline_total_bases",
        "oracle_error_survival",
        "empirical_controls",
        "empirical_controls_plus_oracle_error",
    )
    metrics = []
    for seed in (1, 2):
        metrics.extend(
            [
                _metric(seed, "f1", "baseline_total_bases", 0.4, 6, 0),
                _metric(seed, "f1", "oracle_error_survival", 0.3, 7, 0),
                _metric(seed, "f1", "empirical_controls", 0.2, 8, 3),
                _metric(seed, "f1", "empirical_controls_plus_oracle_error", 0.2, 8, 3),
            ]
        )
    executions = [
        {"method": method, "runtime_seconds": index + 1, "peak_rss_mib": 10 + index}
        for index, method in enumerate(methods)
    ]
    decision = build_decision(metrics, executions)
    assert decision["control_plus_oracle_estimates_identical"] is True
    assert decision["posthoc_candidate"]["all_seed_improvements_positive"] is True
    assert decision["posthoc_candidate"]["mean_absolute_relative_error"] == 0.2
    assert decision["promotion"]["passed"] is False
