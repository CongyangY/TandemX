import math

from benchmarks.abundance.evaluate_multik_collapse import (
    _frozen_model,
    calibration_rows,
    score_estimate,
    summarize,
)


def row(estimate, assembly_bp, truth_ratio, method="single_k21", fit_status="ok"):
    return score_estimate(
        family_id="f1", period=100, estimate=estimate, assembly_bp=assembly_bp,
        truth_ratio=truth_ratio, threshold=.6, method=method, fit_status=fit_status,
    )


def test_score_estimate_preserves_missing_and_native_categories() -> None:
    assert row(100, 50 * 100, .5)["outcome"] == "TP"
    assert row(100, 100 * 100, 1)["outcome"] == "TN"
    assert row(100, 0, 0)["native_status"] == "reads_only"
    assembly_only = row(0, 100, 1)
    assert assembly_only["native_status"] == "assembly_only"
    assert assembly_only["predicted_assembly_read_ratio"] == "NA"
    missing = row(None, 100, .5, method="multik_loglinear", fit_status="missing")
    assert missing["outcome"] == "NA" and missing["read_estimated_bp"] == "NA"


def test_score_estimate_rejects_invalid_numeric_inputs() -> None:
    for estimate in (-1, math.inf, math.nan):
        try:
            row(estimate, 100, .5)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid estimate accepted")


def test_summary_uses_only_available_denominators() -> None:
    rows = []
    for outcome in ("TP", "FN", "FP", "TN", "NA"):
        rows.append({"method": "m", "coverage": 1, "substitution_rate": 0,
                     "assembly_fraction": .5, "outcome": outcome})
    summary = summarize(rows)[0]
    assert summary["available"] == 4 and summary["unavailable"] == 1
    assert summary["sensitivity"] == .5
    assert summary["false_positive_rate"] == .5
    assert summary["precision"] == .5


def test_depth_rule_selection_preserves_seed_level_false_positive_constraint() -> None:
    pairs = []
    for seed in (1, 2, 3):
        for positive, single_ratio, multi_ratio in (
            (True, .58, .48), (True, .62, .52), (False, .65, .51), (False, .7, .61)
        ):
            pair = {}
            assembly_bp = single_ratio * 100 * 10
            for method, ratio in (("single_k21", single_ratio), ("multik_loglinear", multi_ratio)):
                estimate = 100.0 if method == "single_k21" else assembly_bp / ratio / 10
                scored = score_estimate(
                    family_id=f"f{len(pairs)}", period=10, estimate=estimate,
                    assembly_bp=assembly_bp,
                    truth_ratio=.5 if positive else 1.0, threshold=.6,
                    method=method, fit_status="ok",
                )
                scored.update(seed=seed, coverage=1, substitution_rate=0,
                              assembly_fraction=.5 if positive else 1,
                              truth_assembly_read_ratio=.5 if positive else 1,
                              truth_underrepresented=positive,
                              assembly_predicted_bp=assembly_bp,
                              period=10, estimated_haploid_depth=1.0)
                pair[method] = scored
            pairs.append(pair)
    selected, folds = calibration_rows(pairs, [1, 2, 3])
    assert selected == .5
    assert len(folds) == 4
    for row in folds:
        assert row["calibrated_FP"] <= row["baseline_FP"]
        assert row["calibrated_TP"] >= row["baseline_TP"]


def test_heldout_model_requires_exact_predeclared_semantics() -> None:
    config = {"seeds": {"development": [1, 2, 3]}, "collapse_model": {
        "method": "multik_loglinear_else_single_k21",
        "k_values": [15, 21, 27, 31],
        "low_depth_cutoff": 2.0,
        "low_depth_threshold": 0.5,
        "standard_depth_threshold": 0.6,
        "unavailable_rule": "fallback_single_k21",
        "calibration_seeds": [1, 2, 3],
        "calibration_validation_sha256": "a" * 64,
        "calibration_table_sha256": "b" * 64,
        "calibration_environment_sha256": "c" * 64,
    }}
    assert _frozen_model(config) == 0.5
    config["collapse_model"]["low_depth_threshold"] = 0.51
    try:
        _frozen_model(config)
    except ValueError:
        pass
    else:
        raise AssertionError("undeclared held-out threshold accepted")
