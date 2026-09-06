import pytest

from benchmarks.scripts.archive_classifier_validation import (
    confusion,
    validate_paired_classifier_rows,
)
from benchmarks.scripts.archive_depth_gated_validation import validate_depth_gated_counts


def _rows(method: str, outcomes: tuple[str, ...]) -> list[dict]:
    return [
        {
            "seed": "1", "unit_substitution_rate": ".03", "array_fragments": "1",
            "coverage": "5", "substitution_rate": ".001", "assembly_fraction": ".5",
            "family_id": f"f{index}", "method": method, "outcome": outcome,
        }
        for index, outcome in enumerate(outcomes, 1)
    ]


def test_classifier_archive_requires_exact_paired_conditions() -> None:
    rows = _rows("single_k21", ("TP", "FN", "TN")) + _rows(
        "blend_a0.5_t0.5", ("TP", "TP", "TN")
    )
    grouped = validate_paired_classifier_rows(
        rows, "single_k21", "blend_a0.5_t0.5", 3
    )
    assert confusion(grouped["single_k21"])["sensitivity"] == 0.5
    assert confusion(grouped["blend_a0.5_t0.5"])["sensitivity"] == 1.0

    broken = rows[:-1]
    with pytest.raises(ValueError, match="incomplete, duplicated or unpaired"):
        validate_paired_classifier_rows(
            broken, "single_k21", "blend_a0.5_t0.5", 3
        )


def test_classifier_archive_rejects_unknown_outcome() -> None:
    with pytest.raises(ValueError, match="Unknown classifier outcome"):
        confusion([{"outcome": "MAYBE"}])


def test_depth_gated_archive_recomputes_strata_and_fallbacks() -> None:
    selected = [
        {**row, "estimated_haploid_depth": depth, "read_estimated_copies": "10"}
        for row, depth in zip(
            _rows("single_k21", ("TN", "TP"))
            + _rows("depth_gated_blend_v3", ("TN", "TP")),
            (1.1, 5.1, 1.1, 5.1),
        )
    ]
    raw = [
        {**row, "estimated_haploid_depth": depth, "read_estimated_copies": estimate}
        for row, depth, estimate in zip(
            _rows("single_k21", ("TN", "TP"))
            + _rows("multik_loglinear", ("NA", "NA")),
            (1.1, 5.1, 1.1, 5.1),
            ("10", "10", "NA", "NA"),
        )
    ]
    receipt = {
        "low_depth_rows": 1,
        "standard_depth_rows": 1,
        "standard_depth_fallback_rows": 1,
    }
    validate_depth_gated_counts(selected, raw, receipt, 2.0)
    receipt["standard_depth_fallback_rows"] = 0
    with pytest.raises(ValueError, match="stratum or fallback"):
        validate_depth_gated_counts(selected, raw, receipt, 2.0)
