import pytest

from benchmarks.scripts.archive_classifier_validation import (
    confusion,
    validate_paired_classifier_rows,
)


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
