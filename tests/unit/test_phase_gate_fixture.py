"""Regression tests for the development-only phase-gate fixture, not holdout data."""

import pytest

from benchmarks.challenge.phase_gate_fixture import CONDITIONS, generate_case


def test_development_case_is_reproducible() -> None:
    first = generate_case("dev-01", "clean", 2)
    second = generate_case("dev-01", "clean", 2)
    assert first == second
    assert first["read_bases"] >= 2 * first["genome_size"]
    assert first["source_hash"] == first["source_hash"]
    assert first["reads_hash"] == second["reads_hash"]


def test_truth_intervals_are_in_read_coordinates_and_counts_match() -> None:
    case = generate_case("dev-02", "clean", 5)
    lengths = {read_id: len(sequence) for read_id, sequence in case["reads"]}
    counted = {"F1": 0, "F2": 0}
    for read_id, start, end, family in case["truth_intervals"]:
        assert 0 <= start < end <= lengths[read_id]
        counted[family] += end - start
    assert counted == case["read_truth_bp"]


def test_background_homology_has_no_complete_family_truth() -> None:
    case = generate_case("dev-03", "background_homology", 2)
    # The added F1[:80]/random-80 background is deliberately unmasked.
    assert case["source_truth_bp"] == {"F1": 720, "F2": 1440}
    assert all(end - start <= 720 for start, end, _ in case["source_truth_intervals"]["F1"])


def test_interruption_bases_are_excluded_from_repeat_truth() -> None:
    case = generate_case("dev-04", "interruptions", 2)
    assert case["source_truth_bp"] == {"F1": 720, "F2": 1440}
    assert case["genome_size"] == 22_000
    assert case["source_truth_bp"]["F1"] < 6 * 120 + 300


def test_all_declared_conditions_generate_without_holdout() -> None:
    for condition in CONDITIONS:
        case = generate_case("dev-05", condition, 2)
        assert case["condition"] == condition
        assert case["coverage"] == 2
    with pytest.raises(ValueError, match="only dev seed"):
        generate_case("futureholdout-01", "clean", 2)
    with pytest.raises(ValueError, match="unknown condition"):
        generate_case("dev-05", "unknown", 2)


def test_catalogue_and_mask_truth_are_present_for_mutation_cases() -> None:
    for condition in ("substitutions", "indels_small", "indels_long", "partial", "heterogeneous"):
        case = generate_case("dev-06", condition, 2)
        assert set(case["catalogue"]) == {"F1", "F2"}
        assert set(case["source_truth_intervals"]) == {"F1", "F2"}
        assert case["source_truth_bp"]["F1"] > 0
        assert case["source_truth_bp"]["F2"] > 0
