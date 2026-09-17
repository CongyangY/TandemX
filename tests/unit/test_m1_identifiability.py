from __future__ import annotations

from benchmarks.m1_identifiability.diagnostics import (
    catalogue_diagnostics, read_assignment_diagnostics,
)


def test_exact_duplicate_refuses_individual_abundance() -> None:
    result = catalogue_diagnostics({"a": "ACGTAC", "b": "ACGTAC"}, 3)
    assert result["signature_rank"] == 1
    assert result["effective_identifiable_family_dimensions"] == 1
    assert all(item["status"] == "NON_IDENTIFIABLE" for item in result["family"].values())
    assert result["pairwise"][0]["identical_kmer_signature"] is True


def test_full_rank_without_unique_word_is_partial() -> None:
    result = catalogue_diagnostics({"a": "AAC", "b": "ACC"}, 1)
    assert result["signature_rank"] == 2
    assert all(item["status"] == "PARTIALLY_IDENTIFIABLE" for item in result["family"].values())


def test_distinct_exclusive_signatures_are_structurally_identifiable() -> None:
    result = catalogue_diagnostics({"a": "AAAAC", "b": "CCCCG"}, 2)
    assert result["signature_rank"] == 2
    assert all(item["status"] == "IDENTIFIABLE" for item in result["family"].values())


def test_known_source_ties_and_rejection_have_full_denominator() -> None:
    result = read_assignment_diagnostics(
        {"a": "AAAA", "b": "AAAA"}, ["AAAA", "CCCC"], ["a", "unknown"], 0
    )
    assert result["a"]["tied_best"] == 1
    assert result["a"]["correct_unique"] == 0
    assert result["a"]["mean_best_hit_tie_entropy_bits_per_accepted"] == 1.0
    assert result["unknown"]["rejected"] == 1
    assert result["unknown"]["total"] == 1
