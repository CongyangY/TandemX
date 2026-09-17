"""Regression checks for the isolated M2 transition-graph research route."""

import pytest

from benchmarks.m2_routes.graph import ReadPath, audit_transition_graph, to_common_prediction
from benchmarks.m2_routes.graph.prototype import transition_counts


MONOMERS = {"A": "ACGTTGCA", "B": "GGAATCCA", "C": "CATGAGTC",
            "A_snp": "ACGTTGCG"}


def reads(labels: str | tuple[str, ...], n: int = 3, *, confidence: float = 1.0,
          haplotype: str | None = None) -> list[ReadPath]:
    return [ReadPath(f"mol{i}", tuple(labels), (confidence,) * len(labels),
                     f"platform{i % 2}", haplotype=haplotype) for i in range(n)]


def test_transition_graph_retains_flank_phase_and_copy_multiplicity() -> None:
    assert transition_counts(tuple("ABAB")) != transition_counts(tuple("BABA"))
    assert transition_counts(tuple("ABAB")) != transition_counts(tuple("ABABAB"))


def test_copy_compression_candidate_requires_cross_profile_support() -> None:
    result = audit_transition_graph(tuple("ABC"), reads("ABCABC"),
                                    monomer_sequences=MONOMERS)
    assert result.status == "DISCORDANT"
    assert result.event_class == "copy_count_or_graph_multiplicity"
    assert result.discordant_count == 3
    assert result.graph_distance and result.graph_distance > 0


def test_unedited_array_supported_only_at_graph_level() -> None:
    result = audit_transition_graph(tuple("ABCABC"), reads("ABCABC"),
                                    monomer_sequences=MONOMERS)
    assert result.status == "SUPPORTED"
    assert result.graph_distance == 0
    assert "not_complete_order" in result.warning


def test_same_edge_spectrum_but_different_order_abstains() -> None:
    assembly = tuple("ABACA")
    observed = tuple("ACABA")
    assert transition_counts(assembly) == transition_counts(observed)
    result = audit_transition_graph(assembly, reads("ACABA"),
                                    monomer_sequences=MONOMERS)
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "order_unidentifiable_from_transitions"


def test_close_monomer_variant_and_systematic_motif_error_abstain() -> None:
    result = audit_transition_graph(tuple("ABA"), reads(("A", "A_snp", "A")),
                                    monomer_sequences=MONOMERS)
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "label_identity_unresolved"


def test_unverified_catalogue_abstains_on_candidate() -> None:
    result = audit_transition_graph(tuple("ABC"), reads("ABCABC"))
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "catalogue_unverified"


def test_low_support_and_low_confidence_error_paths_abstain() -> None:
    one = audit_transition_graph(tuple("ABC"), reads("ABCABC", 1),
                                 monomer_sequences=MONOMERS)
    assert one.status == "INSUFFICIENT_READ_SUPPORT"
    noisy = audit_transition_graph(tuple("ABC"), reads("ABCABC", confidence=0.6),
                                   monomer_sequences=MONOMERS)
    assert noisy.status == "INSUFFICIENT_READ_SUPPORT"
    assert noisy.rejected_count == 3


def test_single_correlated_error_profile_abstains() -> None:
    correlated = [ReadPath(f"mol{i}", tuple("ABCABC"), (1.0,) * 6, "same_basecaller")
                  for i in range(3)]
    result = audit_transition_graph(tuple("ABC"), correlated,
                                    monomer_sequences=MONOMERS)
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "profile_limited"


def test_mixed_haplotypes_abstain_even_when_majority_discordant() -> None:
    observed = reads("ABCABC", 3, haplotype="h1")
    observed += [ReadPath("m4", tuple("ABC"), (1.0,) * 3,
                          "platform0", haplotype="h2")]
    result = audit_transition_graph(tuple("ABC"), observed,
                                    monomer_sequences=MONOMERS)
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "mixed_haplotypes"


def test_unphased_two_molecule_minority_abstains() -> None:
    observed = reads("ABCABC", 3)
    observed += [ReadPath(f"minor{i}", tuple("ABC"), (1.0,) * 3,
                          f"platform{i}") for i in range(2)]
    result = audit_transition_graph(tuple("ABC"), observed,
                                    monomer_sequences=MONOMERS)
    assert result.status == "AMBIGUOUS"
    assert result.event_class == "mixed_transition_graphs"


def test_no_full_flank_support_abstains() -> None:
    observed = [ReadPath(f"mol{i}", tuple("ABCABC"), (1.0,) * 6,
                         f"platform{i % 2}", both_flanks_verified=False)
                for i in range(3)]
    result = audit_transition_graph(tuple("ABC"), observed,
                                    monomer_sequences=MONOMERS)
    assert result.status == "INSUFFICIENT_READ_SUPPORT"


def test_duplicate_molecule_and_missing_assignments_fail_closed() -> None:
    with pytest.raises(ValueError, match="unique molecule"):
        audit_transition_graph(tuple("ABC"), reads("ABC", 1) * 3,
                               monomer_sequences=MONOMERS)
    with pytest.raises(ValueError, match="Each monomer"):
        audit_transition_graph(tuple("ABC"),
                               [ReadPath("mol", tuple("ABC"), (1.0,), "p")],
                               monomer_sequences=MONOMERS)


def test_common_prediction_keeps_abstention_and_uncalibrated_binary_score() -> None:
    call = audit_transition_graph(tuple("ABC"), reads("ABCABC"),
                                  monomer_sequences=MONOMERS)
    prediction = to_common_prediction("opaque-1", call)
    assert prediction["status"] == "ok"
    assert prediction["event_score"] == 1.0
    assert prediction["predicted_signed_bp_delta"] is None
    abstain = audit_transition_graph(tuple("ABC"), reads("ABCABC", 1),
                                     monomer_sequences=MONOMERS)
    prediction = to_common_prediction("opaque-2", abstain)
    assert prediction["status"] == "abstain"
    assert prediction["event_score"] is None
