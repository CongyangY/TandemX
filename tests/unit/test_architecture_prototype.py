from __future__ import annotations

import pytest

from benchmarks.architecture.prototype import compare_anchored_reads, infer_architecture


MONOMERS = {"A": "ACGTTGCA", "B": "GGAATCCA", "C": "CATGAGTC"}


def _array(labels: str) -> str:
    return "".join(MONOMERS[label] for label in labels)


def test_exact_order_and_repeated_unit() -> None:
    result = infer_architecture(_array("ABABAB"), MONOMERS)
    assert tuple(copy.label for copy in result.copies) == tuple("ABABAB")
    assert result.cyclic_unit == ("A", "B")
    assert result.complete_unit_count == 3
    assert result.variant_copy_indices == ()
    assert result.status == "candidate_periodic"


def test_variant_order_is_reported_without_calling_it_validated_hor() -> None:
    result = infer_architecture(_array("ABAC"), MONOMERS)
    assert result.cyclic_unit == ("A", "B")
    assert result.variant_copy_indices == (3,)
    assert result.template_ambiguous is True
    assert result.status == "candidate_periodic_template_ambiguous"


@pytest.mark.parametrize("labels,variant_index", [("CBABABAB", 0), ("ABABABAC", 7)])
def test_first_or_last_variant_does_not_double_hor_unit(labels: str,
                                                          variant_index: int) -> None:
    result = infer_architecture(_array(labels), MONOMERS)
    assert result.cyclic_unit == ("A", "B")
    assert result.variant_copy_indices == (variant_index,)
    assert result.complete_unit_count == 4
    assert result.template_ambiguous is False


def test_equal_length_rearrangement_detected_beyond_length_only() -> None:
    # Equal bp means a plain read-vs-assembly length check is negative.
    assembly = _array("ABAC")
    reads = {"molecule-1": _array("ABAB"), "molecule-2": _array("ABAB")}
    assert len(assembly) == len(reads["molecule-1"])
    result = compare_anchored_reads(assembly, reads, MONOMERS)
    assert result.status == "candidate_discordance"
    assert result.discordant_reads == ("molecule-1", "molecule-2")


def test_copy_deletion_requires_two_independent_anchored_reads() -> None:
    assembly = _array("ABAB")
    one = compare_anchored_reads(assembly, {"r1": _array("ABABAB")}, MONOMERS)
    assert one.status == "unresolved"
    two = compare_anchored_reads(
        assembly, {"r1": _array("ABABAB"), "r2": _array("ABABAB")}, MONOMERS
    )
    assert two.status == "candidate_discordance"


def test_reverse_complement_equivalent_catalogue_is_rejected() -> None:
    motif = MONOMERS["A"]
    with pytest.raises(ValueError, match="not identifiable"):
        infer_architecture(motif * 2, {"A": motif, "A2": motif[2:] + motif[:2]})


def test_non_decomposable_array_is_unresolved() -> None:
    result = infer_architecture("TTTTTTTT", MONOMERS)
    assert result.status == "unresolved_no_full_decomposition"


def test_single_substitution_and_insertion_have_explicit_edit_cost() -> None:
    sequence = _array("ABAB")
    substituted = sequence[:2] + "T" + sequence[3:]
    inserted = sequence[:10] + "C" + sequence[10:]
    sub = infer_architecture(substituted, MONOMERS)
    ins = infer_architecture(inserted, MONOMERS)
    assert tuple(copy.label for copy in sub.copies) == tuple("ABAB")
    assert tuple(copy.label for copy in ins.copies) == tuple("ABAB")
    assert sum(copy.edit_distance for copy in sub.copies) == 1
    assert sum(copy.edit_distance for copy in ins.copies) == 1


def test_mixed_molecule_evidence_does_not_force_discordance() -> None:
    result = compare_anchored_reads(
        _array("ABAB"),
        {"r1": _array("ABABAB"), "r2": _array("ABABAB"),
         "r3": _array("ABAB"), "r4": _array("ABAB")},
        MONOMERS,
    )
    assert result.status == "unresolved_mixed_molecules"


def test_large_exhaustive_search_fails_before_consuming_unbounded_cpu() -> None:
    motif = "ACGTTGCAGGAATCCAGACTTGACCATAGTCA"
    with pytest.raises(ValueError, match="work estimate"):
        infer_architecture(motif * 60, {"A": motif})


def test_cyclic_unit_is_normalized_but_anchored_phase_shift_is_detected() -> None:
    result = compare_anchored_reads(
        _array("ABAB"), {"r1": _array("BABA"), "r2": _array("BABA")}, MONOMERS
    )
    assert result.assembly.cyclic_unit == ("A", "B")
    assert result.status == "candidate_discordance"


def test_variant_position_relative_to_left_flank_is_retained() -> None:
    result = compare_anchored_reads(
        _array("ABAC"), {"r1": _array("ACAB"), "r2": _array("ACAB")}, MONOMERS
    )
    assert result.status == "candidate_discordance"


@pytest.mark.parametrize("assembly", ["AAAAC" * 4, "AAACA" * 4])
def test_ambiguous_assembly_decomposition_cannot_be_scored(assembly: str) -> None:
    monomers = {"A": "AAAAA", "B": "AACCA"}
    result = compare_anchored_reads(
        assembly, {"r1": "AAAAA" * 4, "r2": "AAAAA" * 4}, monomers
    )
    assert result.assembly.status == "unresolved_ambiguous_decomposition"
    assert result.status == "unresolved"
