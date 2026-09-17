from __future__ import annotations

from benchmarks.m2_routes.alignment.adapter import predict_case
from benchmarks.m2_routes.alignment.prototype import audit, decompose


M = {"A": "ACGTTGCA", "B": "GGAATCCA", "C": "CATGAGTC"}


def array(labels: str) -> str:
    return "".join(M[label] for label in labels)


def reads(labels: str, n: int = 3) -> dict[str, str]:
    return {f"molecule-{i}": array(labels) for i in range(n)}


def test_equal_length_order_discordance_and_repetitive_breakpoint_ambiguity() -> None:
    result = audit(array("ABAC"), reads("ACAB"), M, pairing_status="synthetic")
    assert result.state == "DISCORDANT"
    assert result.label_edit_distance is not None and result.label_edit_distance > 0
    assert result.candidate_event == "equal_length_order_or_label_disagreement"
    deletion = audit(array("ABAB"), reads("ABABAB"), M, pairing_status="synthetic")
    assert deletion.state == "DISCORDANT"
    assert deletion.breakpoint_label_interval is not None
    assert deletion.breakpoint_label_interval[1] - deletion.breakpoint_label_interval[0] > 1


def test_substitution_and_indel_errors_preserve_supported_path() -> None:
    original = array("ABAB")
    substitution = original[:3] + ("A" if original[3] != "A" else "T") + original[4:]
    insertion = original[:13] + "C" + original[13:]
    result = audit(original, {"r1": substitution, "r2": insertion,
                              "r3": original}, M, pairing_status="synthetic")
    assert result.state == "SUPPORTED"
    assert result.assembly.score == 0


def test_close_templates_and_correlated_error_abstain() -> None:
    close = {"A": "ACGTTGCA", "B": "CCGTTGCA"}
    assembly = close["A"] + close["B"] + close["A"]
    all_a = close["A"] * 3
    result = audit(assembly, {f"r{i}": all_a for i in range(3)}, close,
                   pairing_status="synthetic")
    assert result.state == "AMBIGUOUS"
    assert "nonidentifiable" in result.reason


def test_mixed_haplotypes_and_low_support_abstain() -> None:
    mixed = {**reads("ABAB", 2), **{f"alt-{i}": array("ABABAB") for i in range(2)}}
    result = audit(array("ABAB"), mixed, M, pairing_status="synthetic")
    assert result.state == "AMBIGUOUS"
    assert result.reason == "mixed_read_architectures_or_haplotypes"
    low = audit(array("ABAB"), reads("ABABAB", 2), M,
                pairing_status="synthetic")
    assert low.state == "INSUFFICIENT_READ_SUPPORT"


def test_many_unresolved_reads_cannot_be_silently_discarded() -> None:
    molecules = {**reads("ABAB", 3), "bad1": "TTTTTTTT",
                 "bad2": "TTTTTTTT"}
    result = audit(array("ABAB"), molecules, M, pairing_status="synthetic")
    assert result.state == "AMBIGUOUS"
    assert result.reason == "too_many_unresolved_read_paths"
    assert result.ambiguous_count == 2


def test_unverified_pairing_cannot_be_called() -> None:
    result = audit(array("ABAB"), reads("ABABAB"), M,
                   pairing_status="unverified")
    assert result.state == "AMBIGUOUS"
    assert result.reason == "read_pairing_not_verified"


def test_complete_deletion_with_exact_synthetic_flanks() -> None:
    left, right = "TTAAGGCC", "CCTTAAGG"
    row = {"case_id": "opaque_001", "input_status": "ok",
           "array_window_policy": "between_exact_unique_synthetic_flanks",
           "left_flank_sequence": left, "right_flank_sequence": right,
           "assembly_sequence": left + right,
           "raw_read_sequences": {f"r{i}": left + array("ABAB") + right for i in range(3)},
           "candidate_monomers": M, "read_pairing_status": "synthetic_simulated"}
    prediction = predict_case(row)
    assert prediction["status"] == "ok"
    assert prediction["event_score"] == 1.0
    assert prediction["predicted_edited_label_path"] == []
    assert prediction["predicted_edited_interval_bp"] == [len(left), len(left)]


def test_nonunique_flank_abstains_before_sequence_comparison() -> None:
    left = "ACGTTGCA"
    row = {"case_id": "opaque_002", "input_status": "ok",
           "array_window_policy": "between_exact_unique_synthetic_flanks",
           "left_flank_sequence": left, "right_flank_sequence": "GGGGTTTT",
           "assembly_sequence": left + array("ABAB") + "GGGGTTTT",
           "raw_read_sequences": reads("ABAB"), "candidate_monomers": M,
           "read_pairing_status": "synthetic_simulated"}
    prediction = predict_case(row)
    assert prediction["status"] == "abstain"
    assert "flank_absent_or_nonunique" in prediction["reason"]


def test_large_nonhomologous_long_monomers_resolve_with_bounded_dp() -> None:
    a = "ACGTTGCAGGAATCCACATGAGTC" * 7 + "ACG"
    b = "GCTATACCGATGTTACCTGCAAAG" * 7 + "GCT"
    assert len(a) == len(b) == 171
    motifs = {"A": a, "B": b}
    result = decompose(a + b + a + b, motifs)
    assert result.state == "RESOLVED"
    assert tuple(c.label for c in result.copies) == tuple("ABAB")


def test_first_monomer_snp_does_not_invent_new_period_or_structure() -> None:
    original = array("ABABAB")
    changed = "T" + original[1:]
    result = audit(original, {f"r{i}": changed for i in range(3)}, M,
                   pairing_status="synthetic")
    assert result.state == "SUPPORTED"
    assert result.label_edit_distance == 0


def test_homopolymer_indel_can_abstain_without_structural_false_call() -> None:
    motifs = {"A": "ACGTTTTTGCAG", "B": "GGAACCCCTTCA"}
    original = "".join(motifs[x] for x in "ABAB")
    noisy = original.replace("TTTTT", "TTTTTT", 1)
    result = audit(original, {f"r{i}": noisy for i in range(3)}, motifs,
                   pairing_status="synthetic")
    assert result.state in {"SUPPORTED", "AMBIGUOUS", "INSUFFICIENT_READ_SUPPORT"}
    assert result.state != "DISCORDANT"


def test_orientation_inversion_is_separate_from_label_rearrangement() -> None:
    rev_b = M["B"].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    inverted = M["A"] + rev_b + array("AB")
    result = audit(array("ABAB"), {f"r{i}": inverted for i in range(3)}, M,
                   pairing_status="synthetic")
    assert result.state == "DISCORDANT"
    assert result.candidate_event == "orientation_disagreement"


def test_input_indel_and_monomer_variant_do_not_force_exact_bp_estimate() -> None:
    original = array("ABAB")
    changed = original[:12] + "T" + original[13:]
    result = audit(original, {f"r{i}": changed for i in range(3)}, M,
                   pairing_status="synthetic")
    assert result.state in {"SUPPORTED", "AMBIGUOUS", "INSUFFICIENT_READ_SUPPORT"}
    assert result.state != "DISCORDANT"
