import pytest

from benchmarks.challenge.unified_correspondence import match_native_catalogue


def test_rotation_and_reverse_complement_are_one_unique_match():
    truth = {"F1": "AACCGT"}
    native = {"native": "ACGGTT"}  # rotation of the reverse complement.

    result = match_native_catalogue(native, truth)

    assert result["native"].status == "unique"
    assert result["native"].matches == ("F1",)
    assert result["native"].identities == (("F1", 1.0),)


@pytest.mark.parametrize("copies", [2, 3])
def test_pure_integer_multiple_hor_matches_truth_unit(copies):
    truth = {"F1": "AACCGGTT"}

    result = match_native_catalogue({"native": truth["F1"] * copies}, truth)

    assert result["native"].status == "unique"
    assert result["native"].matches == ("F1",)


def test_small_indel_can_meet_declared_global_identity_threshold():
    truth = {"F1": "AACCGGTTAACCGGTTAACC"}
    native = {"native": "AACCGGTTAACGGTTAACC"}

    result = match_native_catalogue(native, truth, threshold=0.9)

    assert result["native"].status == "unique"
    assert result["native"].identities[0][1] == pytest.approx(19 / 20)


def test_unrelated_and_detectable_mixed_hor_are_not_forced_to_a_family():
    truth = {"F1": "AAAA", "F2": "CCCC"}
    result = match_native_catalogue({"unrelated": "ATAT", "mixed": "AAAACCCC"}, truth)

    assert result["unrelated"].status == "unmatched"
    assert result["mixed"].status == "unmatched_or_composite"
    assert result["mixed"].matches == ()


def test_multiple_truth_matches_are_ambiguous_without_best_score_selection():
    truth = {"F1": "AACCGGTTAACC", "F2": "AACCGGTTAACG"}

    result = match_native_catalogue({"native": "AACCGGTTAACC"}, truth, threshold=0.9)

    assert result["native"].status == "ambiguous"
    assert result["native"].matches == ("F1", "F2")


def test_65x_repeat_is_out_of_scope_at_64x_limit():
    truth = {"F1": "ACGT"}

    result = match_native_catalogue({"native": truth["F1"] * 65}, truth, max_multiple=64)

    assert result["native"].status == "out_of_scope"
    assert result["native"].matches == ()
