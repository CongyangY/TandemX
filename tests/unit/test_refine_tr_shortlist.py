from benchmarks.scripts.refine_tr_shortlist import classify_additional_exclusion


def test_additional_exclusion_keeps_only_explicit_no_match():
    assert classify_additional_exclusion("no_match_in_limited_library") == "retained_preliminary_candidate"
    assert classify_additional_exclusion("strong_known_clone_match") == "excluded_known_or_related_sequence"
    assert (
        classify_additional_exclusion("insufficient_length_for_exclusion")
        == "excluded_unresolved_additional_library"
    )
