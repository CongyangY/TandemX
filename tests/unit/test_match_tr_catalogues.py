from benchmarks.scripts.match_tr_catalogues import classify_match


def test_classify_direct_recurrence_requires_identity_and_similar_length():
    state, multiple, error = classify_match(0.97, 1.04, 0.95, 0.90, 1.10, 0.90, 0.05)
    assert state == "direct_catalogue_recurrence"
    assert multiple == 1
    assert abs(error - 0.04) < 1e-12


def test_classify_related_integer_multiple_separately():
    state, multiple, error = classify_match(0.93, 2.02, 0.95, 0.90, 1.10, 0.90, 0.05)
    assert state == "related_period_multiple"
    assert multiple == 2
    assert error < 0.011


def test_classify_nonmatching_sequence():
    state, _, _ = classify_match(0.81, 1.0, 0.95, 0.90, 1.10, 0.90, 0.05)
    assert state == "no_qualifying_match"
