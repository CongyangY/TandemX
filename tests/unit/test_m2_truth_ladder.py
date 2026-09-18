from benchmarks.m2_routes.truth_ladder import assemble, mutate_copy


def test_fixed_noise_preserves_length_and_changes_expected_bases():
    original = "A" * 96
    changed = mutate_copy(original)
    assert len(changed) == 96
    assert changed != original
    assert changed[15] == "C"
    assert changed[70] == "C"


def test_assemble_reverses_only_declared_copy():
    sequence, labels = assemble("A B- A", {"A": "ACGT", "B": "AAAA"}, False)
    assert sequence == "ACGTTTTTACGT"
    assert labels == ["A+", "B-", "A+"]
