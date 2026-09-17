import pytest

from benchmarks.scripts.screen_mo17_native_candidates import shift_identity


def test_shift_identity_detects_exact_period_and_checks_bounds():
    assert shift_identity("ACG" * 6, 3) == 1.0
    assert shift_identity("ACG" * 6, 2) < 0.5
    with pytest.raises(ValueError, match="period out of range"):
        shift_identity("ACG", 3)
