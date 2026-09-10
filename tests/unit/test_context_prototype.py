"""Unit checks for the development-only context adversarial fixture.

These are fixture semantics and regression checks, not a holdout benchmark.
"""

import random

import pytest

from benchmarks.challenge.context_prototype import (
    ContextConfig,
    PeriodicContextPrototype,
    merge_intervals,
    partition_evidence,
)
from tandemx.utils.kmers import reverse_complement


def _monomer(seed: int = 17, length: int = 120) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def test_empty_and_invalid_catalogues_are_rejected() -> None:
    with pytest.raises(ValueError, match="Catalogue is empty"):
        PeriodicContextPrototype({})
    with pytest.raises(ValueError, match="Catalogue requires"):
        PeriodicContextPrototype({"bad": "ACGTN"})
    with pytest.raises(ValueError, match="Catalogue requires"):
        PeriodicContextPrototype({"short": "ACGT"}, ContextConfig(k=11))


def test_complete_six_copy_array_and_reverse_complement_are_detected() -> None:
    unit = _monomer()
    detector = PeriodicContextPrototype({"m120": unit})
    assert detector.intervals(unit * 6) == {"m120": [(0, 720)]}
    assert detector.intervals(reverse_complement(unit * 6)) == {
        "m120": [(0, 720)]
    }


def test_short_local_tandem_decoys_are_not_complete_120bp_family() -> None:
    unit = _monomer()
    detector = PeriodicContextPrototype({"m120": unit})
    # Both are genuine local tandem sequences, but neither spans enough phases
    # of the complete 120-bp catalogue unit to be an m120 call.
    assert detector.intervals(unit[:24] * 10) == {}
    assert detector.intervals(unit[:80] * 4) == {}


def test_nonoverlapping_families_in_one_read_are_counted_separately() -> None:
    a = _monomer(17)
    b = _monomer(31)
    detector = PeriodicContextPrototype({"a": a, "b": b})
    read = a * 3 + "N" * 20 + b * 3
    intervals = detector.intervals(read)
    assert intervals["a"] == [(0, 360)]
    assert intervals["b"] == [(380, 740)]
    counts = partition_evidence(intervals, len(read))
    assert counts["unique_bp"] == {"a": 360, "b": 360}
    assert counts["ambiguous_bp"] == 0


def test_g_flank_retains_known_one_base_background_boundary() -> None:
    a = _monomer(17)
    b = _monomer(31)
    detector = PeriodicContextPrototype({"a": a, "b": b})
    intervals = detector.intervals(a * 3 + "G" * 20 + b * 3)
    # Development observation: one G-flank base is admitted to b's span.
    # Preserve this boundary evidence rather than hiding it with a separator.
    assert intervals["b"] == [(379, 740)]


def test_near_identical_families_retain_overlap_as_ambiguous() -> None:
    a = _monomer()
    b = list(a)
    b[20] = "A" if b[20] != "A" else "C"
    b = "".join(b)
    detector = PeriodicContextPrototype({"a": a, "b": b})
    intervals = detector.intervals(a * 3)
    assert set(intervals) == {"a", "b"}
    counts = partition_evidence(intervals, 360)
    assert counts["ambiguous_bp"] > 0
    assert counts["eligible_bp"]["a"] == 360
    assert counts["eligible_bp"]["b"] > 0


def test_partition_mass_is_conserved_without_double_counting_overlap() -> None:
    intervals = {"a": [(0, 200)], "b": [(100, 300)]}
    counts = partition_evidence(intervals, 400)
    assert counts["unique_bp"] == {"a": 100, "b": 100}
    assert counts["ambiguous_bp"] == 100
    assert counts["unassigned_bp"] == 100
    assert (
        sum(counts["unique_bp"].values())
        + counts["ambiguous_bp"]
        + counts["unassigned_bp"]
        == 400
    )


def test_interval_bounds_are_checked() -> None:
    with pytest.raises(ValueError, match="Interval exceeds read length"):
        partition_evidence({"m120": [(0, 121)]}, 120)
    with pytest.raises(ValueError, match="Invalid half-open interval"):
        merge_intervals([(3, 3)])
    with pytest.raises(ValueError, match="Invalid half-open interval"):
        merge_intervals([(-1, 2)])
