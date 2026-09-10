"""Development-only A3 local-phase semantics; no held-out performance claim."""

import random

import pytest

from benchmarks.challenge.local_phase import LocalPhaseConfig, LocalPhasePrototype
from tandemx.utils.kmers import reverse_complement


def _unit(seed: int = 9) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(120))


def test_one_base_insertion_relocks_local_phase() -> None:
    unit = _unit()
    read = unit * 2 + unit[:60] + "A" + unit[60:] + unit
    observed = LocalPhasePrototype({"f": unit}).intervals(read)
    assert observed == {"f": [(0, len(read))]}


def test_twenty_base_insertion_is_split_not_whole_read_drop() -> None:
    unit = _unit()
    read = unit * 2 + unit[:60] + "A" * 20 + unit[60:] + unit
    observed = LocalPhasePrototype({"f": unit}).intervals(read)
    assert observed["f"]
    assert len(observed["f"]) > 1
    assert observed["f"][0][0] == 0


def test_strict_zero_drift_is_a_conservative_control() -> None:
    unit = _unit()
    read = unit * 2 + unit[:60] + "A" + unit[60:] + unit
    relaxed = LocalPhasePrototype({"f": unit}).intervals(read)
    strict = LocalPhasePrototype({"f": unit}, LocalPhaseConfig(maximum_drift=0)).intervals(read)
    assert relaxed["f"] == [(0, len(read))]
    assert strict["f"] != relaxed["f"]


def test_shared_80bp_background_is_not_a_complete_family() -> None:
    unit = _unit()
    read = unit[:80] * 30
    assert LocalPhasePrototype({"f": unit}).intervals(read) == {}


def test_reverse_complement_is_strand_invariant() -> None:
    unit = _unit()
    detector = LocalPhasePrototype({"f": unit})
    assert detector.intervals(reverse_complement(unit * 4)) == {"f": [(0, 480)]}


def test_close_families_keep_ambiguous_candidate_intervals() -> None:
    unit = _unit()
    near = list(unit)
    for pos in range(10):
        near[pos] = "A" if near[pos] != "A" else "C"
    near = "".join(near)
    observed = LocalPhasePrototype({"a": unit, "b": near}).intervals(unit * 3)
    assert set(observed) == {"a", "b"}
    assert observed["a"] == [(0, 360)]
    assert observed["b"]


@pytest.mark.parametrize(
    "config",
    [
        LocalPhaseConfig(maximum_seed_gap=0),
        LocalPhaseConfig(maximum_drift=-1),
        LocalPhaseConfig(predecessor_cap=0),
        LocalPhaseConfig(reset_penalty=-1),
        LocalPhaseConfig(phase_fraction=1.1),
        LocalPhaseConfig(minimum_units=-1),
    ],
)
def test_invalid_config_fails_fast(config: LocalPhaseConfig) -> None:
    with pytest.raises(ValueError):
        LocalPhasePrototype({"f": _unit()}, config)

