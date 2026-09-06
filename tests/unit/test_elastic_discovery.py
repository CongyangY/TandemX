"""Alignment and consensus checks independent of the benchmark simulator."""
import random

import pytest

from tandemx.discover.alignment import banded_self_align, banded_self_align_many, global_align_ops
from tandemx.discover.consensus import aligned_unit_consensus
from tandemx.discover.elastic import discover_elastic_arrays
from tandemx.discover.rust_backend import rust_backend_available
from tandemx.simulate.toy import reverse_complement


def sequence(seed, length):
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


@pytest.mark.parametrize("backend", ["python", "rust"])
def test_insertions_deletions_and_two_arrays(backend):
    if backend == "rust" and not rust_backend_available():
        pytest.skip("compiled extension unavailable")
    unit = sequence(7, 71)
    units = [unit, unit[:23] + "G" + unit[23:], unit, unit[:43] + unit[44:], unit] * 2
    first = "".join(units)
    second = sequence(9, 53) * 8
    read = sequence(3, 200) + first + "N" * 150 + second + sequence(5, 200)
    calls, _ = discover_elastic_arrays(read, min_period=30, max_period=200, min_span=100, backend=backend)
    assert len(calls) == 2
    hit, consensus, _ = calls[0]
    assert abs(hit.start - 200) <= 5
    assert abs(hit.end - 200 - len(first)) <= 5
    assert hit.gaps > 0
    assert consensus in unit * 2 or unit in consensus * 2
    assert len(consensus) == len(unit)
    assert calls[1][0].start >= 200 + len(first) + 140


@pytest.mark.parametrize("backend", ["python", "rust"])
def test_negative_and_short_inputs(backend):
    if backend == "rust" and not rust_backend_available():
        pytest.skip("compiled extension unavailable")
    for read in ["", "N" * 800, sequence(18, 900), "AC" * 400]:
        assert discover_elastic_arrays(read, min_period=30, max_period=200, min_span=100, backend=backend)[0] == []


def test_parameter_validation():
    for kwargs in [dict(period=0), dict(period=50, band=50), dict(period=50, band=-1), dict(period=50, x_drop=0)]:
        with pytest.raises(ValueError):
            banded_self_align("ACGT" * 100, min_span=100, **kwargs)
    with pytest.raises(ValueError, match="ASCII"):
        banded_self_align("非DNA", 10, 20)


def test_composition_filter_keeps_true_at_rich_array():
    rng = random.Random(76)
    unit = "".join(rng.choices("AAAAAAAATTTTTTTTCG", k=87))
    read = sequence(11, 200) + unit * 8 + sequence(13, 150)
    hits, _ = discover_elastic_arrays(read, min_period=30, max_period=200, min_span=100)
    assert len(hits) == 1
    assert len(hits[0][1]) == 87
    # A few matching flank bases have no uniquely identifiable biological edge.
    hit = hits[0][0]
    overlap = min(hit.end, 200 + len(unit) * 8) - max(hit.start, 200)
    union = max(hit.end, 200 + len(unit) * 8) - min(hit.start, 200)
    assert overlap / union >= 0.98


def test_short_toy_repeat_is_supported_with_explicit_span():
    calls, _ = discover_elastic_arrays("ACGT" * 4, min_period=4, max_period=8, min_span=8)
    assert len(calls) == 1
    assert calls[0][0].start == 0
    assert calls[0][0].end == 16
    assert calls[0][1] == "ACGT"


@pytest.mark.skipif(not rust_backend_available(), reason="compiled extension unavailable")
def test_native_reference_parity_and_reverse_complement():
    unit = sequence(14, 67)
    read = sequence(4, 80) + (unit + unit[:21] + "T" + unit[21:] + unit) * 3 + sequence(3, 70)
    python = banded_self_align(read, 65, 100)
    assert python == banded_self_align(read, 65, 100, backend="rust")
    assert python
    reverse = banded_self_align(reverse_complement(read), 65, 100, backend="rust")
    assert len(reverse) == 1
    assert abs(reverse[0].start - (len(read) - python[0].end)) <= 3
    a = aligned_unit_consensus(read, python[0])[0]
    b = aligned_unit_consensus(read, python[0], backend="rust")[0]
    assert a == b


@pytest.mark.skipif(not rust_backend_available(), reason="compiled extension unavailable")
def test_native_batched_alignment_matches_individual_period_calls():
    first = sequence(70, 67)
    second = sequence(71, 101)
    read = sequence(72, 90) + first * 6 + "N" * 50 + second * 5 + sequence(73, 80)
    periods = [65, 67, 99, 101]
    expected = [
        hit
        for period in periods
        for hit in banded_self_align(read, period, 100, backend="rust")
    ]
    assert banded_self_align_many(read, periods, 100, backend="rust") == expected


@pytest.mark.skipif(not rust_backend_available(), reason="compiled extension unavailable")
@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("period", [1, 7, 31, 67])
def test_native_reused_rows_match_reference_across_path_breaks(seed, period):
    """Check entire tracebacks, including rejected paths and stale row hazards."""
    rng = random.Random(9000 + seed + period)
    unit = "".join(rng.choices("ACGT", k=period))
    bases = list(unit * 7)
    for index in range(seed, len(bases), 13):
        bases[index] = rng.choice("ACGTNRY")
    bases[len(bases) // 2:len(bases) // 2] = list("NryN" if seed % 2 else "ACG")
    read = sequence(seed + 100, 23) + "".join(bases) + sequence(seed + 200, 19)
    if seed % 2:
        read = read.lower()
    for band in sorted({0, min(4, period - 1), period - 1}):
        for x_drop in (1, 40):
            expected = banded_self_align(read, period, max(2, period * 2),
                                         band=band, x_drop=x_drop)
            assert banded_self_align(read, period, max(2, period * 2), band=band,
                                     x_drop=x_drop, backend="rust") == expected


@pytest.mark.parametrize("reference,query", [("", ""), ("A", ""), ("", "AC"),
    ("ACGTACGT", "ACGTTACGT"), ("ACGTACGT", "ACTACGT"), ("AGCGATTCGAT", "AGGATTCAAT")])
def test_global_alignment_consumes_both_sequences(reference, query):
    ops = global_align_ops(reference, query)
    assert sum(op in "MD" for op in ops) == len(reference)
    assert sum(op in "MI" for op in ops) == len(query)
    if rust_backend_available():
        assert ops == global_align_ops(reference, query, backend="rust")
