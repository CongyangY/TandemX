import math

from benchmarks.scripts.independent_verify_tidecluster_factorial import (
    base_union_metrics,
    canonical_monomer,
    maximum_matching,
    reaches_cyclic_threshold,
    same_value,
)


def test_independent_interval_union_and_matching() -> None:
    predicted = [
        {"chrom": "chr1", "start": 10, "end": 30},
        {"chrom": "chr1", "start": 20, "end": 40},
    ]
    truth = [{"chrom": "chr1", "start": 15, "end": 35}]
    recall, precision = base_union_metrics(predicted, truth)
    assert recall == 1.0
    assert precision == 20 / 30
    assert maximum_matching([[0], [0]]) == {0: 0}


def test_independent_cyclic_threshold_handles_rotation_and_reverse_complement() -> None:
    truth = "AACCGT"
    assert canonical_monomer(truth) == canonical_monomer("CGTAAC")
    assert reaches_cyclic_threshold(truth, "CGTAAC")
    assert reaches_cyclic_threshold(truth, "ACGGTT")
    assert not reaches_cyclic_threshold(truth, "TTTTTT")


def test_same_value_handles_json_null_as_nan() -> None:
    assert same_value(None, math.nan)
    assert same_value("0.5", 0.5)
    assert not same_value("0.6", 0.5)
