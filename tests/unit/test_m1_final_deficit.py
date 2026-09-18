"""Checks for the bounded final M1 scoring boundary."""
import pytest

from benchmarks.m1_final_deficit.run import read_derived_estimates, summarize
from tandemx.io.sequences import SequenceRecord


def test_exact_flank_span_requires_three_distinct_molecules() -> None:
    unit = "ACGTTGCA" * 10
    left, right = "T" * 32, "C" * 32
    genome = "A" * 64 + left + unit * 2 + right + "G" * 64
    start = 64 + 32
    truth = [dict(family_id="f1", start=start, end=start + len(unit) * 2)]
    records = [SequenceRecord(id=f"r{i}", sequence=genome) for i in range(3)]
    estimate, support = read_derived_estimates(records, {"f1": unit}, 3.0, genome, truth)
    assert estimate["robust_read_span"]["f1"] == 160
    assert support["f1"]["distinct_read_count"] == 3
    estimate, _ = read_derived_estimates(records[:2], {"f1": unit}, 2.0, genome, truth)
    assert estimate["robust_read_span"]["f1"] is None
    with pytest.raises(ValueError, match="Duplicate read identifier"):
        read_derived_estimates(records[:1] * 3, {"f1": unit}, 3.0, genome, truth)


def test_fixed_read_edit_cannot_establish_calibration() -> None:
    source, estimate = 200, 120
    rows = []
    for fraction in (1.0, .5, 0.0):
        assembly = int(source*fraction)
        injected = source-assembly
        prediction = estimate-assembly
        rows.append(dict(method="production_quantify", panel_id="s1_v00_c30",
                         family="f1", assembly_fraction=fraction,
                         source_bp=source, read_estimate_bp=estimate,
                         injected_missing_bp=injected, predicted_signed_deficit_bp=prediction,
                         signed_error_bp=prediction-injected,
                         absolute_error_bp=abs(prediction-injected),
                         relative_error_to_source=abs(prediction-injected)/source))
    summary = summarize(rows, {})
    assert summary["methods"]["production_quantify"]["a_gate_pass"] is False
    assert summary["b_gate_assessable_from_fixed_read_edits"] is False
    assert len({row["signed_error_bp"] for row in rows}) == 1
