"""CIGAR boundary projection protects source-guided flank spans."""

import pytest

from benchmarks.scripts.score_macadamia_native_span import projected_boundaries


def test_primary_context_cigar_projects_both_orientations():
    base = ["read", "1200", "100", "1100", "+", "context", "1000", "0", "1000",
            "1000", "1000", "60", "tp:A:P", "cg:Z:1000M"]
    assert projected_boundaries("\t".join(base), 100, 200) == [200, 300]
    base[4] = "-"
    assert projected_boundaries("\t".join(base), 100, 200) == [200, 300]
    base[12] = "tp:A:S"
    with pytest.raises(ValueError, match="primary"):
        projected_boundaries("\t".join(base), 100, 200)
