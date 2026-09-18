from benchmarks.m2_routes.diagnose_macadamia_abstention import array_diagnostic


def test_exact_tiled_input_and_off_phase_diagnostic():
    motif = "ACGT" * 36
    exact = array_diagnostic(motif * 2, motif)
    assert exact["full_tiling_possible"]
    assert exact["first_tile_min_edit"] == 0
    empty = array_diagnostic("", motif)
    assert empty["full_tiling_possible"]
    assert empty["no_full_tiling_interpretation"] == "empty_interval"
