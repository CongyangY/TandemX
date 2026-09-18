from benchmarks.m2_routes.macadamia_phase_probe import select_phase


def test_phase_probe_finds_rotated_exact_template():
    motif = ("ACGT" * 35) + "AGTC"
    rotated = motif[13:] + motif[:13]
    selected, result = select_phase(rotated * 2, motif)
    assert result["first_144bp_edit_distance"] == 0
    assert selected == rotated
