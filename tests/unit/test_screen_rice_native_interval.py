from benchmarks.scripts.screen_rice_native_interval import (
    candidate_scores,
    entropy,
    flank_unique_fraction,
    select,
)


def test_candidate_scores_identifies_exact_periodic_window() -> None:
    sequence = b"ACGTACGT" + b"AACCGGTT" * 4 + b"TGCATGCA"
    scores = list(candidate_scores(sequence, period=8, length=32, stride=4, flank=8))
    assert scores == [(8, 24, 1.0)]


def test_entropy_and_local_flank_uniqueness() -> None:
    assert entropy(b"ACGT" * 10) == 2.0
    context = b"ACGTCCGATTCGATGC"
    assert flank_unique_fraction(context, 4, 12, k=2) > 0.0


def test_select_keeps_failed_top_candidate_and_uses_fixed_rules() -> None:
    sequence = b"ACGT" * 40
    config = {"interval_selection_before_read_mapping": {
        "flank_bp_each_side": 40,
        "array_length_bp": 40,
        "minimum_period_shift_identity": 0.8,
        "minimum_base_shannon_entropy_bits_per_base": 1.5,
        "minimum_flank_31mer_uniqueness_fraction_within_context": 0.0,
        "max_ranked_candidates_to_screen": 100,
    }}
    checked, context = select([(0.9, "chr1", 40, 30)], {"chr1": sequence}, config)
    assert checked[0][7]
    assert context == sequence[0:120]
