from benchmarks.discovery.saturation import saturation_decision, transition_metrics
from benchmarks.scripts.plot_discovery_saturation import grouped


def test_transition_metrics_count_new_lost_and_jaccard():
    lower = ["ACGT", "AAAAC"]
    higher = ["CGTA", "AAAAC", "CCCCG"]
    observed = transition_metrics(lower, higher, 1, 2)
    assert observed["matched_family_count"] == 2
    assert observed["new_family_count"] == 1
    assert observed["lost_family_count"] == 0
    assert observed["family_jaccard"] == 2 / 3
    assert observed["new_families_per_added_x"] == 1
    assert observed["transition_pass"] is False


def test_saturation_requires_two_transitions_and_every_seed():
    seeds = [1, 2, 3]
    coverages = [1, 2, 5, 10]
    rows = []
    for low, high, passed in [(1, 2, False), (2, 5, True), (5, 10, True)]:
        for seed in seeds:
            rows.append(
                {
                    "seed": seed,
                    "lower_depth": low,
                    "higher_depth": high,
                    "transition_pass": passed,
                }
            )
    decision = saturation_decision(rows, seeds, coverages)
    assert decision["saturation_reached"] is True
    assert decision["saturation_depth"] == 10


def test_one_seed_failure_prevents_saturation():
    seeds = [1, 2, 3]
    coverages = [1, 2, 5]
    rows = []
    for low, high in [(1, 2), (2, 5)]:
        for seed in seeds:
            rows.append(
                {
                    "seed": seed,
                    "lower_depth": low,
                    "higher_depth": high,
                    "transition_pass": not (seed == 3 and high == 5),
                }
            )
    decision = saturation_decision(rows, seeds, coverages)
    assert decision["saturation_reached"] is False
    assert decision["saturation_depth"] is None


def test_curve_summary_uses_median_and_range():
    rows = [
        {"coverage": 1, "value": value} for value in (1, 2, 100)
    ] + [{"coverage": 2, "value": value} for value in (4, 5, 6)]
    depth, center, low, high = grouped(rows, "value")
    assert depth == [1.0, 2.0]
    assert center == [2.0, 5.0]
    assert low == [1.0, 4.0]
    assert high == [100.0, 6.0]
