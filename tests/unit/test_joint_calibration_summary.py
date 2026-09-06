from benchmarks.scripts.evaluate_joint_multik import summarize


def test_missing_intervals_are_distinct_from_failed_coverage_and_seeds_remain_separate():
    common = dict(coverage='20', error_model='error_free', unit_substitution_rate='0', array_scope='factorial')
    rows = [dict(common, seed=6301, sampling_interval_low=1, truth_covered=True, relative_interval_width=.2),
            dict(common, seed=6301, sampling_interval_low=1, truth_covered=False, relative_interval_width=.4),
            dict(common, seed=6301, sampling_interval_low=None, truth_covered=None, relative_interval_width=None),
            dict(common, seed=6302, sampling_interval_low=None, truth_covered=None, relative_interval_width=None)]
    a, b = summarize(rows)
    assert a['family_conditions'] == 3 and a['intervals_available'] == 2
    assert a['conditional_coverage_fraction'] == .5
    assert a['available_and_covering_fraction'] == 1/3
    assert b['conditional_coverage_fraction'] is None
    assert b['intervals_available'] == 0
    assert 'dependent_families' in a['warning']
