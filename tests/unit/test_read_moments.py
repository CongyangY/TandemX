from collections import Counter
import json
import math
import random

import pytest

from tandemx.quantify.read_moments import FamilyMoments, ReadMoments, collect_read_moments
from tandemx.utils.kmers import reverse_complement


def test_cluster_standard_error_matches_independent_residual_calculation():
    x, y = [10, 20, 40, 30], [1, 2, 3, 0]
    state = ReadMoments(families={'f': FamilyMoments()})
    for exposure, value in zip(x, y):
        state.add(exposure+2, 3, {'f': value}, exposure)
    estimate = state.estimates(1000, {'f': 1}, minimum_effective_reads=2)[0]
    ratio = sum(y)/sum(x)
    residuals = [yi-ratio*xi for xi, yi in zip(x, y)]
    independent_se = 1000*math.sqrt(sum(r*r for r in residuals)/(len(x)-1)*len(x))/sum(x)
    assert estimate.estimated_copy_number == 60
    assert estimate.sampling_standard_error == pytest.approx(independent_se)
    assert estimate.effective_positive_reads == pytest.approx(36/14)
    assert estimate.sampling_interval_low == pytest.approx(max(0, 60-1.959963984540054*independent_se))
    assert estimate.sampling_interval_high == pytest.approx(60+1.959963984540054*independent_se)


def test_counts_match_independent_naive_scan_with_multiplicity_and_reverse_strand():
    sequences = ['ACGACGCTACTANACG', 'TAGCGTCGT', 'AC', 'NNNN', 'AAAAAAA']
    diagnostic = {'f': {'ACG': 2, 'CTA': 1}, 'absent': {'AGC': 1}, 'undetermined': {}}
    state = collect_read_moments(iter(sequences), diagnostic, 3)
    values = []
    for sequence in sequences:
        counts = Counter()
        for i in range(len(sequence)-2):
            word = sequence[i:i+3]
            if set(word) <= set('ACGT'):
                rc = word.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                counts[min(word, rc)] += 1
        values.append((counts['ACG']/2+counts['CTA'])/2)
    assert state.families['f'].sum_y == sum(values)
    assert state.families['f'].sum_y2 == sum(v*v for v in values)
    assert state.sum_x == sum(max(0, len(s)-2) for s in sequences)
    assert state.valid_windows < state.sum_x and state.short_reads == 1
    other = collect_read_moments((reverse_complement(s) for s in sequences), diagnostic, 3)
    assert state == other
    estimates = {e.family_id: e for e in state.estimates(1000, {f: len(d) for f, d in diagnostic.items()})}
    assert estimates['undetermined'].estimated_copy_number is None
    assert estimates['undetermined'].sampling_interval_high is None
    assert 'ambiguous_base_windows' in estimates['f'].warning


def test_sparse_zero_and_no_exposure_do_not_get_false_precise_intervals():
    state = collect_read_moments(['AAAA']*30, {'f': {'ACG': 1}}, 3)
    absent = state.estimates(1000, {'f': 1})[0]
    assert absent.estimated_copy_number == 0
    assert absent.interval_status == 'no_observed_support_not_absence'
    assert absent.sampling_interval_high is None
    uniform = collect_read_moments(['ACG']*30, {'f': {'ACG': 1}}, 3).estimates(1000, {'f': 1})[0]
    assert uniform.interval_status == 'zero_observed_variance_uncalibrated'
    assert uniform.sampling_interval_high is None
    sparse = collect_read_moments(['AAAA']*99+['ACG'], {'f': {'ACG': 1}}, 3).estimates(1000, {'f': 1})[0]
    assert sparse.estimated_copy_number > 0 and sparse.interval_status == 'insufficient_effective_reads'
    assert sparse.sampling_standard_error is None
    for sequences in ([], ['AC']*10):
        with pytest.raises(ValueError, match='usable'):
            collect_read_moments(sequences, {'f': {'ACG': 1}}, 3).estimates(1000, {'f': 1})


def test_invalid_targets_and_input():
    for diagnostic in ({'a': {'ACG': 0}}, {'a': {'ACG': 1.2}}, {'a': {'CGT': 1}},
                       {'a': {'ACN': 1}}, {'a': {'ACG': 1}, 'b': {'ACG': 1}}):
        with pytest.raises(ValueError):
            collect_read_moments(['ACG'], diagnostic, 3)
    for sequences in ([''], ['ACZ'], ['ACG', 'Aé']):
        with pytest.raises(ValueError):
            collect_read_moments(sequences, {'a': {'ACG': 1}}, 3)


def test_read_cluster_interval_calibration_under_independent_bernoulli_model():
    # This checks the statistical formula in its own assumptions, not biological calibration.
    generator = random.Random(341)
    intervals = []
    for _ in range(500):
        state = ReadMoments(families={'f': FamilyMoments()})
        for _ in range(400):
            state.add(1002, 3, {'f': 10 if generator.random() < .25 else 0}, 1000)
        estimate = state.estimates(10_000, {'f': 1})[0]
        assert estimate.sampling_interval_high is not None
        intervals.append(estimate.sampling_interval_low <= 25 <= estimate.sampling_interval_high)
    assert .90 <= sum(intervals)/len(intervals) <= .99


def test_replay_independent_generated_data_and_reject_changed_reads(tmp_path):
    from benchmarks.abundance.simulate import GenomeSpec, build_genome, write_genome, sample_reads
    from benchmarks.challenge.schema import read_table
    from benchmarks.scripts.evaluate_read_moments import worker

    baseline = tmp_path/'baseline'
    genome_dir = baseline/'genomes/s4101'
    spec = GenomeSpec(seed=4101, periods=(31,), copies=(20,), flank_bp=1000)
    write_genome(spec, genome_dir, (1,))
    genome, _, truth = build_genome(spec)
    for error in (0, .01):
        sample_reads(genome, truth, baseline/'reads/s4101'/f'c10_e{error}', seed=1004104,
                     coverage=10, read_length=260, substitution_rate=error)
    (baseline/'environment.json').write_text('{"split":"development"}')
    (baseline/'run_config.json').write_text(json.dumps(dict(
        k=13, seeds={'development': [4101]}, coverages=[10], substitution_rates=[0, .01])))
    output = tmp_path/'output'
    output.mkdir()
    worker(baseline, output)
    rows = read_table(output/'metrics.tsv')
    assert len(rows) == 2 and {r['truth_copies'] for r in rows} == {'20'}
    summary = read_table(output/'summary.tsv')
    assert len(summary) == 2
    assert all(int(r['intervals_available'])+int(r['intervals_unavailable']) == 1 for r in summary)
    source = baseline/'reads/s4101/c10_e0/reads.fa'
    source.write_text(source.read_text()+'\n')
    with pytest.raises(ValueError, match='Read hash'):
        worker(baseline, output)
    (baseline/'environment.json').write_text('{"split":"heldout"}')
    with pytest.raises(ValueError, match='held-out'):
        worker(baseline, output)
