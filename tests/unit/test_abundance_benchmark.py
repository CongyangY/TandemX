"""Check truth generation and score semantics without invoking detector internals."""
import json
from pathlib import Path

import pytest

from benchmarks.abundance.evaluate import finite_nonnegative, score_copy_number, score_localization, score_comparison
from benchmarks.abundance.simulate import GenomeSpec, build_genome, sample_reads
from benchmarks.abundance.run import (
    aggregate,
    run,
    validate_frozen_classifier,
    validate_frozen_localizer,
)
from benchmarks.abundance.evaluate_multik_collapse import run as run_multik_collapse
from benchmarks.challenge.schema import digest_file, read_table


def test_uniform_sampler_reproducible_and_repeat_occupancy_matches_base_oracle(tmp_path: Path):
    spec = GenomeSpec(17, (11, 17), (4, 6), 31)
    genome, monomers, truth = build_genome(spec)
    assert len(genome) == 3*31 + 11*4 + 17*6
    for row in truth:
        assert genome[row['start']:row['end']] == monomers[row['family_id']] * row['copies']
    manifests = []
    for name, error in [('a', 0), ('b', 0), ('c', .1)]:
        manifests.append(sample_reads(genome, truth, tmp_path/name, seed=8, coverage=3, read_length=50, substitution_rate=error))
    assert (tmp_path/'a'/'reads.fa').read_bytes() == (tmp_path/'b'/'reads.fa').read_bytes()
    coords = [read_table(tmp_path/name/'sampling.tsv') for name in ('a','c')]
    assert [(r['genome_start'],r['strand']) for r in coords[0]] == [(r['genome_start'],r['strand']) for r in coords[1]]
    for row in truth:
        positions = set(range(row['start'], row['end']))
        expected = sum((int(r['genome_start'])+j) % len(genome) in positions for r in coords[0] for j in range(50))
        assert manifests[0]['sampled_repeat_bp'][row['family_id']] == expected
    assert manifests[0]['observed_substitutions'] == 0 < manifests[2]['observed_substitutions']
    assert {r['strand'] for r in coords[0]} == {'+','-'}
    assert any(int(r['genome_start'])+50 > len(genome) for r in coords[0])


def test_assemblies_preserve_flanks_and_mononers(tmp_path: Path):
    spec = GenomeSpec(3, (11,17), (4,8), 30)
    full, motifs, truth = build_genome(spec)
    collapsed, same, retained = build_genome(spec, .5)
    assert motifs == same
    assert [r['copies'] for r in retained] == [2,4]
    for first, second, seq in ((truth, retained, collapsed),):
        assert full[:first[0]['start']] == seq[:second[0]['start']]
        assert full[first[0]['end']:first[1]['start']] == seq[second[0]['end']:second[1]['start']]
        assert full[first[1]['end']:] == seq[second[1]['end']:]
    with pytest.raises(ValueError, match='limit'):
        GenomeSpec(1, max_genome_bp=100).validate()
    with pytest.raises(ValueError):
        sample_reads(full, truth, tmp_path/'bad', seed=1, coverage=float('nan'), read_length=10, substitution_rate=0)


def test_divergent_fragmented_arrays_are_deterministic_and_keep_copy_truth():
    spec = GenomeSpec(
        23, (11, 17), (9, 12), 30,
        unit_substitution_rate=.2, array_fragments=3, fragment_gap_bp=7,
    )
    genome, founders, truth = build_genome(spec)
    repeated, same_founders, same_truth = build_genome(spec)
    assert (genome, founders, truth) == (repeated, same_founders, same_truth)
    assert len(truth) == 6
    for family_id, copies in [('f1', 9), ('f2', 12)]:
        rows = [row for row in truth if row['family_id'] == family_id]
        assert len(rows) == 3
        assert sum(row['copies'] for row in rows) == copies
        assert sum(row['repeat_bp'] for row in rows) == copies * len(founders[family_id])
        assert all(rows[index+1]['start'] - rows[index]['end'] == 7 for index in range(2))
        observed = ''.join(genome[row['start']:row['end']] for row in rows)
        assert observed != founders[family_id] * copies

    collapsed, collapsed_founders, retained = build_genome(spec, .5)
    assert collapsed_founders == founders
    assert len(collapsed) < len(genome)
    assert {family: sum(row['copies'] for row in retained if row['family_id'] == family)
            for family in founders} == {'f1': 5, 'f2': 6}

    absent, _, missing = build_genome(spec, 0)
    assert len(missing) == 2
    assert all(row['start'] == row['end'] and row['copies'] == 0 for row in missing)
    assert len(absent) == 3 * spec.flank_bp

    with pytest.raises(ValueError, match='fragmentation'):
        GenomeSpec(1, unit_substitution_rate=1).validate()


def test_scoring_does_not_confuse_sampling_oracle_or_zero_denominators(tmp_path: Path):
    truth = [dict(chrom='chr_sim', family_id='f1', start=10, end=30, copies=2, repeat_bp=20, period=10)]
    path = tmp_path/'cn.tsv'
    path.write_text('family_id\testimated_copy_number\tcopy_number_interval_low\tcopy_number_interval_high\n'
                    'f1\t1\t0.9\t1.1\n')
    scored = score_copy_number(path, truth, dict(sampled_repeat_bp={'f1':10},actual_base_coverage=1))
    assert scored[0]['signed_relative_error'] == -.5
    assert scored[0]['estimator_minus_sampling_oracle'] == 0
    assert not scored[0]['interval_contains_truth']
    assert aggregate([dict(coverage=1, **scored[0])], ['coverage'], 'copy')[0]['interval_empirical_coverage'] == 0
    bed = tmp_path/'arrays.bed'
    bed.write_text('chr_sim\t10\t25\tf1\nchr_sim\t20\t30\tf1\n')
    assert score_localization(bed, truth, 50)[0]['base_precision'] == 1
    assert score_localization(bed, truth, 50)[0]['base_recall'] == 1
    assert score_localization(bed, truth, 50)[0]['fragments'] == 1
    comp = tmp_path/'comp.tsv'
    comp.write_text('family_id\tassembly_read_ratio\tstatus\nf1\t0\treads_only\n')
    missing = [dict(truth[0], repeat_bp=0, copies=0, end=10)]
    row = score_comparison(comp, truth, missing)[0]
    assert row['outcome'] == 'TP' and row['native_status'] == 'reads_only'
    summary = aggregate([row], [], 'compare')[0]
    assert summary['false_positive_rate'] is None
    for value in ('nan','inf','-1'):
        with pytest.raises(ValueError):
            finite_nonnegative(value, 'test')
    path.write_text(path.read_text() + 'f1\t1\t0.9\t1.1\n')
    with pytest.raises(ValueError, match='families'):
        score_copy_number(path, truth, {})


def test_scoring_aggregates_multiple_truth_intervals_per_family(tmp_path: Path):
    truth = [
        dict(chrom='chr_sim', family_id='f1', start=10, end=30, copies=2, repeat_bp=20, period=10),
        dict(chrom='chr_sim', family_id='f1', start=40, end=70, copies=3, repeat_bp=30, period=10),
    ]
    bed = tmp_path/'arrays.bed'
    bed.write_text('chr_sim\t10\t30\tf1\nchr_sim\t40\t70\tf1\n')
    localization = score_localization(bed, truth, 100)[0]
    assert localization['true_assembly_bp'] == localization['predicted_assembly_bp'] == 50
    assert localization['base_recall'] == localization['base_precision'] == 1
    assert localization['truth_fragments'] == localization['fragments'] == 2

    copy = tmp_path/'copy.tsv'
    copy.write_text('family_id\testimated_copy_number\tcopy_number_interval_low\tcopy_number_interval_high\n'
                    'f1\t5\t4\t6\n')
    sampling = dict(sampled_repeat_bp={'f1':50}, actual_base_coverage=1)
    assert score_copy_number(copy, truth, sampling)[0]['truth_copies'] == 5

    comp = tmp_path/'comp.tsv'
    comp.write_text('family_id\tassembly_read_ratio\tstatus\nf1\t0.4\tpossible_collapse\n')
    retained = [dict(row, copies=1, repeat_bp=10, end=row['start']+10) for row in truth]
    scored = score_comparison(comp, truth, retained)[0]
    assert scored['truth_assembly_read_ratio'] == .4
    assert scored['outcome'] == 'TP'


def test_abundance_runner_small_cli_workflow(tmp_path: Path):
    config = dict(seeds={'development':[19], 'heldout':[20]}, periods=[31], copies=[10], flank_bp=300,
                  coverages=[5], read_length=100, substitution_rates=[0], assembly_fractions=[1,0],
                  collapse_threshold=.6, k=11, timeout_seconds=30)
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    run(config_path, tmp_path/'result', 'development')
    validation = json.loads((tmp_path/'result'/'validation.json').read_text())
    assert validation['successful'] == validation['executions'] == 5
    assert validation['copy_number_family_rows'] == 1
    assert validation['comparison_family_rows'] == 2
    assert json.loads((tmp_path/'result'/'environment.json').read_text())['benchmark_source_sha256']


def test_abundance_runner_crosses_divergence_and_fragmentation_factors(tmp_path: Path):
    config = dict(
        seeds={'development':[29], 'heldout':[30]}, periods=[31], copies=[10], flank_bp=300,
        unit_substitution_rates=[0, .1], array_fragment_counts=[1, 2], fragment_gap_bp=17,
        coverages=[5], read_length=100, substitution_rates=[0], assembly_fractions=[1, 0],
        collapse_threshold=.6, k=11, timeout_seconds=30,
    )
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    run(config_path, tmp_path/'result', 'development')
    validation = json.loads((tmp_path/'result'/'validation.json').read_text())
    assert validation['challenge_scenarios'] == 4
    assert validation['successful'] == validation['executions'] == 20
    assert validation['copy_number_family_rows'] == 4
    assert validation['localization_family_rows'] == 8
    assert validation['comparison_family_rows'] == 8
    rows = read_table(tmp_path/'result'/'comparison_metrics.tsv')
    assert {(row['unit_substitution_rate'], row['array_fragments']) for row in rows} == {
        ('0.0', '1'), ('0.0', '2'), ('0.1', '1'), ('0.1', '2'),
    }


def test_abundance_runner_localization_only_uses_selected_identity_model(tmp_path: Path):
    config = dict(
        seeds={'development':[31]}, periods=[31], copies=[10], flank_bp=300,
        unit_substitution_rates=[.05], array_fragment_counts=[1], fragment_gap_bp=0,
        coverages=[5], read_length=100, substitution_rates=[0], assembly_fractions=[1, 0],
        collapse_threshold=.6, k=11, timeout_seconds=30,
        locate_identity_model='iid_base', locate_min_identity=.9,
    )
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    result = tmp_path/'result'
    run(config_path, result, 'development', localization_only=True)
    validation = json.loads((result/'validation.json').read_text())
    assert validation['mode'] == 'localization_only'
    assert validation['successful'] == validation['executions'] == 2
    assert validation['localization_family_rows'] == 2
    assert validation['copy_number_family_rows'] == 0
    assert validation['comparison_family_rows'] == 0
    assert not (result/'reads').exists()
    environment = json.loads((result/'environment.json').read_text())
    assert environment['locate_identity_model'] == 'iid_base'
    assert environment['locate_min_identity'] == .9


def test_iid_heldout_requires_matching_development_evidence(tmp_path: Path):
    development = tmp_path/'development'
    development.mkdir()
    rows = (
        'assembly_fraction\tbase_recall\tbase_precision\tpredicted_assembly_bp\n'
        '1\t0.98\t0.99\t100\n'
        '0\t\t\t0\n'
    )
    (development/'localization_metrics.tsv').write_text(rows)
    development_config = {
        'seeds': {'development': [31]},
        'locate_identity_model': 'iid_base',
        'locate_min_identity': .9,
    }
    (development/'run_config.json').write_text(json.dumps(development_config))
    config_hash = digest_file(development/'run_config.json')
    (development/'validation.json').write_text(json.dumps({
        'complete': True, 'mode': 'localization_only', 'executions': 2,
        'successful': 2, 'localization_family_rows': 2,
    }))
    (development/'environment.json').write_text(json.dumps({
        'split': 'development', 'localization_only': True,
        'locate_identity_model': 'iid_base', 'locate_min_identity': .9,
        'config_sha256': config_hash,
    }))
    model = {
        'method': 'iid_base_exact_anchor_monomer_bridge',
        'identity_model': 'iid_base', 'k': 21, 'min_identity': .9,
        'array_merge_gap_rule': 'max(2*k,monomer_length)',
        'development_seeds': [31], 'development_result': str(development),
        'selection_gates': {
            'full_assembly_mean_base_recall_min': .95,
            'positive_assembly_mean_base_precision_min': .95,
            'absent_family_false_positive_rate_max': 0,
        },
        'observed_development_metrics': {
            'full_assembly_mean_base_recall': .98,
            'positive_assembly_mean_base_precision': .99,
            'absent_family_false_positive_rate': 0,
        },
    }
    files = {
        'development_validation_sha256': 'validation.json',
        'development_metrics_sha256': 'localization_metrics.tsv',
        'development_environment_sha256': 'environment.json',
        'development_config_sha256': 'run_config.json',
    }
    for field, filename in files.items():
        model[field] = digest_file(development/filename)
    config = {
        'seeds': {'development': [39], 'heldout': [40]}, 'k': 21,
        'locate_identity_model': 'iid_base', 'locate_min_identity': .9,
        'localizer_model': model,
    }
    assert validate_frozen_localizer(config, 'heldout') == {
        'full_assembly_mean_base_recall': .98,
        'positive_assembly_mean_base_precision': .99,
        'absent_family_false_positive_rate': 0,
    }
    (development/'localization_metrics.tsv').write_text(rows.replace('0.98', '0.97'))
    with pytest.raises(ValueError, match='evidence differs'):
        validate_frozen_localizer(config, 'heldout')


def test_classifier_heldout_requires_matching_robust_evidence(tmp_path: Path):
    development = tmp_path/'classifier_development'
    development.mkdir()
    development_config = {
        'development_seeds': [31], 'reserved_heldout_seeds': [40],
    }
    (development/'run_config.json').write_text(json.dumps(development_config))
    config_hash = digest_file(development/'run_config.json')
    metrics = {
        'minimum_seed_sensitivity_delta': .08,
        'maximum_seed_false_positive_rate_delta': 0,
        'minimum_seed_precision_delta': .01,
        'full_sensitivity_delta': .1,
        'full_false_positive_rate_delta': -.01,
        'full_precision_delta': .02,
    }
    (development/'validation.json').write_text(json.dumps({
        'complete': True, 'development_only': True,
        'development_seeds': [31], 'reserved_heldout_seeds': [40],
        'selected_candidate': 'blend_a0.5_t0.5',
        'acceptance': {**metrics, 'passed': True},
    }))
    (development/'environment.json').write_text(json.dumps({
        'development_seeds': [31], 'reserved_heldout_seeds': [40],
        'config_sha256': config_hash,
    }))
    (development/'selection.tsv').write_text(
        'candidate_id\tblend_alpha\tdecision_threshold\tselection_method\n'
        'blend_a0.5_t0.5\t0.5\t0.5\tseed_robust_minimax\n'
    )
    (development/'candidate_robust_summary.tsv').write_text('candidate_id\nblend_a0.5_t0.5\n')
    (development/'selected_metrics.tsv').write_text('seed\n31\n')
    model = {
        'method': 'log_space_single_multik_blend',
        'selection_method': 'seed_robust_minimax',
        'blend_alpha': .5, 'decision_threshold': .5,
        'k_values': [15, 21, 27, 31], 'fallback_rule': 'fallback_single_k21',
        'development_seeds': [31], 'development_result': str(development),
        'selection_gates': {
            'minimum_seed_sensitivity_delta': .05,
            'maximum_seed_false_positive_rate_delta': 0,
            'minimum_seed_precision_delta': 0,
            'full_sensitivity_delta': .05,
            'full_false_positive_rate_delta': 0,
            'full_precision_delta': 0,
        },
        'observed_development_metrics': metrics,
    }
    files = {
        'development_validation_sha256': 'validation.json',
        'development_selection_sha256': 'selection.tsv',
        'development_candidate_summary_sha256': 'candidate_robust_summary.tsv',
        'development_selected_metrics_sha256': 'selected_metrics.tsv',
        'development_environment_sha256': 'environment.json',
        'development_config_sha256': 'run_config.json',
    }
    for field, filename in files.items():
        model[field] = digest_file(development/filename)
    config = {
        'seeds': {'development': [31], 'heldout': [40]},
        'classifier_development_rule': {
            'method': 'log_space_single_multik_blend',
            'alpha_grid': [0, .25, .5, .75, 1],
            'decision_threshold_grid': [.45, .5, .55, .6],
            'multik_k_values': [15, 21, 27, 31],
            'unavailable_or_nonpositive_rule': 'fallback_single_k21',
        },
        'classifier_model': model,
    }
    assert validate_frozen_classifier(config, 'heldout') == metrics
    (development/'selected_metrics.tsv').write_text('seed\n32\n')
    with pytest.raises(ValueError, match='evidence differs'):
        validate_frozen_classifier(config, 'heldout')


def test_frozen_multik_evaluator_replays_every_challenge_scenario(tmp_path: Path):
    config = dict(
        seeds={'development':[39], 'heldout':[40]}, periods=[41], copies=[10], flank_bp=300,
        unit_substitution_rates=[0, .05], array_fragment_counts=[2], fragment_gap_bp=17,
        coverages=[5], read_length=100, substitution_rates=[0], assembly_fractions=[1, .5],
        collapse_threshold=.6, k=11, timeout_seconds=30,
        collapse_model={
            'method':'multik_loglinear_else_single_k21', 'k_values':[15, 21, 27, 31],
            'low_depth_cutoff':2.0, 'low_depth_threshold':.5,
            'standard_depth_threshold':.6, 'unavailable_rule':'fallback_single_k21',
            'calibration_seeds':[39], 'calibration_validation_sha256':'a'*64,
            'calibration_table_sha256':'b'*64, 'calibration_environment_sha256':'c'*64,
        },
    )
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    baseline = tmp_path/'baseline'
    run(config_path, baseline, 'heldout')
    evaluated = tmp_path/'evaluated'
    run_multik_collapse(baseline, evaluated, backend='python', split='heldout')
    validation = json.loads((evaluated/'validation.json').read_text())
    assert validation['complete'] is True
    assert validation['challenge_scenarios'] == 2
    assert validation['multik_read_conditions'] == 2
    assert validation['comparison_family_rows'] == 12
    assert validation['expected_comparison_family_rows'] == 12
