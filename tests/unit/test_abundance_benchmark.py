"""Check truth generation and score semantics without invoking detector internals."""
import json
from pathlib import Path

import pytest

from benchmarks.abundance.evaluate import finite_nonnegative, score_copy_number, score_localization, score_comparison
from benchmarks.abundance.simulate import GenomeSpec, build_genome, sample_reads
from benchmarks.abundance.run import aggregate, run
from benchmarks.challenge.schema import read_table


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
