"""Protocol and bookkeeping checks without consuming formal seeds."""
import random

import pytest

from benchmarks.challenge.schema import ArrayRecord
from benchmarks.challenge.simulate import Scenario
from benchmarks.challenge.unified_simulate import shared_fragment_read
from benchmarks.challenge.unified_correspondence import CatalogueCorrespondence
from benchmarks.scripts.run_srf_unified_comparison import abundance_metrics, conditions


def test_conditions_match_frozen_design():
    items = {c.name: c for c in conditions()}
    assert len(items) == 6
    assert items['clean'].positive_fraction == .7
    assert items['low_abundance'].positive_fraction == .1
    assert items['indel'].insertion_rate == items['indel'].deletion_rate == .001
    assert items['shared_fragment'].period == 120
    for c in items.values():
        if c.negative_kind != "shared_fragment":
            c.validate()


def test_background_partial_fragment_deterministic_and_length():
    c = Scenario('test', period=120, negative_kind='shared_fragment')
    monomer = 'A' * 100 + 'C' * 20
    seq = shared_fragment_read(c.read_length, monomer, random.Random(17))
    assert len(seq) == 5000
    assert monomer not in seq
    assert seq.count('A' * 100) > 20
    assert seq == shared_fragment_read(c.read_length, monomer, random.Random(17))
    with pytest.raises(ValueError, match='period >100'):
        shared_fragment_read(5000, 'A'*100, random.Random(17))


def test_missing_family_and_unassigned_mass_remain():
    correspondence = {'n1': CatalogueCorrespondence('unique', ('f1',), (('f1', 1.0),)),
                      'n2': CatalogueCorrespondence('unmatched', (), ())}
    truth = [ArrayRecord('r1', 0, 100, 10, 'A'*10, 'f1'),
             ArrayRecord('r2', 0, 100, 10, 'C'*10, 'f2')]
    got = abundance_metrics({'n1': 100, 'n2': 77}, correspondence, truth, {'f1': 'A'*10, 'f2': 'C'*10})
    assert got['MARE'] == .5
    assert got['unassigned_abundance_bp'] == 77
    assert got['per_family'][1]['estimated_bp'] == 0


def test_overlay_preserves_positive_arrays_and_truth_hashes(tmp_path):
    from dataclasses import replace
    from benchmarks.challenge.unified_simulate import generate_unified_dataset
    from benchmarks.challenge.simulate import generate_dataset
    from benchmarks.challenge.adapters import read_fasta
    from benchmarks.challenge.schema import digest_file, read_table
    config = Scenario('smoke', period=120, negative_kind='shared_fragment', read_count=10)
    original = tmp_path / 'original'
    overlay = tmp_path / 'overlay'
    generate_dataset(replace(config, negative_kind='random'), 17, original)
    manifest = generate_unified_dataset(config, 17, overlay)
    assert digest_file(original / 'truth_arrays.tsv') == digest_file(overlay / 'truth_arrays.tsv')
    before, after = read_fasta(original / 'reads.fa'), read_fasta(overlay / 'reads.fa')
    for r in read_table(overlay / 'truth_reads.tsv'):
        if int(r['truth_positive']):
            assert before[r['read_id']] == after[r['read_id']]
        else:
            assert r['negative_kind'] == 'shared_fragment'
            assert before[r['read_id']] != after[r['read_id']]
    for name, item in manifest['files'].items():
        assert digest_file(overlay / name) == item['sha256']
