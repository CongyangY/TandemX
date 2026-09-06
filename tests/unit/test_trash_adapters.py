import csv

import pytest

from benchmarks.challenge.trash_adapters import ARRAY_FILES, UNIT_FILES, iter_regions, iter_units, iter_unit_extents
from benchmarks.challenge.trash_adapters import region_catalogue


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_array_window_and_extraction_interpretations_are_separate(tmp_path):
    write_csv(tmp_path/ARRAY_FILES['trash'], [dict(start=0, end=99, **{
        'fasta.name': 'chrA', 'most.freq.value.N': 4, 'consensus.primary': 'acg'})])
    grid = list(iter_regions(tmp_path, 'trash', 'window_grid'))[0]
    extract = list(iter_regions(tmp_path, 'trash', 'r_extraction', 'consensus_length'))[0]
    assert (grid.start, grid.end, grid.period, grid.sequence) == (0, 100, 4, 'ACG')
    assert (extract.start, extract.end, extract.period) == (0, 99, 3)
    with pytest.raises(ValueError, match='TRASH1'):
        list(iter_regions(tmp_path, 'trash', 'one_based'))


def test_trash2_arrays_and_units_include_last_base_and_strand(tmp_path):
    write_csv(tmp_path/ARRAY_FILES['trash2'], [dict(start=1, end=20, seqID='chrB', top_N=5,
                                                 representative='atgcn')])
    write_csv(tmp_path/UNIT_FILES['trash2'], [dict(start=16, end=20, seqID='chrB', width=5, strand='-')])
    region = list(iter_regions(tmp_path, 'trash2', 'one_based'))[0]
    unit = list(iter_units(tmp_path, 'trash2'))[0]
    assert (region.start, region.end, region.period) == (0, 20, 5)
    assert (unit.start, unit.end) == (15, 20)
    with pytest.raises(ValueError, match='policy'):
        list(iter_units(tmp_path, 'trash2', -1))


def test_trash1_offsets_are_explicit_and_never_clip_invalid_units(tmp_path):
    path = tmp_path/UNIT_FILES['trash']
    write_csv(path, [dict(start=2, end=6, width=5, strand='+', **{'seq.name': 'chrA'})])
    assert list(iter_units(tmp_path, 'trash', -1))[0].start == 0
    assert list(iter_units(tmp_path, 'trash', 0))[0].start == 1
    write_csv(path, [dict(start=1, end=5, width=5, strand='+', **{'seq.name': 'chrA'})])
    with pytest.raises(ValueError, match='Invalid array record'):
        list(iter_units(tmp_path, 'trash', -1))


def test_missing_consensus_and_late_malformed_rows_are_not_omitted(tmp_path):
    path = tmp_path/ARRAY_FILES['trash2']
    row = dict(start=1, end=20, seqID='chrB', top_N=5, representative='NA')
    write_csv(path, [row])
    with pytest.raises(ValueError, match='Unresolved'):
        list(iter_regions(tmp_path, 'trash2', 'one_based'))
    write_csv(path, [{**row, 'representative': 'ACGTA'}])
    with path.open('a') as handle:
        handle.write('1,20,chrB\n')
    with pytest.raises(ValueError, match='Malformed'):
        list(iter_regions(tmp_path, 'trash2', 'one_based'))


def test_unit_extent_uses_native_membership_and_preserves_internal_gaps(tmp_path):
    write_csv(tmp_path/ARRAY_FILES['trash2'], [dict(start=1, end=100, seqID='chrB', top_N=5,
                                                 representative='ACGTA', array_num_ID=8)])
    units = [dict(start=20, end=24, seqID='chrB', width=5, strand='-', arrayID=8),
             dict(start=30, end=34, seqID='chrB', width=5, strand='+', arrayID=8)]
    write_csv(tmp_path/UNIT_FILES['trash2'], units)
    region = list(iter_unit_extents(tmp_path))[0]
    assert (region.start, region.end, region.period) == (19, 34, 5)
    write_csv(tmp_path/UNIT_FILES['trash2'], [{**units[0], 'arrayID': 9}])
    with pytest.raises(ValueError, match='unknown array'):
        list(iter_unit_extents(tmp_path))
    # An empty main repeat table cannot make a native array disappear.
    with (tmp_path/UNIT_FILES['trash2']).open('w') as handle:
        handle.write(','.join(units[0])+'\n')
    with pytest.raises(ValueError, match='without observed'):
        list(iter_unit_extents(tmp_path))


def test_fractional_native_period_is_preserved_not_rounded_or_dropped(tmp_path):
    from benchmarks.challenge.evaluate import score_arrays
    from benchmarks.challenge.schema import ArrayRecord
    row = dict(start=0, end=999, **{'fasta.name': 'chrA', 'most.freq.value.N': '171.5',
                                   'consensus.primary': 'ACGTA'})
    path = tmp_path/ARRAY_FILES['trash']
    write_csv(path, [row])
    native = list(iter_regions(tmp_path, 'trash', 'window_grid'))
    assert native[0].period == 171.5
    metrics, _ = score_arrays(native, [ArrayRecord('chrA', 0, 1000, 171)], {'chrA': 1000})
    assert metrics['matched_period_mae_bp'] == .5
    for invalid in ['nan', 'inf', '-1', '0']:
        write_csv(path, [{**row, 'most.freq.value.N': invalid}])
        with pytest.raises(ValueError, match='period'):
            list(iter_regions(tmp_path, 'trash', 'window_grid'))


def test_secondary_consensus_is_retained_with_explicit_missing_and_scope_counts(tmp_path):
    rows = [{'consensus.primary': 'ACGTA'*7, 'consensus.secondary': 'GCTAA'*8},
            {'consensus.primary': 'ACG', 'consensus.secondary': 'NA'}]
    path = tmp_path/ARRAY_FILES['trash']
    write_csv(path, rows)
    seq, counts = region_catalogue(tmp_path, 'trash', True)
    assert seq == ['ACGTA'*7, 'GCTAA'*8]
    assert counts == dict(native_region_count=2, secondary_missing_regions=1,
                          excluded_sequence_length_count=1, in_scope_native_sequences=2)
    assert region_catalogue(tmp_path, 'trash')[0] == ['ACGTA'*7]
    write_csv(path, [{**rows[0], 'consensus.secondary': 'ACRY'}])
    with pytest.raises(ValueError, match='DNA'):
        region_catalogue(tmp_path, 'trash', True)
