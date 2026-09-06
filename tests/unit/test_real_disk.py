import gzip
import json
import random
import sqlite3

import pytest

from benchmarks.challenge.schema import ArrayRecord, digest_file, iter_table, write_table
from benchmarks.challenge.adapters import iter_arrays, parse_arrays
from benchmarks.scripts.real_disk import ARRAY_FIELDS, normalize_disk, prepare_disk_input
from benchmarks.scripts.run_real_comparators import describe_arrays, prepare_input


def selected(tmp_path, entries):
    path = tmp_path/'sample_001.fastq.gz'
    with gzip.open(path, 'wt') as handle:
        reads = bases = 0
        for name, seq in entries:
            handle.write(f'@{name}\n{seq}\n+\n'+('I'*len(seq))+'\n')
            reads += 1
            bases += len(seq)
    receipt = tmp_path/'sampling_receipt.json'
    receipt.write_text(json.dumps(dict(complete=True, samples=[dict(
        sample_id='sample_001', status='ok', read_count=reads, total_bases=bases,
        fastq_sha256=digest_file(path))])))
    return receipt


def test_disk_input_parity_and_refusal_to_overwrite(tmp_path):
    receipt = selected(tmp_path, [('r1 description', 'ACGT'), ('r2', 'NNAC')])
    old, lengths = prepare_input(receipt, 'sample_001', tmp_path/'memory.fa')
    new = prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'index.sqlite')
    assert old == new and (tmp_path/'memory.fa').read_bytes() == (tmp_path/'disk.fa').read_bytes()
    with sqlite3.connect(tmp_path/'index.sqlite') as db:
        assert dict(db.execute('SELECT * FROM reads')) == lengths
    with pytest.raises(ValueError, match='overwrite'):
        prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'new.sqlite')


def test_disk_accepts_more_than_old_read_cap(tmp_path):
    # Generated stress input, not a checked-in toy fixture. Single-base reads
    # independently exercise the controller cap without a multi-GB test dataset.
    receipt = selected(tmp_path, ((f'r{i}', 'A') for i in range(100_001)))
    with pytest.raises(ValueError, match='cap'):
        prepare_input(receipt, 'sample_001', tmp_path/'memory.fa')
    result = prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'index.sqlite')
    assert result['read_count'] == result['total_bases'] == 100_001


@pytest.mark.parametrize('fault', ['duplicate', 'hash', 'truncated', 'empty'])
def test_disk_invalid_input_has_no_published_fasta(tmp_path, fault):
    entries = [('r1', 'ACGT'), ('r1' if fault == 'duplicate' else 'r2', 'NNAC')]
    receipt = selected(tmp_path, [] if fault == 'empty' else entries)
    path = tmp_path/'sample_001.fastq.gz'
    if fault == 'hash':
        metadata = json.loads(receipt.read_text())
        metadata['samples'][0]['fastq_sha256'] = '0'*64
        receipt.write_text(json.dumps(metadata))
    if fault == 'truncated':
        path.write_bytes(path.read_bytes()[:-5])
    with pytest.raises((ValueError, EOFError)):
        prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'index.sqlite')
    assert not (tmp_path/'disk.fa').exists()
    with sqlite3.connect(tmp_path/'index.sqlite') as db:
        assert db.execute("SELECT value FROM metadata WHERE name='ready'").fetchone() is None


def native_file(path, tool, arrays):
    if tool == 'tandemx':
        write_table(path, [dict(read_id=a.read_id, read_start=a.start, read_end=a.end, period_bp=a.period)
                          for a in arrays], ['read_id', 'read_start', 'read_end', 'period_bp'])
    elif tool == 'ultra':
        write_table(path, [dict(SeqID=a.read_id, Start=a.start, End=a.end, Period=a.period, Consensus='ACGT')
                          for a in arrays], ['SeqID', 'Start', 'End', 'Period', 'Consensus'])
    elif tool == 'trf':
        path.write_text(''.join(f'@{a.read_id}\n{a.start+1} {a.end} {a.period} 2 10 90 0 100 25 25 25 25 2 ACGT\n'
                                for a in arrays))
    else:
        path.write_text(''.join(f'{a.read_id} 1 500 2 {a.start+1} {a.end} {a.period} 2 0 0 ACGT\n' for a in arrays))


@pytest.mark.parametrize('tool', ['tandemx', 'trf', 'tidehunter', 'ultra'])
def test_disk_union_and_native_order_match_independent_base_mask(tmp_path, tool):
    rng = random.Random(431)
    lengths = {f'r{i}': 500 for i in range(23)}
    receipt = selected(tmp_path, ((name, 'A'*length) for name, length in lengths.items()))
    prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'index.sqlite')
    arrays = []
    for _ in range(170):
        start = rng.randrange(300)
        arrays.append(ArrayRecord(rng.choice(list(lengths)), start, rng.randrange(start+1, 501), rng.choice([20, 100, 1001])))
    arrays.extend([ArrayRecord('r0', 0, 500, 100)]*3)
    native = tmp_path/'native.txt'
    native_file(native, tool, arrays)
    parsed = parse_arrays(tool, native, 30, 1000, 100)
    assert list(iter_arrays(tool, native, 30, 1000, 100)) == parsed
    result = normalize_disk(tool, native, tmp_path/'index.sqlite', tmp_path/'normalized.tsv')
    mask = {(a.read_id, b) for a in parsed for b in range(a.start, a.end)}
    assert result == describe_arrays(parsed, lengths)
    assert result['observed_union_bp'] == len(mask)
    from dataclasses import asdict
    write_table(tmp_path/'expected.tsv', (asdict(a) for a in parsed), ARRAY_FIELDS)
    assert (tmp_path/'normalized.tsv').read_bytes() == (tmp_path/'expected.tsv').read_bytes()


@pytest.mark.parametrize('fault', ['unknown', 'past_end', 'late_malformed', 'no_calls'])
def test_disk_normalization_failure_is_not_partial_success(tmp_path, fault):
    receipt = selected(tmp_path, [('r1', 'ACGT'*100)])
    prepare_disk_input(receipt, 'sample_001', tmp_path/'disk.fa', tmp_path/'index.sqlite')
    arrays = [ArrayRecord('r1', 0, 200, 100)]
    if fault == 'unknown':
        arrays.append(ArrayRecord('missing', 0, 200, 100))
    if fault == 'past_end':
        arrays.append(ArrayRecord('r1', 0, 401, 100))
    native = tmp_path/'native.tsv'
    native_file(native, 'tandemx', [] if fault == 'no_calls' else arrays)
    if fault == 'late_malformed':
        with native.open('a') as f:
            f.write('r1\t0\t200\n')
    if fault == 'no_calls':
        result = normalize_disk('tandemx', native, tmp_path/'index.sqlite', tmp_path/'out.tsv')
        assert result['observed_in_scope_calls'] == result['observed_union_bp'] == 0
    else:
        with pytest.raises(ValueError):
            normalize_disk('tandemx', native, tmp_path/'index.sqlite', tmp_path/'out.tsv')
        assert not (tmp_path/'out.tsv').exists()
        with sqlite3.connect(tmp_path/'index.sqlite') as db:
            assert db.execute('SELECT count(*) FROM arrays').fetchone()[0] == 0


def test_tsv_iterator_is_lazy_and_checks_late_rows(tmp_path):
    path = tmp_path/'table.tsv'
    path.write_text('a\tb\n1\t2\n3\n')
    rows = iter_table(path, {'a', 'b'})
    assert next(rows) == {'a': '1', 'b': '2'}
    with pytest.raises(ValueError, match='Malformed'):
        next(rows)


def test_replay_controller_checks_all_tools_and_input_hashes(tmp_path):
    from dataclasses import asdict
    from benchmarks.scripts.replay_real_normalization import replay
    receipt = selected(tmp_path, [('r1', 'ACGT'*100), ('r2', 'TGCA'*100)])
    previous = tmp_path/'previous'
    previous.mkdir()
    metadata, lengths = prepare_input(receipt, 'sample_001', previous/'reads.fa')
    (previous/'environment.json').write_text(json.dumps({'input': metadata}))
    summaries = []
    paths = {'tandemx':'discover/candidate_reads.tsv', 'trf':'trf.txt', 'tidehunter':'tidehunter.tsv'}
    for tool, name in paths.items():
        path = previous/tool/name
        path.parent.mkdir(parents=True)
        native_file(path, tool, [ArrayRecord('r1', 0, 200, 100), ArrayRecord('r1', 100, 300, 100)])
        arrays = parse_arrays(tool, path, 30, 1000, 100)
        write_table(previous/tool/'normalized_arrays.tsv', (asdict(a) for a in arrays), ARRAY_FIELDS)
        summaries.append(dict(tool=tool, normalization='ok', **describe_arrays(arrays, lengths)))
    write_table(previous/'summary.tsv', summaries, list(summaries[0]))
    result = replay(receipt, 'sample_001', tmp_path/'replay', previous)
    assert result['complete'] and len(result['normalized_tools']) == 3
    assert all(r['byte_identical'] and r['metrics_identical'] for r in result['normalized_tools'])
    (previous/'reads.fa').write_text('changed')
    with pytest.raises(ValueError, match='FASTA'):
        replay(receipt, 'sample_001', tmp_path/'invalid', previous)
    failure = json.loads((tmp_path/'invalid/replay_receipt.json').read_text())
    assert not failure['complete'] and failure['error']
