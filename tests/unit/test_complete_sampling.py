import gzip
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.fastq_stream import hashed_fastq, records
from benchmarks.scripts.qc_complete_fastq import qc
from benchmarks.scripts.sample_complete_fastq import fraction_threshold, read_hash, sample


def make_reads(path: Path, reverse: bool = False):
    rows = [f'@r{i} description\n'+('ACGT'*(i%20+1))+f'\n+r{i}\n'+('I'*(4*(i%20+1)))+'\n'
            for i in range(1000)]
    path.write_bytes(gzip.compress(''.join(reversed(rows) if reverse else rows).encode()))


def read_ids(path: Path) -> set[bytes]:
    with hashed_fastq(path) as (handle, digest):
        identifiers = {r.identifier for r in records(handle)}
    assert digest.hexdigest() == hashlib.sha256(path.read_bytes()).hexdigest()
    return identifiers


def test_nested_order_independent_membership_and_reproducible_bytes(tmp_path: Path):
    path = tmp_path/'reads.fq.gz'
    make_reads(path)
    verified = qc(path, tmp_path/'qc')
    receipt = tmp_path/'qc/qc.json'
    first = sample(path, receipt, tmp_path/'a', ['1', '.1', '.5'], 'ERR123', 123, 1000)
    second = sample(path, receipt, tmp_path/'b', ['1', '.1', '.5'], 'ERR123', 123, 1000)
    assert first == second
    sets = [read_ids(tmp_path/'a'/f'sample_{i:03d}.fastq.gz') for i in range(1, 4)]
    assert sets[0] < sets[1] < sets[2] and len(sets[2]) == 1000
    assert 400 < len(sets[1]) < 600
    assert first['samples'][2]['total_bases'] == verified['total_bases']
    assert first['samples'][2]['read_n50'] == verified['read_n50']
    assert first['samples'][2]['gc_fraction'] == .5
    assert first['samples'][2]['mean_reported_error_probability'] == pytest.approx(.0001)
    assert first['samples'][2]['nominal_total_base_coverage'] == verified['total_bases']/1000
    for i in range(1, 4):
        file = f'sample_{i:03d}.fastq.gz'
        assert (tmp_path/'a'/file).read_bytes() == (tmp_path/'b'/file).read_bytes()
    other = tmp_path/'reverse.fq.gz'
    make_reads(other, True)
    qc(other, tmp_path/'qc2')
    sample(other, tmp_path/'qc2/qc.json', tmp_path/'c', ['.1', '.5', '1'], 'ERR123', 123)
    assert [read_ids(tmp_path/'c'/f'sample_{i:03d}.fastq.gz') for i in range(1, 4)] == sets
    # Seeds and libraries change membership; identical sequences do not force co-selection.
    assert read_hash(b'r1', 'ERR123', 123) != read_hash(b'r1', 'ERR123', 124)
    assert read_hash(b'r1', 'ERR123', 123) != read_hash(b'r1', 'ERR124', 123)


def test_exact_threshold_and_invalid_sampling_requests(tmp_path: Path):
    assert fraction_threshold('.5')[1] == 1 << 127
    assert fraction_threshold('1')[1] == 1 << 128
    assert fraction_threshold('.1')[1] == (1 << 128)//10
    for invalid in ['0', '-.1', 'nan', 'Infinity', '1.1', 'garbage', '1e-100']:
        with pytest.raises(ValueError):
            fraction_threshold(invalid)
    path = tmp_path/'reads.fq.gz'
    make_reads(path)
    qc(path, tmp_path/'qc')
    receipt = tmp_path/'qc/qc.json'
    for fractions in [[], ['.5', '.50'], [str(i/100) for i in range(1, 14)]]:
        with pytest.raises(ValueError, match='distinct'):
            sample(path, receipt, tmp_path/'bad', fractions, 'ERR123', 0)
    for namespace, seed in [('wrong namespace', 1), ('ERR123', -1), ('ERR123', 1 << 64)]:
        with pytest.raises(ValueError):
            sample(path, receipt, tmp_path/'bad', ['.1'], namespace, seed)
    result = sample(path, receipt, tmp_path/'empty', ['1e-30'], 'ERR123', 0)
    assert result['complete'] and result['samples'][0]['status'] == 'empty_random_sample'
    assert result['samples'][0]['read_n50'] is None
    assert read_ids(tmp_path/'empty/sample_001.fastq.gz') == set()


def test_stale_qc_never_accepts_output(tmp_path: Path):
    path = tmp_path/'reads.fq.gz'
    make_reads(path)
    qc(path, tmp_path/'qc')
    receipt = tmp_path/'qc/qc.json'
    path.write_bytes(gzip.compress(b'@different\nAC\n+\nII\n'))
    with pytest.raises(ValueError, match='SHA256'):
        sample(path, receipt, tmp_path/'bad', ['1'], 'ERR123', 0)
    assert not (tmp_path/'bad/sample_001.fastq.gz').exists()
    assert (tmp_path/'bad/sample_001.fastq.gz.partial').exists()
    assert not json.loads((tmp_path/'bad/sampling_receipt.json').read_text())['complete']
    failed = tmp_path/'failed.json'
    failed.write_text('{"complete":false}')
    with pytest.raises(ValueError, match='successful'):
        sample(path, failed, tmp_path/'bad2', ['1'], 'ERR123', 0)


def test_hashed_parser_all_gzip_members_and_plain_crlf(tmp_path: Path):
    first = b'@r1\r\nAC\r\n+\r\nII\r\n'
    second = b'@r2\nGT\n+\nII'  # Valid final quality line need not have a newline.
    path = tmp_path/'all.fq.gz'
    path.write_bytes(gzip.compress(first)+gzip.compress(second))
    assert read_ids(path) == {b'r1', b'r2'}
    plain = tmp_path/'plain.fq'
    plain.write_bytes(first+second)
    assert read_ids(plain) == {b'r1', b'r2'}
    result = qc(path, tmp_path/'qc')
    assert result['input_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    path.write_bytes(path.read_bytes()[:-3])
    with pytest.raises(EOFError):
        read_ids(path)
