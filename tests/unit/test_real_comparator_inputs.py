import gzip

import pytest

from benchmarks.challenge.schema import ArrayRecord
from benchmarks.scripts.qc_complete_fastq import qc
from benchmarks.scripts.sample_complete_fastq import sample
from benchmarks.scripts.run_real_comparators import prepare_input, describe_arrays


def test_real_input_conversion_retains_identical_sequence_and_checks_receipt(tmp_path):
    path = tmp_path/'reads.fq.gz'
    path.write_bytes(gzip.compress(b'@r1 description\nACGT\n+\nIIII\n@r2\nNNAC\n+\n!!!!\n'))
    qc(path, tmp_path/'qc')
    sample(path, tmp_path/'qc/qc.json', tmp_path/'samples', ['1'], 'ERR123', 6101)
    receipt = tmp_path/'samples/sampling_receipt.json'
    row, lengths = prepare_input(receipt, 'sample_001', tmp_path/'out.fa')
    assert (tmp_path/'out.fa').read_text() == '>r1 description\nACGT\n>r2\nNNAC\n'
    assert lengths == {'r1': 4, 'r2': 4} and row['total_bases'] == 8
    changed = tmp_path/'samples/sample_001.fastq.gz'
    changed.write_bytes(gzip.compress(b'@r1\nACGT\n+\nIIII\n'))
    with pytest.raises(ValueError, match='receipt'):
        prepare_input(receipt, 'sample_001', tmp_path/'changed.fa')


def test_descriptive_union_is_not_duplicate_base_accuracy():
    arrays = [ArrayRecord('a', 2, 8, 3), ArrayRecord('a', 4, 10, 3), ArrayRecord('b', 0, 2, 1)]
    result = describe_arrays(arrays, {'a': 12, 'b': 8})
    assert result['observed_in_scope_calls'] == 3 and result['observed_positive_reads'] == 2
    assert result['observed_union_bp'] == 10 and result['observed_union_base_fraction'] == .5
    assert describe_arrays([], {'a': 12})['observed_union_bp'] == 0
    with pytest.raises(ValueError, match='coordinates'):
        describe_arrays([ArrayRecord('a', 0, 13, 3)], {'a': 12})
    with pytest.raises(ValueError, match='unknown'):
        describe_arrays([ArrayRecord('c', 0, 2, 1)], {'a': 12})
