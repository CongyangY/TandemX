import gzip
import hashlib

import pytest

from benchmarks.scripts.qc_reference_fasta import reference_qc


def test_reference_qc_gzip_streaming_and_iupac_counts(tmp_path):
    source = tmp_path/'ref.fa.gz'
    source.write_bytes(gzip.compress(b'>chr1 description\nacgt\nNN\n>ChrM\nRY\n'))
    result = reference_qc(source)
    assert result['total_bases'] == 8 and result['contig_count'] == 2
    assert result['n_bases'] == 2 and result['other_ambiguous_bases'] == 2
    assert result['contigs'][0]['length_bp'] == 6
    assert result['input_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    source.write_bytes(source.read_bytes()[:-4])
    with pytest.raises(EOFError):
        reference_qc(source)


def test_reference_qc_streams_unwrapped_chromosome_lines(tmp_path):
    source = tmp_path / 'unwrapped.fa'
    sequence = b'ACGT' * 300_000
    source.write_bytes(b'>Chr1\n' + sequence + b'\n>Chr2\nNN\n')
    result = reference_qc(source)
    assert result['contig_count'] == 2
    assert result['total_bases'] == len(sequence) + 2
    assert result['contigs'][0]['length_bp'] == len(sequence)


def test_reference_qc_still_bounds_headers(tmp_path):
    source = tmp_path / 'long_header.fa'
    source.write_bytes(b'>' + b'x' * (1024 * 1024) + b'\nAC\n')
    with pytest.raises(ValueError, match='header exceeds'):
        reference_qc(source)


@pytest.mark.parametrize('content', [b'', b'>r\n', b'>r\n>r2\nAC\n', b'>r\nAC\n>r\nAC\n',
                                     b'>\nAC\n', b'AC\n', b'>r\nACX\n'])
def test_reference_invalid_inputs(tmp_path, content):
    source = tmp_path/'ref.fa'
    source.write_bytes(content)
    with pytest.raises(ValueError):
        reference_qc(source)
