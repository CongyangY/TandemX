from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path

import pytest

from benchmarks.scripts.fetch_ena_complete import files_for_run, download, verify_file
from benchmarks.scripts.qc_complete_fastq import qc, n50


def test_complete_qc_distribution_and_denominators(tmp_path: Path):
    path=tmp_path/'reads.fastq.gz'
    content=b'@r1\nACGT\n+\nIIII\n@r2\nNN\n+\n!!\n@r3\nCG\n+\n55\n'
    with gzip.open(path,'wb') as handle: handle.write(content)
    result=qc(path,tmp_path/'qc',3,8)
    assert result['read_n50']==4 and result['median_length']==2
    assert n50(Counter({2:2,4:1}))==4
    assert result['gc_fraction']==.5 and result['n_fraction']==.25
    assert result['mean_reported_error_probability']==pytest.approx((4*.0001+2+2*.01)/8)
    assert result['gzip_trailer_checked']
    assert (tmp_path/'qc'/'read_ids.sqlite').exists()


@pytest.mark.parametrize('data,error',[
    (b'', 'Empty'),
    (b'@r\nAC\n+\nII\n@r\nAC\n+\nII\n','Duplicate'),
    (b'@r\nAC\n+\nI\n','Invalid'),
    (b'@r\nAC\n+r2\nII\n','identifier'),
    (b'@r\nAC\n+\n','Incomplete'),
    (b'@r\nAZ\n+\nII\n','Invalid'),
])
def test_bad_fastq_retains_failure_receipt(tmp_path: Path,data: bytes,error: str):
    source=tmp_path/'reads.fq';source.write_bytes(data)
    with pytest.raises(ValueError,match=error):qc(source,tmp_path/'qc')
    assert not json.loads((tmp_path/'qc'/'qc.json').read_text())['complete']


def test_truncated_gzip_is_not_complete(tmp_path: Path):
    path=tmp_path/'reads.fq.gz'
    path.write_bytes(gzip.compress(b'@r\nACGT\n+\nIIII\n')[:-6])
    with pytest.raises(EOFError): qc(path,tmp_path/'qc')
    assert not json.loads((tmp_path/'qc'/'qc.json').read_text())['complete']


def test_download_metadata_and_resume(tmp_path: Path,monkeypatch):
    payload=gzip.compress(b'@r\nACGT\n+\nIIII\n')
    md5=hashlib.md5(payload).hexdigest()
    meta=tmp_path/'ena.tsv'
    meta.write_text('run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n'
                    f'ERR123\tftp.sra.ebi.ac.uk/a/ERR123.fastq.gz\t{md5}\t{len(payload)}\n')
    _, entries=files_for_run(meta,'ERR123')
    entry=entries[0]
    partial=tmp_path/(entry['filename']+'.partial');partial.write_bytes(payload[:10])
    class Response(io.BytesIO):
        status=206
        headers={'Content-Range':f'bytes 10-{len(payload)-1}/{len(payload)}'}
    def open_response(request,timeout):
        assert request.get_header('Range')=='bytes=10-'
        return Response(payload[10:])
    monkeypatch.setattr('urllib.request.urlopen',open_response)
    assert download(entry,tmp_path)['source_md5_verified']
    assert (tmp_path/entry['filename']).read_bytes()==payload
    assert download(entry,tmp_path)['transfer']=='existing_reverified'
    with pytest.raises(ValueError,match='MD5'):verify_file(tmp_path/entry['filename'],len(payload),'0'*32)
    meta.write_text(meta.read_text().replace('ftp.sra.ebi.ac.uk','example.com'))
    with pytest.raises(ValueError):files_for_run(meta,'ERR123')


def test_count_mismatch_and_ignored_resume_are_explicit(tmp_path: Path,monkeypatch):
    path=tmp_path/'reads.fq';path.write_bytes(b'@r\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError,match='Read count'):qc(path,tmp_path/'qc',2,4)
    partial=tmp_path/'a.fq.gz.partial';partial.write_bytes(b'ab')
    class Response(io.BytesIO):
        status=200
        headers={}
    monkeypatch.setattr('urllib.request.urlopen',lambda *args,**kwargs:Response(b'abcd'))
    with pytest.raises(ValueError,match='ignored resume'):
        download(dict(filename='a.fq.gz',url='https://ftp.sra.ebi.ac.uk/a.fq.gz',expected_bytes=4,expected_md5='0'*32),tmp_path)
    assert partial.read_bytes()==b'ab'
