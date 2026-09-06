from io import BytesIO
import hashlib
from urllib.error import HTTPError

import pytest

from benchmarks.scripts import fetch_morex_reference as module


class Response(BytesIO):
    def __init__(self, content, status=200, headers=None):
        super().__init__(content)
        self.status = status
        self.headers = headers or {}


def test_sha_download_and_exact_resume(tmp_path, monkeypatch):
    content = b'>chr1\nACGTN\n'
    sha = hashlib.sha256(content).hexdigest()
    path = tmp_path/'reference.fa'
    path.with_suffix('.fa.partial').write_bytes(content[:5])
    def response(request, timeout):
        assert request.get_header('Range') == 'bytes=5-'
        return Response(content[5:], 206, {'Content-Range': 'bytes 5-11/12'})
    monkeypatch.setattr(module, 'urlopen', response)
    result = module.download_sha256('https://example.org/reference', path, sha, 100)
    assert result['bytes'] == len(content) and path.read_bytes() == content
    assert module.download_sha256('https://example.org/reference', path, sha, 100)['transfer'] == 'existing_reverified'


@pytest.mark.parametrize('fault', ['wrong_range', 'checksum', 'budget', 'forbidden'])
def test_sha_download_failures_preserve_partial_and_do_not_publish(tmp_path, monkeypatch, fault):
    content = b'ACGTN'
    path = tmp_path/'reference.fa'
    if 'range' in fault:
        path.with_suffix('.fa.partial').write_bytes(b'A')
    calls = []
    def response(request, timeout):
        calls.append(request)
        if fault == 'forbidden':
            raise HTTPError(request.full_url, 403, 'Forbidden', {}, None)
        if fault == 'wrong_range':
            return Response(content, 206, {'Content-Range': 'bytes 0-4/5'})
        return Response(content)
    monkeypatch.setattr(module, 'urlopen', response)
    sha = '0'*64 if fault == 'checksum' else hashlib.sha256(content).hexdigest()
    with pytest.raises((ValueError, HTTPError)):
        module.download_sha256('https://example.org/reference', path, sha, 4 if fault == 'budget' else 100)
    assert not path.exists() and len(calls) == 1


def test_ignored_range_uses_only_a_verified_full_response_prefix(tmp_path, monkeypatch):
    content = b'>chr1\nACGTN\n'
    path = tmp_path/'reference.fa'
    partial = path.with_suffix('.fa.partial')
    partial.write_bytes(content[:7])
    monkeypatch.setattr(module, 'urlopen', lambda request, timeout: Response(content, 200))
    result = module.download_sha256(
        'https://example.org/reference', path, hashlib.sha256(content).hexdigest(), 100
    )
    assert result['transfer'] == 'downloaded_complete_after_verified_prefix'
    assert path.read_bytes() == content and not partial.exists()

    path.unlink()
    partial.write_bytes(content[:7])
    corrupted = b'X' + content[1:]
    monkeypatch.setattr(module, 'urlopen', lambda request, timeout: Response(corrupted, 200))
    with pytest.raises(ValueError, match='prefix differs'):
        module.download_sha256(
            'https://example.org/reference', path, hashlib.sha256(content).hexdigest(), 100
        )
    assert partial.read_bytes() == content[:7] and not path.exists()


def test_complete_partial_is_reverified_without_network(tmp_path, monkeypatch):
    path = tmp_path/'reference.fa'
    path.with_suffix('.fa.partial').write_bytes(b'ACGT')
    monkeypatch.setattr(module, 'urlopen', lambda *a, **kw: pytest.fail('Unexpected network'))
    result = module.download_sha256('https://example.org/ref', path, hashlib.sha256(b'ACGT').hexdigest(), 4)
    assert result['transfer'] == 'complete_partial_reverified'


def test_fetch_archives_metadata_and_checks_fasta_to_completion(tmp_path, monkeypatch):
    content = b'>chr1\nACGTN\n>chrUn\nRY\n'
    sha = hashlib.sha256(content).hexdigest()
    monkeypatch.setattr(module, 'EXPECTED_SHA256', sha)
    def response(request, timeout):
        if isinstance(request, str):
            return Response((sha+' '+module.FILENAME).encode())
        return Response(content)
    monkeypatch.setattr(module, 'urlopen', response)
    result = module.fetch(tmp_path/'reference')
    assert result['complete'] and result['qc']['total_bases'] == 7
    assert result['qc']['n_bases'] == 1 and result['qc']['other_ambiguous_bases'] == 2
    assert result['qc']['input_sha256'] == sha
    monkeypatch.setattr(module, 'urlopen', lambda *a, **kw: pytest.fail('Unexpected network'))
    assert module.fetch(tmp_path/'reference')['complete']
