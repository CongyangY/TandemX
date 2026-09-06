"""Acquire MorexV3 from its original, SHA-256-pinned IPK data publication."""
from __future__ import annotations

import argparse
import hashlib
import http.client
import json
from pathlib import Path
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.qc_reference_fasta import reference_qc


LANDING = ('https://doi.ipk-gatersleben.de/DOI/b2f47dfb-47ff-4114-89ae-bad8dcc515a1/'
           'b6e6a2e5-2746-4522-8465-019c8f56df7f/1')
FILENAME = 'Barley_MorexV3_pseudomolecules.fasta'
EXPECTED_SHA256 = '54c98a04d13ff97350f5f3a5bfa45ac395ad640df8bb1f7598eca4e7edb437c1'
MAX_BYTES = 4_500_000_000


class DivergentPartialError(ValueError):
    """The server ignored Range and no longer matches the retained partial."""

    def __init__(self, first_difference_offset: int, partial_sha256: str):
        super().__init__(
            f'Full-response prefix differs from partial at byte {first_difference_offset}; '
            'partial retained'
        )
        self.first_difference_offset = first_difference_offset
        self.partial_sha256 = partial_sha256


def download_sha256(url: str, path: Path, expected_sha256: str, max_bytes: int) -> dict:
    """Bound a chunked transfer and verify source SHA before renaming.

    A partial can be resumed only when the server acknowledges its exact range.
    HTTP access-denied responses are not retried. No partial is deleted on error.
    """
    if max_bytes <= 0 or len(expected_sha256) != 64 or set(expected_sha256)-set('0123456789abcdef'):
        raise ValueError('Require a positive byte budget and lowercase SHA-256')
    partial = path.with_name(path.name+'.partial')
    if path.exists():
        if not 0 < path.stat().st_size <= max_bytes or digest_file(path) != expected_sha256:
            raise ValueError('Existing reference differs from the pinned source')
        return dict(bytes=path.stat().st_size, sha256=expected_sha256, transfer='existing_reverified')
    transfer = 'downloaded_complete'
    for attempt in range(1, 6):
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > max_bytes:
            raise ValueError('Partial reference exceeds byte budget')
        sha = hashlib.sha256()
        if offset:
            with partial.open('rb') as handle:
                for block in iter(lambda: handle.read(8*1024*1024), b''):
                    sha.update(block)
        if offset and sha.hexdigest() == expected_sha256:
            partial.rename(path)
            return dict(bytes=offset, sha256=expected_sha256, transfer='complete_partial_reverified')
        request = Request(url, headers={'Range': f'bytes={offset}-'} if offset else {})
        try:
            with urlopen(request, timeout=60) as response:
                if response.status not in {200, 206}:
                    raise ValueError(f'Unexpected reference response: {response.status}')
                if offset and response.status == 200:
                    # Some publication repositories ignore Range. Retain local
                    # progress only after the complete response prefix is proven
                    # byte-identical to the partial already on disk.
                    compared = 0
                    with partial.open('rb') as existing:
                        while expected := existing.read(1024*1024):
                            observed = response.read(len(expected))
                            if observed != expected:
                                difference = next(
                                    (index for index, pair in enumerate(zip(expected, observed))
                                     if pair[0] != pair[1]),
                                    min(len(expected), len(observed)),
                                )
                                raise DivergentPartialError(compared + difference, sha.hexdigest())
                            compared += len(expected)
                    transfer = 'downloaded_complete_after_verified_prefix'
                if response.status == 206 and not response.headers.get('Content-Range', '').startswith(f'bytes {offset}-'):
                    raise ValueError('Incorrect reference resume range')
                with partial.open('ab' if offset else 'wb') as handle:
                    while block := response.read(1024*1024):
                        if offset+len(block) > max_bytes:
                            raise ValueError('Reference transfer exceeds byte budget')
                        handle.write(block)
                        sha.update(block)
                        offset += len(block)
            if not offset or sha.hexdigest() != expected_sha256:
                raise ValueError('Reference SHA-256 differs from the published source; partial retained')
            partial.rename(path)
            return dict(bytes=offset, sha256=expected_sha256, transfer=transfer)
        except (OSError, http.client.HTTPException) as exc:
            with path.with_name(path.name+'.transfer.log').open('a') as log:
                log.write(f'attempt={attempt}\tbytes={offset}\terror={exc}\n')
            if attempt == 5 or isinstance(exc, HTTPError) and exc.code not in {408, 429, 500, 502, 503, 504}:
                raise
            time.sleep(2)
    raise RuntimeError('Reference acquisition exhausted retries')


def fetch(outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    plan_path = outdir/'reference_plan.json'
    identity = dict(filename=FILENAME, landing_url=LANDING, url=LANDING+'/DOWNLOAD',
                    expected_sha256=EXPECTED_SHA256, max_bytes=MAX_BYTES,
                    doi='10.5447/ipk/2021/3', license='CC BY 4.0')
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        if any(plan.get(k) != v for k, v in identity.items()):
            raise ValueError('Existing Morex reference plan differs')
        if digest_file(outdir/'source_metadata.html') != plan['metadata_sha256']:
            raise ValueError('Archived Morex source metadata changed')
    else:
        with urlopen(LANDING, timeout=45) as response:
            content = response.read(2_000_001)
        if len(content) > 2_000_000 or EXPECTED_SHA256.encode() not in content or FILENAME.encode() not in content:
            raise ValueError('Morex metadata does not establish the pinned filename/checksum')
        (outdir/'source_metadata.html').write_bytes(content)
        plan = dict(**identity, metadata_sha256=hashlib.sha256(content).hexdigest(),
                    script_sha256=digest_file(Path(__file__)))
        plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    receipt = dict(complete=False, scope='original_published_MorexV3_pseudomolecules',
                   plan_sha256=digest_file(plan_path), script_sha256=digest_file(Path(__file__)),
                   qc_script_sha256=digest_file(Path(__file__).with_name('qc_reference_fasta.py')),
                   parser_sha256=digest_file(Path(__file__).with_name('fastq_stream.py')),
                   warning='same_cultivar_study_context_not_verified_identical_donor_or_satellite_copy_truth')
    try:
        reference = outdir/FILENAME
        try:
            receipt['transfer'] = download_sha256(plan['url'], reference, EXPECTED_SHA256, MAX_BYTES)
        except DivergentPartialError as exc:
            partial = reference.with_name(reference.name+'.partial')
            receipt['preserved_divergent_partial'] = {
                'path': str(partial.resolve()),
                'bytes': partial.stat().st_size,
                'sha256': exc.partial_sha256,
                'first_difference_offset': exc.first_difference_offset,
            }
            clean = reference.with_name(reference.name+'.clean')
            clean_transfer = download_sha256(plan['url'], clean, EXPECTED_SHA256, MAX_BYTES)
            clean.rename(reference)
            receipt['transfer'] = {
                **clean_transfer,
                'recovery': 'clean_download_after_divergent_partial',
                'divergent_partial_retained': str(partial.resolve()),
            }
        receipt['qc'] = reference_qc(reference)
        if receipt['qc']['input_sha256'] != EXPECTED_SHA256:
            raise ValueError('Reference changed between download verification and FASTA QC')
        receipt['complete'] = True
        return receipt
    except Exception as exc:
        receipt['error'] = str(exc)
        raise
    finally:
        (outdir/'reference_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    fetch(parser.parse_args().outdir)
