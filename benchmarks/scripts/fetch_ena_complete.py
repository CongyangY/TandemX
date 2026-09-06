"""Download complete ENA FASTQ files with resumable transfer and exact checksums."""
from __future__ import annotations

import argparse
import hashlib
import http.client
import json
from pathlib import Path
import re
import time
import urllib.request
from urllib.parse import urlsplit

from benchmarks.challenge.schema import read_table, digest_file


def files_for_run(metadata: Path, accession: str) -> tuple[dict, list[dict]]:
    rows = [r for r in read_table(metadata, {"run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes"}) if r["run_accession"] == accession]
    if len(rows) != 1 or not re.fullmatch(r"[ESD]RR\d+", accession):
        raise ValueError("Metadata must identify exactly one requested ENA run")
    row = rows[0]
    groups = [row[key].split(";") for key in ("fastq_ftp", "fastq_md5", "fastq_bytes")]
    if len({len(group) for group in groups}) != 1:
        raise ValueError("FASTQ URL/hash/size cardinality differs")
    files = []
    for path, md5, size in zip(*groups):
        url = path if path.startswith("https://") else "https://" + path
        parsed = urlsplit(url)
        name = Path(parsed.path).name
        if (parsed.scheme != "https" or parsed.hostname != "ftp.sra.ebi.ac.uk" or parsed.query
            or not re.fullmatch(r"[A-Za-z0-9_.-]+\.f(?:ast)?q\.gz", name)
            or not re.fullmatch(r"[0-9a-f]{32}", md5) or not size.isdigit() or int(size) <= 0):
            raise ValueError("Invalid complete ENA FASTQ metadata")
        files.append(dict(url=url, filename=name, expected_md5=md5, expected_bytes=int(size)))
    if len({r['filename'] for r in files}) != len(files):
        raise ValueError("Duplicate target filename")
    return row, files


def verify_file(path: Path, expected_bytes: int, expected_md5: str) -> dict:
    if path.stat().st_size != expected_bytes:
        raise ValueError(f"Size mismatch: {path}")
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            md5.update(chunk); sha.update(chunk)
    if md5.hexdigest() != expected_md5:
        raise ValueError(f"MD5 mismatch: {path}; retain file for investigation")
    return dict(bytes=expected_bytes, md5=md5.hexdigest(), sha256=sha.hexdigest(), source_md5_verified=True)


def download(entry: dict, outdir: Path) -> dict:
    final = outdir / entry['filename']
    partial = final.with_name(final.name+'.partial')
    if final.exists():
        return dict(entry, **verify_file(final, entry['expected_bytes'], entry['expected_md5']), transfer='existing_reverified')
    for attempt in range(1, 6):
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > entry['expected_bytes']:
            raise ValueError(f"Partial file exceeds expected size: {partial}")
        if offset == entry['expected_bytes']:
            break
        request = urllib.request.Request(entry['url'], headers={'Range':f'bytes={offset}-'} if offset else {})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if offset and response.status == 200:
                    raise ValueError('Server ignored resume range; retain partial instead of appending duplicate data')
                if response.status == 206 and not response.headers.get('Content-Range','').startswith(f'bytes {offset}-'):
                    raise ValueError('Resume response has incorrect byte range')
                with partial.open('ab' if offset else 'wb') as handle:
                    while chunk := response.read(1024*1024):
                        if offset + len(chunk) > entry['expected_bytes']:
                            raise ValueError('Transfer exceeds the metadata byte budget')
                        handle.write(chunk)
                        offset += len(chunk)
            if offset != entry['expected_bytes']:
                raise OSError('Incomplete transfer; retrying at last byte')
            break
        except (OSError, TimeoutError, http.client.HTTPException) as exc:
            with (outdir/'transfer.log').open('a') as log:
                log.write(f"attempt={attempt}\tfile={final.name}\tbytes={offset}\terror={exc}\n")
            if attempt == 5:
                raise
            time.sleep(2)
    verified = verify_file(partial, entry['expected_bytes'], entry['expected_md5'])
    partial.rename(final)
    return dict(entry, **verified, transfer='downloaded_complete')


def fetch(metadata: Path, accession: str, outdir: Path, max_gb: float) -> None:
    row, files = files_for_run(metadata, accession)
    if not 0 < max_gb < float('inf') or sum(r['expected_bytes'] for r in files) > max_gb*1e9:
        raise ValueError('Requested complete files exceed the explicit download budget')
    outdir.mkdir(parents=True, exist_ok=True)
    plan = dict(accession=accession, metadata=row, metadata_sha256=digest_file(metadata), files=files,
                max_download_gb=max_gb, sampling='complete_ENA_run', script_sha256=digest_file(Path(__file__)))
    plan_path = outdir/'download_plan.json'
    if plan_path.exists():
        previous=json.loads(plan_path.read_text())
        if previous['files'] != files or previous['accession'] != accession:
            raise ValueError('Existing download plan differs; choose a new directory')
    else:
        plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    receipts = []
    for entry in files:
        receipts.append(download(entry,outdir))
        (outdir/'download_receipt.json').write_text(json.dumps(dict(plan, completed_files=receipts,
            complete=len(receipts)==len(files), fastq_qc_status='not_run'),indent=2)+'\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--accession', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--max-download-gb', type=float, required=True)
    args=parser.parse_args()
    fetch(args.metadata,args.accession,args.outdir,args.max_download_gb)
