"""Acquire a pinned GenBank assembly with official checksums and streaming QC."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import urllib.request
from urllib.parse import urljoin, urlsplit

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.fetch_ena_complete import download
from benchmarks.scripts.qc_reference_fasta import reference_qc


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.links.extend(v for k, v in attrs if k == 'href' and v)


def assembly_parent(accession: str) -> str:
    match = re.fullmatch(r'GCA_(\d{3})(\d{3})(\d{3})\.[1-9]\d*', accession)
    if not match:
        raise ValueError('Require an explicitly versioned GCA accession')
    return 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/'+'/'.join(match.groups())+'/'


def resolve_directory(parent: str, accession: str, listing: str) -> str:
    links = Links()
    links.feed(listing)
    matches = {urljoin(parent, link) for link in links.links
               if link.startswith(accession+'_') and link.endswith('/')
               and '/' not in link[:-1] and '?' not in link and '#' not in link}
    if len(matches) != 1:
        raise ValueError('No unique versioned assembly directory in official NCBI listing')
    return matches.pop()


def fetch(accession: str, bioproject: str, outdir: Path, source_paper: str, max_bytes: int) -> dict:
    parent = assembly_parent(accession)
    if not re.fullmatch(r'PRJ[A-Z]{2}\d+', bioproject) or max_bytes <= 0:
        raise ValueError('Require a BioProject and positive compressed-size budget')
    outdir.mkdir(parents=True, exist_ok=True)
    plan_path = outdir/'reference_plan.json'
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        if any(plan[k] != v for k, v in dict(accession=accession, bioproject=bioproject,
                                            source_paper=source_paper, max_bytes=max_bytes).items()):
            raise ValueError('Existing reference plan differs from request')
        if any(digest_file(outdir/n) != m['sha256'] for n, m in plan['metadata'].items()):
            raise ValueError('Archived reference metadata changed')
    else:
        metadata = {}

        def retrieve(name, url):
            with urllib.request.urlopen(url, timeout=45) as response:
                content = response.read(2_000_001)
            if len(content) > 2_000_000:
                raise ValueError('Reference metadata exceeds 2 MB cap')
            (outdir/name).write_bytes(content)
            metadata[name] = dict(url=url, bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
            return content

        api = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+accession+'/dataset_report'
        report = json.loads(retrieve('dataset_report.json', api))['reports']
        if len(report) != 1 or report[0]['accession'] != accession or report[0]['assembly_info']['bioproject_accession'] != bioproject:
            raise ValueError('Assembly accession or source project disagrees with NCBI')
        report = report[0]
        listing = retrieve('assembly_directory.html', parent).decode()
        base = resolve_directory(parent, accession, listing)
        stem = urlsplit(base).path.rstrip('/').rsplit('/', 1)[1]
        checksums = {}
        for line in retrieve('md5checksums.txt', base+'md5checksums.txt').decode().splitlines():
            value, name = line.split(maxsplit=1)
            name = name.removeprefix('./')
            if name in checksums or not re.fullmatch(r'[a-fA-F0-9]{32}', value):
                raise ValueError('Malformed official checksum manifest')
            checksums[name] = value.lower()
        optional = {}
        for suffix in ('_assembly_report.txt', '_assembly_stats.txt', '_fcs_report.txt'):
            name = stem+suffix
            if name not in checksums and suffix == '_fcs_report.txt':
                optional[name] = 'not_listed_in_official_checksum_manifest'
                continue
            content = retrieve(name, base+name)
            if hashlib.md5(content).hexdigest() != checksums[name]:
                raise ValueError('Reference metadata MD5 mismatch: '+name)
        filename = stem+'_genomic.fna.gz'
        with urllib.request.urlopen(urllib.request.Request(base+filename, method='HEAD'), timeout=45) as response:
            size = int(response.headers['Content-Length'])
        if not 0 < size <= max_bytes:
            raise ValueError('Reference exceeds explicit compressed-size budget')
        plan = dict(accession=accession, bioproject=bioproject, source_paper=source_paper,
                    max_bytes=max_bytes, ftp_directory=base, metadata=metadata, unavailable_optional_metadata=optional,
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    expected_total_bases=int(report['assembly_stats']['total_sequence_length']),
                    expected_contigs=int(report['assembly_stats']['number_of_contigs']),
                    file=dict(filename=filename, url=base+filename, expected_bytes=size, expected_md5=checksums[filename]))
        plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    sources = {str(p): digest_file(p) for p in (Path(__file__), Path(__file__).with_name('fetch_ena_complete.py'),
               Path(__file__).with_name('qc_reference_fasta.py'), Path(__file__).with_name('fastq_stream.py'))}
    receipt = dict(complete=False, accession=accession, plan_sha256=digest_file(plan_path), source_hashes=sources)
    try:
        receipt['transfer'] = download(plan['file'], outdir)
        receipt['qc'] = reference_qc(outdir/plan['file']['filename'])
        if (receipt['qc']['total_bases'] != plan['expected_total_bases']
                or receipt['qc']['contig_count'] != plan['expected_contigs']):
            raise ValueError('Complete FASTA totals disagree with official NCBI metadata')
        if any(digest_file(Path(p)) != h for p, h in sources.items()):
            raise ValueError('Acquisition/QC source changed during execution')
        receipt['complete'] = True
        return receipt
    except Exception as exc:
        receipt['error'] = str(exc)
        raise
    finally:
        (outdir/'reference_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accession', required=True)
    parser.add_argument('--bioproject', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--source-paper', required=True)
    parser.add_argument('--max-download-gb', type=float, required=True)
    args = parser.parse_args()
    fetch(args.accession, args.bioproject, args.outdir, args.source_paper, int(args.max_download_gb*1e9))
