"""Acquire the versioned Mo17 T2T reference from its verified NCBI record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.fetch_ena_complete import download
from benchmarks.scripts.qc_reference_fasta import reference_qc


def fetch(outdir: Path) -> dict:
    accession = 'GCA_022117705.1'
    stem = accession+'_Zm-Mo17-REFERENCE-CAU-T2T-assembly'
    base = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/117/705/'+stem+'/'
    api = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+accession+'/dataset_report'
    outdir.mkdir(parents=True, exist_ok=True)
    plan_path = outdir/'reference_plan.json'
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        if plan['accession'] != accession or plan['ftp_directory'] != base:
            raise ValueError('Existing reference plan differs')
    else:
        metadata = {}
        for name, url in [('dataset_report.json', api), ('md5checksums.txt', base+'md5checksums.txt'),
                          (stem+'_assembly_report.txt', base+stem+'_assembly_report.txt'),
                          (stem+'_assembly_stats.txt', base+stem+'_assembly_stats.txt'),
                          (stem+'_fcs_report.txt', base+stem+'_fcs_report.txt')]:
            with urllib.request.urlopen(url, timeout=45) as response:
                content = response.read(2_000_001)
            if len(content) > 2_000_000:
                raise ValueError('Reference metadata exceeds size budget')
            (outdir/name).write_bytes(content)
            metadata[name] = dict(url=url, bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
        reports = json.loads((outdir/'dataset_report.json').read_text())['reports']
        if len(reports) != 1 or reports[0]['accession'] != accession or reports[0]['assembly_info']['bioproject_accession'] != 'PRJNA751841':
            raise ValueError('NCBI reference accession/project does not match source study')
        checksums = {}
        for line in (outdir/'md5checksums.txt').read_text().splitlines():
            digest, filename = line.split(maxsplit=1)
            checksums[filename.removeprefix('./')] = digest
        for name in metadata:
            if name in checksums and hashlib.md5((outdir/name).read_bytes()).hexdigest() != checksums[name]:
                raise ValueError('Metadata does not match official MD5: '+name)
        filename = stem+'_genomic.fna.gz'
        with urllib.request.urlopen(urllib.request.Request(base+filename, method='HEAD'), timeout=45) as response:
            size = int(response.headers['Content-Length'])
        if not 0 < size <= 1_000_000_000:
            raise ValueError('Reference exceeds one-GB compressed download budget')
        plan = dict(accession=accession, ftp_directory=base, metadata=metadata,
                    source_paper='https://doi.org/10.1038/s41588-023-01419-6',
                    expected_total_bases=int(reports[0]['assembly_stats']['total_sequence_length']),
                    expected_contigs=int(reports[0]['assembly_stats']['number_of_contigs']),
                    source_script_sha256=digest_file(Path(__file__)),
                    file=dict(filename=filename, url=base+filename, expected_bytes=size, expected_md5=checksums[filename]))
        plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    transfer = download(plan['file'], outdir)
    result = reference_qc(outdir/plan['file']['filename'])
    if result['total_bases'] != plan['expected_total_bases'] or result['contig_count'] != plan['expected_contigs']:
        raise ValueError('Reference FASTA totals differ from archived NCBI metadata')
    receipt = dict(complete=True, accession=accession, transfer=transfer, qc=result,
                   plan_sha256=digest_file(plan_path), script_sha256=digest_file(Path(__file__)),
                   qc_script_sha256=digest_file(Path(__file__).with_name('qc_reference_fasta.py')),
                   parser_sha256=digest_file(Path(__file__).with_name('fastq_stream.py')),
                   downloader_sha256=digest_file(Path(__file__).with_name('fetch_ena_complete.py')))
    (outdir/'reference_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    fetch(parser.parse_args().outdir)
