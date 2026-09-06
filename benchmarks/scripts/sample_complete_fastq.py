"""Nested deterministic Bernoulli sampling across a QC-verified complete FASTQ."""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import logging
import math
from pathlib import Path
import re

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.fastq_stream import hashed_fastq, records
from benchmarks.scripts.qc_complete_fastq import length_quantile, n50


def fraction_threshold(value: str) -> tuple[str, int]:
    try:
        fraction = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError('Invalid sampling fraction') from exc
    if not fraction.is_finite() or not 0 < fraction <= 1:
        raise ValueError('Sampling fractions must be finite and in (0, 1]')
    numerator, denominator = fraction.as_integer_ratio()
    threshold = numerator * (1 << 128) // denominator
    if threshold == 0:
        raise ValueError('Sampling fraction is below 128-bit resolution')
    return str(fraction), threshold


def read_hash(identifier: bytes, namespace: str, seed: int) -> int:
    source = namespace.encode('ascii')
    payload = len(source).to_bytes(4, 'big') + source + identifier
    return int.from_bytes(hashlib.blake2b(
        payload, digest_size=16, key=seed.to_bytes(8, 'big'),
        person=b'TandemX-smp-v1').digest(), 'big')


def sample(path: Path, qc_receipt: Path, outdir: Path, fractions: list[str],
           namespace: str, seed: int, genome_size: int | None = None) -> dict:
    """Membership is order independent; emitted record order follows the source."""
    if not re.fullmatch(r'[A-Za-z0-9_.-]{1,128}', namespace):
        raise ValueError('Namespace must be 1-128 ASCII accession/filename-safe characters')
    if not 0 <= seed < 1 << 64:
        raise ValueError('Seed must be an unsigned 64-bit integer')
    if genome_size is not None and genome_size <= 0:
        raise ValueError('Genome size must be positive')
    thresholds = sorted((fraction_threshold(f) for f in fractions), key=lambda x: x[1])
    if not 1 <= len(thresholds) <= 12 or len({t for _, t in thresholds}) != len(thresholds):
        raise ValueError('Provide 1-12 distinct sampling fractions')
    verified = json.loads(qc_receipt.read_text())
    if (verified.get('complete') is not True or verified.get('fastq_records_valid') is not True
            or verified.get('exact_duplicate_read_ids') != 0
            or not re.fullmatch('[0-9a-f]{64}', verified.get('input_sha256', ''))
            or verified.get('read_count', 0) <= 0 or verified.get('total_bases', 0) <= 0):
        raise ValueError('A successful full-file QC receipt with unique IDs is required')
    outdir.mkdir(parents=True, exist_ok=False)
    states = [dict(sample_id=f'sample_{i:03d}', fraction=f, hash_threshold=str(t),
                   threshold=t, reads=0, bases=0, lengths=Counter(), base_counts=Counter(),
                   quality_counts=Counter(), joint=Counter())
              for i, (f, t) in enumerate(thresholds, 1)]
    plan = dict(method='blake2b_128_read_bernoulli_v1', seed=seed, namespace=namespace,
                fractions=[f for f, _ in thresholds], source_path=str(path.resolve()),
                expected_input_sha256=verified['input_sha256'], qc_receipt_sha256=digest_file(qc_receipt),
                genome_size_denominator=genome_size, script_sha256=digest_file(Path(__file__)),
                parser_sha256=digest_file(Path(__file__).with_name('fastq_stream.py')),
                uncertainty='Independent read inclusion, not fixed bases or independent specimens; library bias persists',
                coverage_note='Total selected bases / supplied size is nominal coverage, not measured nuclear depth')
    (outdir/'sampling_plan.json').write_text(json.dumps(plan, indent=2)+'\n')
    read_count = total_bases = 0
    try:
        with ExitStack() as stack:
            for state in states:
                raw = stack.enter_context((outdir/(state['sample_id']+'.fastq.gz.partial')).open('wb'))
                state['writer'] = stack.enter_context(gzip.GzipFile(filename='', fileobj=raw, mode='wb',
                                                                   compresslevel=1, mtime=0))
                state['ids'] = stack.enter_context((outdir/(state['sample_id']+'.ids.tsv.partial')).open('wb'))
                state['ids'].write(b'read_id_hex\thash128_hex\tlength_bp\n')
            with hashed_fastq(path) as (handle, input_digest):
                for record in records(handle):
                    read_count += 1
                    total_bases += len(record.sequence)
                    rank = read_hash(record.identifier, namespace, seed)
                    if rank < states[-1]['threshold']:
                        seq_counts = Counter(record.sequence.upper())
                        q_counts = Counter(record.quality)
                        gc = (seq_counts[ord('G')]+seq_counts[ord('C')])/len(record.sequence)
                        error = sum(n*10**(-(q-33)/10) for q, n in q_counts.items())/len(record.sequence)
                        joint = (len(record.sequence)//1000, min(100, int(gc*100)),
                                 int((-10*math.log10(error))//5)*5)
                        encoded = record.as_bytes()
                        id_row = f'{record.identifier.hex()}\t{rank:032x}\t{len(record.sequence)}\n'.encode('ascii')
                        for state in states:
                            if rank < state['threshold']:
                                state['writer'].write(encoded)
                                state['ids'].write(id_row)
                                state['reads'] += 1
                                state['bases'] += len(record.sequence)
                                state['lengths'][len(record.sequence)] += 1
                                state['base_counts'].update(seq_counts)
                                state['quality_counts'].update(q_counts)
                                state['joint'][joint] += 1
                    if read_count % 100_000 == 0:
                        logging.info('Scanned %s reads / %s bases', read_count, total_bases)
            if input_digest.hexdigest() != verified['input_sha256']:
                raise ValueError('Full input SHA256 differs from the supplied QC receipt')
            if read_count != verified['read_count'] or total_bases != verified['total_bases']:
                raise ValueError('Observed input totals differ from full QC')
        summaries = []
        for state in states:
            prefix = outdir/state['sample_id']
            for suffix in ('.fastq.gz', '.ids.tsv'):
                Path(str(prefix)+suffix+'.partial').rename(str(prefix)+suffix)
            count, bases = state['reads'], state['bases']
            summary = dict(sample_id=state['sample_id'], inclusion_probability=state['fraction'],
                           hash_threshold=state['hash_threshold'], read_count=count, total_bases=bases,
                           expected_read_count=float(Decimal(state['fraction'])*read_count),
                           expected_total_bases=float(Decimal(state['fraction'])*total_bases),
                           nominal_total_base_coverage=bases/genome_size if genome_size else None,
                           read_n50=n50(state['lengths']) if count else None,
                           median_length=length_quantile(state['lengths'], .5) if count else None,
                           gc_fraction=(state['base_counts'][71]+state['base_counts'][67])/bases if bases else None,
                           n_fraction=state['base_counts'][78]/bases if bases else None,
                           mean_reported_error_probability=sum(n*10**(-(q-33)/10) for q, n in state['quality_counts'].items())/bases if bases else None,
                           status='ok' if count else 'empty_random_sample',
                           fastq_sha256=digest_file(Path(str(prefix)+'.fastq.gz')),
                           ids_sha256=digest_file(Path(str(prefix)+'.ids.tsv')))
            summaries.append(summary)
            write_table(Path(str(prefix)+'.length_histogram.tsv'),
                        [dict(length_bp=l, read_count=n) for l, n in sorted(state['lengths'].items())],
                        ['length_bp', 'read_count'])
            write_table(Path(str(prefix)+'.joint_distribution.tsv'),
                        [dict(length_bin_kb=l, gc_bin_percent=g, mean_quality_bin_phred=q, read_count=n)
                         for (l, g, q), n in sorted(state['joint'].items())],
                        ['length_bin_kb', 'gc_bin_percent', 'mean_quality_bin_phred', 'read_count'])
        result = dict(complete=True, source_sha256=input_digest.hexdigest(),
                      source_read_count=read_count, source_total_bases=total_bases,
                      plan_sha256=digest_file(outdir/'sampling_plan.json'), samples=summaries)
        (outdir/'sampling_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
        return result
    except Exception as exc:
        (outdir/'sampling_receipt.json').write_text(json.dumps(
            dict(complete=False, error=str(exc), scanned_reads=read_count, scanned_bases=total_bases), indent=2)+'\n')
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq', type=Path, required=True)
    parser.add_argument('--qc-receipt', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--fractions', nargs='+', required=True, help='Read inclusion probabilities in (0,1]')
    parser.add_argument('--namespace', required=True, help='Stable included-library accession')
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--genome-size', type=int)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    sample(args.fastq, args.qc_receipt, args.outdir, args.fractions, args.namespace, args.seed, args.genome_size)
