"""Full-file FASTQ integrity and distribution QC with disk-backed ID validation."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import sqlite3

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.fastq_stream import hashed_fastq, records


def length_quantile(histogram: Counter[int], fraction: float, weighted: bool=False) -> int:
    total=sum(length*count if weighted else count for length,count in histogram.items())
    cumulative=0
    for length,count in sorted(histogram.items()):
        cumulative += length*count if weighted else count
        if cumulative >= fraction*total:
            return length
    raise ValueError('Empty length distribution')


def n50(histogram: Counter[int]) -> int:
    target=sum(length*count for length,count in histogram.items())/2
    cumulative=0
    for length,count in sorted(histogram.items(),reverse=True):
        cumulative+=length*count
        if cumulative>=target:
            return length
    raise ValueError('Empty length distribution')


def qc(path: Path, outdir: Path, expected_reads: int | None=None, expected_bases: int | None=None) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    lengths, bases, qualities = Counter(), Counter(), Counter()
    joint = Counter()
    ids=sqlite3.connect(outdir/'read_ids.sqlite')
    ids.execute('PRAGMA cache_size=-16384')
    ids.execute('CREATE TABLE ids (name BLOB PRIMARY KEY) WITHOUT ROWID')
    count=total=0
    try:
        with hashed_fastq(path) as (handle, input_digest):
            for record in records(handle):
                seq, quality, identifier = record.sequence, record.quality, record.identifier
                try:
                    ids.execute('INSERT INTO ids VALUES (?)',(identifier,))
                except sqlite3.IntegrityError as exc:
                    raise ValueError(f'Duplicate FASTQ identifier: {identifier!r}') from exc
                count+=1;total+=len(seq)
                seq_counts=Counter(seq.upper()); q_counts=Counter(quality)
                lengths[len(seq)]+=1; bases.update(seq_counts);qualities.update(q_counts)
                gc=(seq_counts[ord('G')]+seq_counts[ord('C')])/len(seq)
                error=sum(n*10**(-(q-33)/10) for q,n in q_counts.items())/len(seq)
                mean_q=-10*math.log10(error)
                joint[(len(seq)//1000, min(100,int(100*gc)), int(mean_q//5)*5)]+=1
                if count%10000==0:
                    ids.commit()
        ids.commit()
        if not count:
            raise ValueError('Empty FASTQ input')
        if expected_reads is not None and count != expected_reads:
            raise ValueError(f'Read count {count} differs from expected {expected_reads}')
        if expected_bases is not None and total != expected_bases:
            raise ValueError(f'Base count {total} differs from expected {expected_bases}')
        mean_error=sum(n*10**(-(q-33)/10) for q,n in qualities.items())/total
        result=dict(complete=True,fastq_records_valid=True,gzip_trailer_checked=path.suffix=='.gz',
                    exact_duplicate_read_ids=0, read_count=count,total_bases=total,
                    min_length=min(lengths),max_length=max(lengths),median_length=length_quantile(lengths,.5),
                    # N50 uses descending cumulative bases, unlike a lower-tail quantile.
                    read_n50=n50(lengths),
                    gc_fraction=(bases[ord('G')]+bases[ord('C')])/total,n_fraction=bases[ord('N')]/total,
                    mean_reported_error_probability=mean_error,phred_from_mean_reported_error=-10*math.log10(mean_error),
                    quality_note='Phred+33 reported base-quality probabilities, not empirically measured read accuracy',
                    input_sha256=input_digest.hexdigest(),script_sha256=digest_file(Path(__file__)),
                    parser_sha256=digest_file(Path(__file__).with_name('fastq_stream.py')),
                    biological_qc='species/material/ploidy/background contamination/mapping bias not established by file QC')
        write_table(outdir/'length_histogram.tsv', [dict(length_bp=l,read_count=n) for l,n in sorted(lengths.items())],['length_bp','read_count'])
        write_table(outdir/'base_quality_histogram.tsv',[dict(phred=q-33,base_count=n) for q,n in sorted(qualities.items())],['phred','base_count'])
        write_table(outdir/'joint_distribution.tsv',[dict(length_bin_kb=l,gc_bin_percent=g,mean_quality_bin_phred=q,read_count=n)
                    for (l,g,q),n in sorted(joint.items())],['length_bin_kb','gc_bin_percent','mean_quality_bin_phred','read_count'])
        (outdir/'qc.json').write_text(json.dumps(result,indent=2)+'\n')
        return result
    except Exception as exc:
        (outdir/'qc.json').write_text(json.dumps(dict(complete=False,error=str(exc),records_before_failure=count),indent=2)+'\n')
        raise
    finally:
        ids.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--expected-reads',type=int)
    parser.add_argument('--expected-bases',type=int)
    args=parser.parse_args()
    qc(args.fastq,args.outdir,args.expected_reads,args.expected_bases)
