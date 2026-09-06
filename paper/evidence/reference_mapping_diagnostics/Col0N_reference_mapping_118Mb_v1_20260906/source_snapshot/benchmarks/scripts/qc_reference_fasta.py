"""Streaming reference FASTA QC without retaining complete chromosomes."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from benchmarks.scripts.fastq_stream import hashed_fastq


def reference_qc(path: Path) -> dict:
    contigs, identifiers = [], set()
    current = None
    total = Counter()
    counts = Counter()
    with hashed_fastq(path) as (handle, digest):
        while line := handle.readline(5_000_002):
            if len(line) >= 5_000_002:
                raise ValueError('Reference line exceeds bounded parser limit')
            line = line.strip()
            if not line:
                continue
            if line.startswith(b'>'):
                if current is not None:
                    if not counts:
                        raise ValueError('Empty reference contig')
                    contigs.append(dict(contig=current, length_bp=sum(counts.values()), base_counts=dict(counts)))
                fields = line[1:].split()
                if not fields:
                    raise ValueError('Empty FASTA identifier')
                current = fields[0].decode('ascii')
                if current in identifiers:
                    raise ValueError('Duplicate reference identifier')
                identifiers.add(current)
                counts = Counter()
            else:
                if current is None or set(line.upper())-set(b'ACGTRYSWKMBDHVN'):
                    raise ValueError('Invalid reference sequence or missing header')
                count = Counter(line.upper().decode('ascii'))
                counts.update(count)
                total.update(count)
        if current is None or not counts:
            raise ValueError('Empty reference input or terminal contig')
        contigs.append(dict(contig=current, length_bp=sum(counts.values()), base_counts=dict(counts)))
    return dict(complete=True, contig_count=len(contigs), total_bases=sum(total.values()),
                n_bases=total['N'], other_ambiguous_bases=sum(n for b, n in total.items() if b not in 'ACGTN'),
                base_counts=dict(total), input_sha256=digest.hexdigest(), contigs=contigs,
                biological_qc='file content only; assembly completeness and satellite copy truth not established')
