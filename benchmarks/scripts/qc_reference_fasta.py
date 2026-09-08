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
        at_line_start = True
        while fragment := handle.readline(1024 * 1024):
            line_terminated = fragment.endswith(b'\n')
            if at_line_start and fragment.startswith(b'>'):
                if not line_terminated:
                    raise ValueError('Reference header exceeds bounded parser limit')
                line = fragment.rstrip(b'\n\r')
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
                at_line_start = True
                continue
            line = fragment.rstrip(b'\n\r') if line_terminated else fragment
            if line:
                if current is None or set(line.upper()) - set(b'ACGTRYSWKMBDHVN'):
                    raise ValueError('Invalid reference sequence or missing header')
                count = Counter(line.upper().decode('ascii'))
                counts.update(count)
                total.update(count)
            at_line_start = line_terminated
        if current is None or not counts:
            raise ValueError('Empty reference input or terminal contig')
        contigs.append(dict(contig=current, length_bp=sum(counts.values()), base_counts=dict(counts)))
    return dict(complete=True, contig_count=len(contigs), total_bases=sum(total.values()),
                n_bases=total['N'], other_ambiguous_bases=sum(n for b, n in total.items() if b not in 'ACGTN'),
                base_counts=dict(total), input_sha256=digest.hexdigest(), contigs=contigs,
                biological_qc='file content only; assembly completeness and satellite copy truth not established')
