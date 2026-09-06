"""Disk-backed, alignment-aware reference QC; mapping is not biological truth."""
from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path
import re
import sqlite3

from benchmarks.challenge.schema import write_table
from benchmarks.scripts.fastq_stream import hashed_fastq, records
from benchmarks.scripts.real_disk import connect_index


def prepare_queries(fastq: Path, database: Path, expected: dict) -> dict:
    """Index IDs, lengths and composition without retaining sequences."""
    connection = connect_index(database, create=True)
    count = bases = 0
    try:
        connection.execute('CREATE TABLE reads (name TEXT PRIMARY KEY, length INTEGER, gc_bin INTEGER, purine_bin INTEGER) WITHOUT ROWID')
        with hashed_fastq(fastq) as (handle, digest):
            for record in records(handle):
                sequence = record.sequence.upper()
                length = len(sequence)
                gc = sequence.count(b'G')+sequence.count(b'C')
                purine = sequence.count(b'G')+sequence.count(b'A')
                pyrimidine = sequence.count(b'C')+sequence.count(b'T')
                try:
                    connection.execute('INSERT INTO reads VALUES (?,?,?,?)',
                        (record.identifier.decode('ascii'), length, min(19, 20*gc//length),
                         min(19, 20*max(purine, pyrimidine)//length)))
                except sqlite3.IntegrityError as exc:
                    raise ValueError('Duplicate query read ID') from exc
                count += 1
                bases += length
                if count % 10000 == 0:
                    connection.commit()
        if count != expected['read_count'] or bases != expected['total_bases'] or digest.hexdigest() != expected['fastq_sha256']:
            raise ValueError('Mapping queries differ from completed sample receipt')
        connection.commit()
        return dict(read_count=count, total_bases=bases, fastq_sha256=digest.hexdigest())
    finally:
        connection.close()


def cigar_blocks(cigar: str, query_span: int, target_start: int, target_end: int) -> list[tuple[int, int]]:
    """Validate PAF CIGAR and retain target M/= /X blocks, excluding deletions."""
    target, query, previous = target_start, 0, 0
    blocks = []
    for match in re.finditer(r'([1-9][0-9]*)([MIDN=X])', cigar):
        if match.start() != previous:
            raise ValueError('Invalid or unsupported PAF CIGAR')
        previous = match.end()
        length, op = int(match[1]), match[2]
        if op in 'M=X':
            blocks.append((target, target+length))
        if op in 'MDN=X':
            target += length
        if op in 'MI=X':
            query += length
    if previous != len(cigar) or not cigar or target != target_end or query != query_span:
        raise ValueError('PAF CIGAR consumption disagrees with coordinates')
    return blocks


def load_alignments(paf: Path, database: Path, reference_lengths: dict[str, int]) -> dict:
    """Strict full PAF validation, retaining native primary and secondary rows."""
    connection = connect_index(database)
    rows = 0
    kinds = defaultdict(int)
    try:
        connection.execute('CREATE TABLE alignments (name TEXT, start INTEGER, end INTEGER, target TEXT, kind TEXT, mapq INTEGER, matches INTEGER, block INTEGER)')
        connection.execute('CREATE TABLE target_blocks (target TEXT, start INTEGER, end INTEGER, mapq INTEGER)')
        with paf.open() as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 12:
                    raise ValueError('Malformed PAF line')
                name, qlen, start, end, strand, target, tlen, ts, te, matches, block, mapq = fields[:12]
                qlen, start, end, tlen, ts, te, matches, block, mapq = map(int, (qlen, start, end, tlen, ts, te, matches, block, mapq))
                tags = {}
                for tag in fields[12:]:
                    parts = tag.split(':', 2)
                    if len(parts) != 3 or parts[0] in tags:
                        raise ValueError('Malformed or duplicate PAF tag')
                    tags[parts[0]] = (parts[1], parts[2])
                length = connection.execute('SELECT length FROM reads WHERE name=?', (name,)).fetchone()
                if (length != (qlen,) or reference_lengths.get(target) != tlen or strand not in {'+', '-'}
                        or not 0 <= start < end <= qlen or not 0 <= ts < te <= tlen
                        or not 0 <= matches <= block or block < 1 or not 0 <= mapq <= 255):
                    raise ValueError('Unknown query/target or invalid PAF coordinates/counts')
                if tags.get('tp', ('', '')) not in {('A', x) for x in ('P', 'S', 'I', 'i')} or tags.get('cg', ('', ''))[0] != 'Z':
                    raise ValueError('Require typed minimap2 alignment and CIGAR tags')
                kind = tags['tp'][1]
                blocks = cigar_blocks(tags['cg'][1], end-start, ts, te)
                if matches > sum(b-a for a, b in blocks):
                    raise ValueError('PAF match count exceeds CIGAR-aligned bases')
                connection.execute('INSERT INTO alignments VALUES (?,?,?,?,?,?,?,?)', (name, start, end, target, kind, mapq, matches, block))
                if kind in {'P', 'I'}:
                    connection.executemany('INSERT INTO target_blocks VALUES (?,?,?,?)', [(target, a, b, mapq) for a, b in blocks])
                rows += 1
                kinds[kind] += 1
                if rows % 10000 == 0:
                    connection.commit()
        connection.execute('CREATE INDEX alignment_query ON alignments(name,start,end)')
        connection.execute('CREATE INDEX block_target ON target_blocks(target,start,end)')
        connection.commit()
        return dict(alignment_rows=rows, alignment_types=dict(kinds))
    finally:
        connection.close()


def summarize_mapping(database: Path, reference_lengths: dict[str, int], organelles: set[str], outdir: Path) -> dict:
    """Query-span unions and CIGAR-aligned target bases with all-read denominators."""
    if organelles-set(reference_lengths):
        raise ValueError('Unknown organellar reference contig')
    connection = connect_index(database)
    groups = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    try:
        fields = ['read_id', 'length_bp', 'gc_bin_5pct', 'max_purine_pyrimidine_bin_5pct',
                  'any_alignment_span_bp', 'primary_span_bp', 'primary_mapq20_span_bp',
                  'primary_organelle_span_bp', 'primary_other_reference_span_bp',
                  'primary_rows', 'secondary_rows', 'organelle_primary_rows', 'primary_compartment_overlap_bp']
        query = connection.execute('SELECT r.name,r.length,r.gc_bin,r.purine_bin,a.start,a.end,a.target,a.kind,a.mapq FROM reads r LEFT JOIN alignments a ON r.name=a.name ORDER BY r.name,a.start,a.end')
        current = None

        def record_summary():
            if current is None:
                return
            name, length, gc_bin, purine_bin = current
            any_bp, primary_bp, high_bp, organelle_bp, other_bp = union
            overlap = organelle_bp+other_bp-primary_bp
            row = [name, length, gc_bin, purine_bin, *union, *counts, overlap]
            writer.writerow(row)
            values = dict(read_count=1, total_bases=length, any_mapped_reads=int(any_bp > 0),
                          primary_mapped_reads=int(primary_bp > 0), primary_mapq20_reads=int(high_bp > 0),
                          any_alignment_span_bp=any_bp, primary_span_bp=primary_bp,
                          primary_mapq20_span_bp=high_bp, organelle_primary_reads=int(counts[2] > 0),
                          primary_organelle_span_bp=organelle_bp, primary_other_reference_span_bp=other_bp,
                          primary_compartment_overlap_bp=overlap)
            for key, value in values.items():
                totals[key] += value
                groups[('gc', gc_bin)][key] += value
                groups[('max_purine_pyrimidine', purine_bin)][key] += value

        with (outdir/'read_mapping.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(fields)
            for name, length, gc_bin, purine_bin, start, end, target, kind, mapq in query:
                if current is None or name != current[0]:
                    record_summary()
                    current = (name, length, gc_bin, purine_bin)
                    union, last, counts = [0]*5, [0]*5, [0]*3
                if start is None:
                    continue
                primary = kind in {'P', 'I'}
                for axis, included in enumerate((True, primary, primary and 20 <= mapq < 255,
                                                primary and target in organelles,
                                                primary and target not in organelles)):
                    if included:
                        union[axis] += max(0, end-max(start, last[axis]))
                        last[axis] = max(last[axis], end)
                counts[0 if primary else 1] += 1
                counts[2] += primary and target in organelles
            record_summary()
        if not totals['read_count']:
            raise ValueError('No indexed query reads')
        composition = [dict(dimension=dimension, bin_start_percent=5*value, bin_end_percent=5*(value+1), **counts)
                       for (dimension, value), counts in sorted(groups.items())]
        write_table(outdir/'composition_mapping.tsv', composition, list(composition[0]))
        targets = []
        for target, length in reference_lengths.items():
            sums, union, last = [0, 0], [0, 0], [0, 0]
            for a, b, mapq in connection.execute('SELECT start,end,mapq FROM target_blocks WHERE target=? ORDER BY start,end', (target,)):
                for axis, included in enumerate((True, 20 <= mapq < 255)):
                    if included:
                        sums[axis] += b-a
                        union[axis] += max(0, b-max(a, last[axis]))
                        last[axis] = max(last[axis], b)
            targets.append(dict(contig=target, reference_length_bp=length, compartment='organelle' if target in organelles else 'other_reference',
                primary_aligned_target_bases=sums[0], primary_reference_union_bp=union[0], primary_aligned_base_depth=sums[0]/length,
                mapq20_aligned_target_bases=sums[1], mapq20_reference_union_bp=union[1], mapq20_aligned_base_depth=sums[1]/length))
        write_table(outdir/'reference_mapping.tsv', targets, list(targets[0]))
        return dict(totals, any_mapped_read_fraction=totals['any_mapped_reads']/totals['read_count'],
                    primary_mapq20_read_fraction=totals['primary_mapq20_reads']/totals['read_count'],
                    primary_query_span_fraction=totals['primary_span_bp']/totals['total_bases'],
                    warning='alignment_concordance_not_accuracy_truth;unmapped_is_not_proven_contamination;mapq20_not_unique_copy_proof;secondary_output_is_capped;composition_association_not_empirical_sequencing_bias')
    finally:
        connection.close()
