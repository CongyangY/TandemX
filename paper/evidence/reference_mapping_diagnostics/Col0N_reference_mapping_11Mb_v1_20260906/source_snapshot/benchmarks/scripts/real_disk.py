"""Disk-backed real-input validation and descriptive array normalization.

The retained SQLite index uses a bounded page cache. Neither read identifiers
nor predictions are accumulated into a Python list or per-read dictionary.
"""
from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path
import sqlite3

from benchmarks.challenge.adapters import iter_arrays
from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.fastq_stream import hashed_fastq, records


ARRAY_FIELDS = ['read_id', 'start', 'end', 'period', 'sequence', 'family_id']


def connect_index(path: Path, *, create: bool = False) -> sqlite3.Connection:
    if create == path.exists():
        raise ValueError('Require a new index for creation, or an existing index for normalization')
    connection = sqlite3.connect(path)
    connection.execute('PRAGMA cache_size=-16384')
    connection.execute('PRAGMA temp_store=FILE')
    connection.execute('PRAGMA mmap_size=0')
    return connection


def prepare_disk_input(receipt: Path, sample_id: str, output: Path, database: Path) -> dict:
    """Verify a complete selected FASTQ before publishing the common FASTA.

    Interrupted/invalid work retains a .partial FASTA and the incomplete index.
    The ready metadata row is written only after checksum and all totals agree.
    """
    sampling = json.loads(receipt.read_text())
    samples = [r for r in sampling['samples'] if r['sample_id'] == sample_id]
    if not sampling.get('complete') or len(samples) != 1 or samples[0]['status'] != 'ok':
        raise ValueError('Require one completed nonempty whole-library sample')
    sample = samples[0]
    path = receipt.parent/(sample_id+'.fastq.gz')
    partial = output.with_name(output.name+'.partial')
    if output.exists() or partial.exists():
        raise ValueError('Refusing to overwrite an existing input FASTA or partial')
    connection = connect_index(database, create=True)
    reads = bases = 0
    try:
        connection.execute('CREATE TABLE reads (read_id TEXT PRIMARY KEY, length INTEGER NOT NULL) WITHOUT ROWID')
        connection.execute('CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID')
        with hashed_fastq(path) as (handle, checksum), partial.open('xb') as fasta:
            for record in records(handle):
                identifier = record.identifier.decode('ascii')
                try:
                    connection.execute('INSERT INTO reads VALUES (?, ?)', (identifier, len(record.sequence)))
                except sqlite3.IntegrityError as exc:
                    raise ValueError(f'Duplicate read ID in selected input: {identifier}') from exc
                reads += 1
                bases += len(record.sequence)
                fasta.write(b'>'+record.header[1:]+b'\n'+record.sequence+b'\n')
                if reads % 10_000 == 0:
                    connection.commit()
        if (not reads or checksum.hexdigest() != sample['fastq_sha256']
                or reads != sample['read_count'] or bases != sample['total_bases']):
            raise ValueError('Selected FASTQ does not match its sampling receipt')
        result = dict(**sample, sampling_receipt_sha256=digest_file(receipt),
                      fasta_sha256=digest_file(partial))
        connection.execute('INSERT INTO metadata VALUES (?, ?)', ('ready', json.dumps(result, sort_keys=True)))
        connection.commit()
        partial.rename(output)
        return result
    finally:
        connection.close()


def normalize_disk(tool: str, native: Path, database: Path, output: Path,
                   min_period: int = 30, max_period: int = 1000, min_span: int = 100) -> dict:
    """Retain call order in TSV; compute union with a disk-sorted interval scan.

    Invalid coordinates, unknown reads and late native parse errors fail the
    entire normalization. A .partial file is not a successful normalized output.
    Output calls are descriptive and provide no accuracy truth.
    """
    partial = output.with_name(output.name+'.partial')
    if output.exists() or partial.exists():
        raise ValueError('Refusing to overwrite normalized output or partial')
    connection = connect_index(database)
    try:
        ready = connection.execute('SELECT value FROM metadata WHERE name=?', ('ready',)).fetchone()
        if ready is None:
            raise ValueError('Input index is not complete')
        sample = json.loads(ready[0])
        connection.execute('CREATE TABLE IF NOT EXISTS arrays '
                           '(tool TEXT, ordinal INTEGER, read_id TEXT, start INTEGER, end INTEGER, '
                           'PRIMARY KEY (tool, ordinal)) WITHOUT ROWID')
        connection.execute('CREATE INDEX IF NOT EXISTS arrays_sorted ON arrays(tool, read_id, start, end)')
        if connection.execute('SELECT 1 FROM arrays WHERE tool=? LIMIT 1', (tool,)).fetchone():
            raise ValueError('Tool already has normalized predictions in this index')
        count = 0
        with connection, partial.open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, ARRAY_FIELDS, delimiter='\t')
            writer.writeheader()
            for array in iter_arrays(tool, native, min_period, max_period, min_span):
                count += 1
                connection.execute('INSERT INTO arrays VALUES (?, ?, ?, ?, ?)',
                                   (tool, count, array.read_id, array.start, array.end))
                writer.writerow(asdict(array))
            invalid = connection.execute(
                'SELECT a.read_id, a.end, r.length FROM arrays a LEFT JOIN reads r '
                'ON a.read_id=r.read_id WHERE a.tool=? AND (r.read_id IS NULL OR a.end>r.length) LIMIT 1',
                (tool,)).fetchone()
            if invalid is not None:
                raise ValueError(f'Native prediction has an unknown read or invalid coordinates: {invalid}')
            covered = positive = 0
            previous, union_end = None, -1
            for read_id, start, end in connection.execute(
                    'SELECT read_id, start, end FROM arrays WHERE tool=? ORDER BY read_id, start, end', (tool,)):
                if read_id != previous:
                    positive += 1
                    previous, union_end = read_id, -1
                covered += max(0, end-max(start, union_end))
                union_end = max(union_end, end)
        partial.rename(output)
        return dict(observed_in_scope_calls=count, observed_positive_reads=positive,
                    observed_union_bp=covered, observed_union_base_fraction=covered/sample['total_bases'])
    finally:
        connection.close()
