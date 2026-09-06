"""Native TRASH CSVs, with explicit coordinate and period sensitivity policies.

TRASH1 initially creates zero-based inclusive windows, then extracts them with
R's one-based str_sub. Its exported monomers add region.start to local one-based
positions. These are different contracts; never silently apply a monomer offset
to array windows. TRASH2 arrays and monomers use one-based inclusive positions.
"""
from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path

from .schema import ArrayRecord


ARRAY_FILES = {'trash': 'Summary.of.repetitive.regions.assembly.fa.csv',
               'trash2': 'assembly.fa_arrays.csv'}
UNIT_FILES = {'trash': 'all.repeats.from.assembly.fa.csv',
              'trash2': 'assembly.fa_repeats.csv'}


def csv_records(path: Path, fields: set[str]) -> Iterator[dict[str, str]]:
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if not fields <= set(names) or len(names) != len(set(names)):
            raise ValueError(f'Unexpected native CSV header: {path}')
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f'Malformed native CSV row: {path}')
            yield row


def iter_regions(native: Path, tool: str, coordinates: str, period_source: str = 'native_peak') -> Iterator[ArrayRecord]:
    """Retain all native regions; period-range filtering belongs to evaluation."""
    if period_source not in {'native_peak', 'consensus_length'}:
        raise ValueError('Unsupported period source')
    if tool == 'trash':
        identifier, period, sequence = 'fasta.name', 'most.freq.value.N', 'consensus.primary'
        if coordinates not in {'window_grid', 'r_extraction'}:
            raise ValueError('TRASH1 requires window_grid or r_extraction coordinates')
    elif tool == 'trash2':
        identifier, period, sequence = 'seqID', 'top_N', 'representative'
        if coordinates != 'one_based':
            raise ValueError('TRASH2 requires one_based coordinates')
    else:
        raise ValueError('Unsupported TRASH version')
    fields = {identifier, period, sequence, 'start', 'end'}
    for row in csv_records(native/ARRAY_FILES[tool], fields):
        a, b = int(row['start']), int(row['end'])
        if a < (1 if tool == 'trash2' else 0) or b < a:
            raise ValueError('Invalid native array coordinates')
        if coordinates == 'window_grid':
            start, end = a, b + 1
        else:
            # In the native R extraction, start=0 selects the first base. This
            # explicitly preserves that special case, not arbitrary clipping.
            start, end = max(1, a) - 1, b
        consensus = row[sequence].upper()
        if consensus in {'', 'NA'}:
            raise ValueError('Unresolved native array consensus; do not omit its row')
        size = int(row[period]) if period_source == 'native_peak' else len(consensus)
        yield ArrayRecord(row[identifier], start, end, size, consensus)


def iter_units(native: Path, tool: str, offset: int = 0) -> Iterator[ArrayRecord]:
    """Offset from one-based interpretation, explicitly supplied by evaluator.

    This yields unit intervals for base-union coverage only. Width is stored as
    the record's period to satisfy the interval container; it is not an inferred
    array period and must not enter array or family recall calculations.
    """
    if tool not in UNIT_FILES or offset not in {-1, 0} or (tool == 'trash2' and offset != 0):
        raise ValueError('Unsupported native monomer coordinate policy')
    identifier = 'seq.name' if tool == 'trash' else 'seqID'
    for row in csv_records(native/UNIT_FILES[tool], {identifier, 'start', 'end', 'width', 'strand'}):
        a, b, width = (int(row[k]) for k in ('start', 'end', 'width'))
        if a < 1 or b < a or width != b-a+1 or row['strand'] not in {'+', '-'}:
            raise ValueError('Invalid native monomer width/strand/coordinates')
        yield ArrayRecord(row[identifier], a-1+offset, b+offset, width)


def iter_unit_extents(native: Path, period_source: str = 'native_peak') -> Iterator[ArrayRecord]:
    """TRASH2 arrays from its explicit unit-to-array IDs and outer unit bounds.

    Author documentation calls the array table approximate. Retain that table
    separately; this derives boundaries from the main repeat output, never from
    truth overlap. Gaps between units remain inside the array interval. An array
    without units is unresolved, not silently dropped from the denominator.
    """
    arrays = {}
    fields = {'seqID', 'array_num_ID'}
    raw_arrays = csv_records(native/ARRAY_FILES['trash2'], fields)
    for region, raw in zip(iter_regions(native, 'trash2', 'one_based', period_source), raw_arrays, strict=True):
        key = raw['seqID'], int(raw['array_num_ID'])
        if key in arrays or key[1] < 1 or len(arrays) >= 100_000:
            raise ValueError('Duplicate/invalid array ID or development array-count limit exceeded')
        arrays[key] = region
    extents = {}
    raw_units = csv_records(native/UNIT_FILES['trash2'], {'seqID', 'arrayID'})
    for unit, raw in zip(iter_units(native, 'trash2'), raw_units, strict=True):
        key = raw['seqID'], int(raw['arrayID'])
        if key not in arrays:
            raise ValueError('Native monomer points to an unknown array ID')
        if key in extents:
            start, end = extents[key]
            extents[key] = min(start, unit.start), max(end, unit.end)
        else:
            extents[key] = unit.start, unit.end
    if arrays.keys() != extents.keys():
        raise ValueError('Native array without observed monomers; do not omit its row')
    for key, region in arrays.items():
        yield ArrayRecord(region.read_id, *extents[key], region.period, region.sequence)
