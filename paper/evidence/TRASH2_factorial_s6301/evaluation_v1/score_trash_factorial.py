"""Score completed TRASH assembly runs against independently planted genomes.

This is a development evaluation, not a comparison to HiFi read-based runs.
Native region, actual monomer coverage, coordinate and period endpoints remain
separate. No founder template or truth is supplied to comparator execution.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil

from benchmarks.abundance.stream_genome import FixedFastaReader
from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays, score_base_coverage
from benchmarks.challenge.run import json_safe, source_manifest
from benchmarks.challenge.schema import ArrayRecord, digest_file, read_table, write_table
from benchmarks.challenge.sequence_metrics import score_threshold_recovery
from benchmarks.challenge.trash_adapters import ARRAY_FILES, UNIT_FILES, csv_records, iter_regions, iter_units


def bounded_records(records):
    result = []
    for record in records:
        if len(result) >= 100_000:
            raise ValueError('Development evaluator is limited to 100,000 interval records')
        result.append(record)
    return result


def unit_coordinate_audit(genome: Path, native: Path, tool: str) -> list[dict]:
    """Exact local sequence checks using bounded FASTA seeks, not planted truth."""
    name = UNIT_FILES['trash'] if tool == 'trash' else 'assembly.fa_repeats_with_seq.csv'
    sequence_key, id_key = ('seq', 'seq.name') if tool == 'trash' else ('sequence', 'seqID')
    index = json.loads((genome/'genome_index.json').read_text())
    rows = [dict(offset_from_native_1based=i, monomer_rows=0, width_discrepancies=0,
                 strand_adjusted_matches=0, outside_reference=0) for i in range(-2, 3)]
    with (genome/'genome.fa').open('rb') as handle:
        fasta = FixedFastaReader(handle, index)
        for record in csv_records(native/name, {sequence_key, id_key, 'start', 'end', 'width', 'strand'}):
            a, b, width = (int(record[k]) for k in ('start', 'end', 'width'))
            sequence = record[sequence_key].upper()
            if (record[id_key] != index['contig'] or record['strand'] not in {'+', '-'}
                    or a < 1 or b < a or not sequence or set(sequence)-set('ACGTN')):
                raise ValueError('Invalid native monomer sequence/coordinate record')
            for row in rows:
                row['monomer_rows'] += 1
                row['width_discrepancies'] += int(width != b-a+1 or width != len(sequence))
                start, end = a-1+row['offset_from_native_1based'], b+row['offset_from_native_1based']
                if not 0 <= start < end <= index['length']:
                    row['outside_reference'] += 1
                    continue
                fragment = fasta.get(start, end)
                if record['strand'] == '-':
                    fragment = fragment.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                row['strand_adjusted_matches'] += int(sequence == fragment)
    return rows


def evaluate(genome: Path, run: Path, outdir: Path) -> dict:
    genome, run, outdir = genome.resolve(), run.resolve(), outdir.resolve()
    manifest = json.loads((genome/'manifest.json').read_text())
    execution = json.loads((run/'execution.json').read_text())
    tool = execution['tool']
    if (tool not in ARRAY_FILES or not execution['complete']
            or manifest.get('generator') != 'streamed_factorial_genome_v1'):
        raise ValueError('Require successful TRASH execution and independent factorial genome')
    index = json.loads((genome/'genome_index.json').read_text())
    if index['length'] != manifest['genome_bp']:
        raise ValueError('Genome index/manifest length disagreement')
    for name in ('catalogue.fa', 'truth_copy_number.tsv'):
        if (genome/name).stat().st_size > 10_000_000:
            raise ValueError('Development truth/catalogue exceeds explicit 10 MB limit')
    inputs = {genome/name: value for name, value in manifest['files'].items()}
    inputs.update({run/'native'/name: value for name, value in execution['native_products'].items()})
    inputs[run/'execution.json'] = digest_file(run/'execution.json')
    inputs[genome/'manifest.json'] = digest_file(genome/'manifest.json')
    random_access_source = Path(__file__).resolve().parents[1]/'abundance/stream_genome.py'
    inputs[random_access_source] = digest_file(random_access_source)
    unit_name = UNIT_FILES[tool] if tool == 'trash' else 'assembly.fa_repeats_with_seq.csv'
    inputs[run/'native'/unit_name] = digest_file(run/'native'/unit_name)
    if execution['input_sha256'] != manifest['files']['genome.fa']:
        raise ValueError('Comparator input differs from generated genome')
    if any(digest_file(path) != expected for path, expected in inputs.items()):
        raise ValueError('Comparator native output or generated truth changed')
    lengths = {index['contig']: index['length']}
    founders = read_fasta(genome/'catalogue.fa')
    truth_rows = read_table(genome/'truth_copy_number.tsv')
    truth = [ArrayRecord(r['chrom'], int(r['start']), int(r['end']), int(r['period']),
                        founders[r['family_id']], r['family_id']) for r in truth_rows]
    if len(truth) != len(founders) or len(founders) != len({r.family_id for r in truth}):
        raise ValueError('Expected one independent planted array per founder')
    outdir.mkdir(parents=True, exist_ok=False)
    receipt = dict(complete=False, tool=tool, seed=manifest['seed'], genome_bp=index['length'],
                   input_hashes={str(p): h for p, h in inputs.items()},
                   source=source_manifest(Path(__file__).resolve().parents[2]),
                   helper_sha256=digest_file(Path(__file__)),
                   parameters=dict(min_period=30, max_period=1000, min_span=100, min_iou=.5, sequence_threshold=.9),
                   warning='development_IID_haploid_genome;one_seed;assembly_only;not_comparable_to_HiFi_read_timings',
                   coordinate_note='all native CSVs unchanged;TRASH1 window-grid and R-extraction alternatives;monomer offsets independently audited')
    shutil.copyfile(Path(__file__), outdir/Path(__file__).name)
    try:
        audit = unit_coordinate_audit(genome, run/'native', tool)
        write_table(outdir/'unit_coordinate_audit.tsv', audit, list(audit[0]))
        rows = []
        coordinates = ('window_grid', 'r_extraction') if tool == 'trash' else ('one_based',)
        for policy in coordinates:
            for period in ('native_peak', 'consensus_length'):
                all_regions = bounded_records(iter_regions(run/'native', tool, policy, period))
                regions = [r for r in all_regions if 30 <= r.period <= 1000 and r.end-r.start >= 100]
                metrics, details = score_arrays(regions, truth, lengths)
                # One chromosome is not an independent read-level specificity experiment.
                metrics = {k: v for k, v in metrics.items() if 'read' not in k}
                family_metrics, family_details = score_threshold_recovery([r.sequence for r in regions], founders, .9)
                metrics.update(family_metrics, raw_region_count=len(all_regions), excluded_scope_count=len(all_regions)-len(regions))
                name = policy+'__'+period
                folder = outdir/name
                folder.mkdir()
                write_table(folder/'normalized_arrays.tsv', [asdict(r) for r in regions],
                            ['read_id', 'start', 'end', 'period', 'sequence', 'family_id'])
                write_table(folder/'array_details.tsv', details,
                            ['prediction_index', 'read_id', 'start', 'end', 'period', 'matched_truth_index', 'truth_family_id', 'iou', 'status'])
                write_table(folder/'family_recovery.tsv', family_details,
                            ['truth_id', 'recovered', 'assigned_sequence_index', 'threshold', 'criterion'])
                rows.append(dict(coordinate_policy=policy, period_source=period, **metrics))
        fields = list(rows[0])
        write_table(outdir/'region_metrics.tsv', rows, fields)
        units = []
        for offset in ((0, -1) if tool == 'trash' else (0,)):
            row = dict(offset_from_native_1based=offset, status='unavailable', warning='unit_coverage_only_not_array_period_or_recall')
            try:
                records = bounded_records(iter_units(run/'native', tool, offset))
                if any(r.read_id not in lengths or r.end > lengths[r.read_id] for r in records):
                    raise ValueError('Native unit coordinate exceeds reference')
                row.update(score_base_coverage(records, truth), unit_count=len(records), status='ok')
            except ValueError as exc:
                # A nonrepresentable sensitivity policy is missing, never zero
                # accuracy. The other policy and unchanged native rows remain.
                row['warning'] += ';'+str(exc)
            units.append(row)
        fields = list(dict.fromkeys(k for row in units for k in row))
        write_table(outdir/'unit_coverage_metrics.tsv', ({k: r.get(k) for k in fields} for r in units), fields)
        if not any(row['status'] == 'ok' for row in units):
            raise ValueError('No valid native monomer coverage endpoint')
        if any(digest_file(path) != expected for path, expected in inputs.items()):
            raise ValueError('Inputs changed during evaluation')
        receipt.update(complete=True, region_metrics=rows, unit_coverage_metrics=units,
                       outputs={str(p.relative_to(outdir)): digest_file(p) for p in outdir.rglob('*.tsv')})
        return receipt
    except Exception as exc:
        receipt['error'] = str(exc)
        raise
    finally:
        (outdir/'evaluation.json').write_text(json.dumps(json_safe(receipt), indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    evaluate(args.genome, args.run, args.outdir)
