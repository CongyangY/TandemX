"""Calibrate joint-read intervals on the frozen development factorial replay."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import sys

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from tandemx.discover.mvp import read_fasta
from tandemx.quantify.joint_multik import estimate_joint_multik
from tandemx.quantify.multik import DEFAULT_K_VALUES
from tandemx.quantify.mvp import read_monomer_fasta


def summarize(rows: list[dict]) -> list[dict]:
    # Report each independent genome seed separately. Families and nested read
    # conditions are dependent; do not produce a binomial CI from pooled rows.
    fields = ('seed', 'coverage', 'error_model', 'unit_substitution_rate', 'array_scope')
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in fields)].append(row)
    output = []
    for key, group in sorted(groups.items()):
        available = [r for r in group if r['sampling_interval_low'] is not None]
        covered = sum(r['truth_covered'] for r in available)
        output.append(dict(zip(fields, key), family_conditions=len(group),
            intervals_available=len(available), truth_covered=covered,
            conditional_coverage_fraction=covered/len(available) if available else None,
            available_and_covering_fraction=covered/len(group),
            mean_relative_interval_width=statistics.mean(r['relative_interval_width'] for r in available) if available else None,
            warning='descriptive_within_genome_dependent_families_not_independent_binomial_trials'))
    return output


def worker(previous: Path, outdir: Path) -> None:
    prior_env = json.loads((previous/'environment.json').read_text())
    prior_valid = json.loads((previous/'validation.json').read_text())
    if not prior_valid['complete'] or prior_valid['heldout_used']:
        raise ValueError('Require complete development multi-k replay')
    baseline = Path(prior_env['baseline'])
    if any(digest_file(baseline/name) != value for name, value in prior_env['baseline_sha256'].items()):
        raise ValueError('Changed baseline')
    datasets = json.loads((baseline/'environment.json').read_text())['datasets']
    prior_rows = [r for r in read_table(previous/'metrics.tsv') if r['method'] == 'multik_loglinear']
    by_key = {(int(r['seed']), r['condition_id'], r['family_id']): r for r in prior_rows}
    if len(by_key) != len(prior_rows):
        raise ValueError('Duplicate previous family conditions')
    scored, rows, points, inputs = set(), [], [], []
    for directory, receipt_hash in datasets.items():
        directory = Path(directory)
        if digest_file(directory/'generation_receipt.json') != receipt_hash:
            raise ValueError('Changed generation receipt')
        generation = json.loads((directory/'generation_receipt.json').read_text())
        if generation['split'] != 'development' or generation['heldout_used'] or not generation['complete']:
            raise ValueError('Only completed development data permitted')
        genome = json.loads((directory/'genome/manifest.json').read_text())
        catalogue_path = directory/'genome/catalogue.fa'
        if digest_file(catalogue_path) != genome['files']['catalogue.fa']:
            raise ValueError('Changed founder catalogue')
        catalogue = list(read_monomer_fasta(catalogue_path))
        for condition in generation['conditions_completed']:
            folder = directory/'reads'/condition['condition_id']
            if digest_file(folder/'manifest.json') != condition['manifest_sha256']:
                raise ValueError('Changed read manifest')
            sampling = json.loads((folder/'manifest.json').read_text())
            if digest_file(folder/'reads.fa') != sampling['files']['reads.fa']:
                raise ValueError('Changed reads')
            # Only observed reads, supplied catalogue, fixed G/k/min-support
            # reach the estimator. Truth is used below, after inference.
            result = estimate_joint_multik((r.sequence for r in read_fasta(folder/'reads.fa')),
                catalogue, genome['genome_bp'], DEFAULT_K_VALUES, backend='rust', minimum_effective_reads=20)
            if result.read_count != sampling['read_count'] or result.total_bases != sampling['total_bases']:
                raise ValueError('Observed totals differ from receipt')
            context = dict(seed=generation['seed'], condition_id=condition['condition_id'])
            inputs.append(dict(**context, reads_sha256=sampling['files']['reads.fa'],
                catalogue_sha256=genome['files']['catalogue.fa'], read_count=result.read_count,
                total_bases=result.total_bases, warning=result.warning))
            points.extend(dict(**context, **p) for p in result.per_k)
            for estimate in result.estimates:
                key = (generation['seed'], condition['condition_id'], estimate.family_id)
                if key in scored or key not in by_key:
                    raise ValueError('Repeated or unexpected family condition')
                scored.add(key)
                prior = by_key[key]
                value = estimate.estimated_copy_number
                old_value = None if prior['estimate'] == 'NA' else float(prior['estimate'])
                if ((value is None) != (old_value is None)
                        or value is not None and not math.isclose(value, old_value, rel_tol=1e-8, abs_tol=1e-10)):
                    raise ValueError('Joint collector changed previous point estimate beyond rounding tolerance')
                truth = float(prior['truth_copies'])
                lo, hi = estimate.sampling_interval_low, estimate.sampling_interval_high
                row = asdict(estimate)
                row['positive_reads_by_k'] = ','.join(map(str, estimate.positive_reads_by_k))
                rows.append(dict(**context, **row,
                    **{k: prior[k] for k in ('coverage', 'error_model', 'unit_substitution_rate', 'array_scope', 'period')},
                    truth_copies=truth, truth_covered=lo <= truth <= hi if lo is not None else None,
                    relative_interval_width=(hi-lo)/truth if lo is not None else None,
                    signed_relative_error=(value-truth)/truth if value is not None else None))
            for name, data in [('metrics.tsv', rows), ('per_k.tsv', points), ('inputs.tsv', inputs),
                               ('summary.tsv', summarize(rows))]:
                write_table(outdir/name, data, list(data[0]))
            with (outdir/'worker.log').open('a') as log:
                log.write(f's{generation["seed"]}/{condition["condition_id"]}\t{result.read_count} reads\n')
    if scored != set(by_key):
        raise ValueError('Omitted previous family conditions')
    (outdir/'validation.json').write_text(json.dumps(dict(complete=True, input_conditions=len(inputs),
        family_conditions=len(rows), point_estimates_agree_previous=True, k_values=list(DEFAULT_K_VALUES),
        minimum_effective_reads=20, heldout_used=False,
        scientific_acceptance='requires_stratified_coverage_review_not_assumed'), indent=2)+'\n')


def run(previous: Path, outdir: Path) -> None:
    previous, outdir = previous.resolve(), outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir/'source_snapshot'
    provenance = source_manifest(root, snapshot)
    target = snapshot/Path(__file__).relative_to(root)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(__file__, target)
    hashes = {name: digest_file(previous/name) for name in ('environment.json', 'validation.json', 'metrics.tsv')}
    provenance.update(script_sha256=digest_file(target), previous=str(previous), previous_sha256=hashes,
        scope='development_known_founders_conditional_sampling_uncertainty',
        resource_note='diagnostic_during_other_jobs_not_isolated_resource_ranking')
    (outdir/'environment.json').write_text(json.dumps(provenance, indent=2)+'\n')
    command = [sys.executable, str(target), '--previous', str(previous), '--outdir', str(outdir), '--worker']
    measured = run_process(command, outdir/'stdout.log', outdir/'stderr.log', 3600,
                           {**os.environ, 'PYTHONPATH': str(snapshot)}, snapshot)
    (outdir/'execution.json').write_text(json.dumps(dict(command=command, **measured), indent=2)+'\n')
    if measured['exit_code'] != 0 or measured['timed_out']:
        raise RuntimeError('Joint multi-k calibration failed; retain logs and missing values')
    if any(digest_file(previous/name) != value for name, value in hashes.items()):
        raise RuntimeError('Previous results changed during calibration')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    (worker if args.worker else run)(args.previous, args.outdir)
