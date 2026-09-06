"""Replay development abundance inputs with the experimental read-cluster model."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import statistics
import sys

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from tandemx.discover.mvp import read_fasta
from tandemx.quantify.mvp import read_monomer_fasta, family_kmer_membership, monomer_kmer_counts
from tandemx.quantify.read_moments import collect_read_moments
from tandemx.utils.kmers import is_low_complexity_kmer


def worker(baseline: Path, outdir: Path) -> None:
    config = json.loads((baseline/'run_config.json').read_text())
    if json.loads((baseline/'environment.json').read_text())['split'] != 'development':
        raise ValueError('This diagnostic replay must not use held-out inputs')
    k = config['k']
    estimates, metrics, inputs = [], [], []
    for seed in config['seeds']['development']:
        folder = baseline/'genomes'/f's{seed}'
        genome = json.loads((folder/'manifest.json').read_text())
        for name in ('catalogue.fa', 'truth_copy_number.tsv'):
            if digest_file(folder/name) != genome['files'][name]:
                raise ValueError(f'Genomic input hash differs: {name}')
        catalogue = list(read_monomer_fasta(folder/'catalogue.fa'))
        membership = family_kmer_membership(catalogue, k)
        diagnostic = {m.family_id: {word: n for word, n in monomer_kmer_counts(m.sequence, k).items()
                                  if len(membership[word]) == 1 and not is_low_complexity_kmer(word)}
                      for m in catalogue}
        truth = {r['family_id']: r for r in read_table(folder/'truth_copy_number.tsv')}
        if len(diagnostic) != len(catalogue) or set(diagnostic) != set(truth):
            raise ValueError('Duplicate or mismatched catalogue/truth families')
        for coverage in config['coverages']:
            for error in config['substitution_rates']:
                reads_dir = baseline/'reads'/f's{seed}'/f'c{coverage}_e{error}'
                sampling = json.loads((reads_dir/'manifest.json').read_text())
                observed_hash = digest_file(reads_dir/'reads.fa')
                if observed_hash != sampling['files']['reads.fa']:
                    raise ValueError(f'Read hash differs: {reads_dir}')
                # No sampled coordinates, actual error rate or planted counts enter the estimator.
                moments = collect_read_moments((r.sequence for r in read_fasta(reads_dir/'reads.fa')), diagnostic, k)
                if moments.read_count != sampling['read_count'] or moments.total_bases != sampling['read_count']*sampling['read_length']:
                    raise ValueError('Read totals differ from generator manifest')
                context = dict(seed=seed, coverage=coverage, substitution_rate=error)
                inputs.append(dict(**context, reads_sha256=observed_hash, catalogue_sha256=genome['files']['catalogue.fa'],
                                   truth_sha256=genome['files']['truth_copy_number.tsv'], sampling_sha256=digest_file(reads_dir/'manifest.json'),
                                   read_count=moments.read_count, total_bases=moments.total_bases, total_exposure=moments.sum_x))
                for estimate in moments.estimates(genome['genome_bp'], {f: len(d) for f, d in diagnostic.items()}):
                    estimates.append(dict(**context, **asdict(estimate)))
                    record = truth[estimate.family_id]
                    copies = int(record['copies'])
                    oracle = sampling['sampled_repeat_bp'][estimate.family_id]/int(record['period'])/sampling['actual_base_coverage']
                    value, lo, hi = estimate.estimated_copy_number, estimate.sampling_interval_low, estimate.sampling_interval_high
                    metrics.append(dict(**context, family_id=estimate.family_id, truth_copies=copies,
                                        estimate=value, signed_relative_error=(value-copies)/copies if value is not None else None,
                                        absolute_relative_error=abs(value-copies)/copies if value is not None else None,
                                        sampling_oracle_copy_estimate=oracle,
                                        estimator_minus_sampling_oracle=(value-oracle)/copies if value is not None else None,
                                        interval_available=lo is not None, interval_contains_truth=lo<=copies<=hi if lo is not None else None,
                                        interval_relative_width=(hi-lo)/copies if lo is not None else None,
                                        interval_status=estimate.interval_status))
                with (outdir/'worker.log').open('a') as log:
                    log.write(f'{context}\t{moments.read_count} reads\n')
    for name, rows in (('estimates.tsv', estimates), ('metrics.tsv', metrics), ('inputs.tsv', inputs),
                       ('summary.tsv', summarize(metrics))):
        write_table(outdir/name, rows, list(rows[0]))


def summarize(rows: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in rows:
        groups[(row['coverage'], row['substitution_rate'])].append(row)
    summaries = []
    for (coverage, error), group in sorted(groups.items()):
        observed = [r for r in group if r['estimate'] is not None]
        intervals = [r for r in group if r['interval_available']]
        summaries.append(dict(coverage=coverage, substitution_rate=error, total_family_conditions=len(group),
                              estimates_available=len(observed), intervals_available=len(intervals),
                              intervals_unavailable=len(group)-len(intervals),
                              mean_signed_relative_error=statistics.mean(r['signed_relative_error'] for r in observed) if observed else None,
                              mean_absolute_relative_error=statistics.mean(r['absolute_relative_error'] for r in observed) if observed else None,
                              mean_estimator_minus_oracle=statistics.mean(r['estimator_minus_sampling_oracle'] for r in observed) if observed else None,
                              conditional_interval_coverage=statistics.mean(r['interval_contains_truth'] for r in intervals) if intervals else None,
                              mean_available_interval_relative_width=statistics.mean(r['interval_relative_width'] for r in intervals) if intervals else None))
    return summaries


def run(baseline: Path, outdir: Path) -> None:
    baseline, outdir = baseline.resolve(), outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir/'source_snapshot'
    source = source_manifest(root, snapshot)
    target = snapshot/Path(__file__).relative_to(root)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(__file__, target)
    source.update(script_sha256=digest_file(target), baseline=str(baseline),
                  baseline_environment_sha256=digest_file(baseline/'environment.json'),
                  baseline_config_sha256=digest_file(baseline/'run_config.json'),
                  scientific_status='experimental Python estimator, conditional known-catalogue development only',
                  resource_note='diagnostic replay during data acquisition; not a publication speed comparison')
    (outdir/'environment.json').write_text(json.dumps(source, indent=2)+'\n')
    command = [sys.executable, str(target), '--baseline', str(baseline), '--outdir', str(outdir), '--worker']
    measured = run_process(command, outdir/'stdout.log', outdir/'stderr.log', 1800,
                           {**os.environ, 'PYTHONPATH': str(snapshot)}, snapshot)
    (outdir/'execution.json').write_text(json.dumps(dict(command=command, **measured), indent=2)+'\n')
    if measured['exit_code'] != 0 or measured['timed_out']:
        raise RuntimeError('Read-cluster diagnostic replay failed; inspect preserved logs')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    (worker if args.worker else run)(args.baseline, args.outdir)
