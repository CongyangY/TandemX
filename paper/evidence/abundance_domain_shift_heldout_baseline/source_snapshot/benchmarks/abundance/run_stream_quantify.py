"""Measure conditional copy-number errors on completed factorial source genomes."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

from benchmarks.abundance.evaluate import score_copy_number
from benchmarks.abundance.run import aggregate
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


def run(datasets: list[Path], outdir: Path, timeout: float = 900) -> None:
    if not datasets or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Require completed factorial datasets and a positive timeout')
    datasets = [p.resolve() for p in datasets]
    generations = [json.loads((p/'generation_receipt.json').read_text()) for p in datasets]
    if (len({r['seed'] for r in generations}) != len(generations)
            or any(not r['complete'] or r['split'] != 'development' or r['heldout_used'] or not r['conditions_completed'] for r in generations)):
        raise ValueError('Require unique completed development genomes, never held-out inputs')
    outdir = outdir.resolve(); outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    source = outdir/'source_snapshot'
    provenance = source_manifest(root, source)
    benchmark_hashes = {}
    for path in sorted(Path(__file__).parent.glob('*.py')):
        target = source/path.relative_to(root); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path,target); benchmark_hashes[str(path.relative_to(root))]=digest_file(target)
    provenance.update(benchmark_sha256=benchmark_hashes,
                      datasets={str(p): digest_file(p/'generation_receipt.json') for p in datasets},
                      scope='known_founder_catalogue_conditional_quantification_not_de_novo_recall',
                      estimator_inputs='observed_reads_founder_catalogue_fixed_genome_size_only',
                      resources='diagnostic_during_data_acquisition_and_other_development_jobs_not_publication_ranking',
                      interval_interpretation='native_kmer_spread_not_a_calibrated_sampling_confidence_interval')
    (outdir/'environment.json').write_text(json.dumps(provenance,indent=2)+'\n')
    metrics, executions = [], []
    for directory, generation in zip(datasets, generations):
        genomic = directory/'genome'
        genome = json.loads((genomic/'manifest.json').read_text())
        for name in ('catalogue.fa','truth_copy_number.tsv'):
            if digest_file(genomic/name) != genome['files'][name]:
                raise ValueError('Genomic truth/catalogue hash mismatch')
        truth = read_table(genomic/'truth_copy_number.tsv')
        factors = {r['family_id']: r for r in truth}
        for condition in generation['conditions_completed']:
            reads_dir = directory/'reads'/condition['condition_id']
            if digest_file(reads_dir/'manifest.json') != condition['manifest_sha256']:
                raise ValueError('Read-condition receipt hash mismatch')
            sampling = json.loads((reads_dir/'manifest.json').read_text())
            reads = reads_dir/'reads.fa'
            if digest_file(reads) != sampling['files']['reads.fa']:
                raise ValueError('Read input hash mismatch')
            folder = outdir/'runs'/f's{generation["seed"]}'/condition['condition_id']
            folder.mkdir(parents=True)
            command = [sys.executable,'-m','tandemx.cli','quantify','--reads',str(reads),
                       '--catalog',str(genomic/'catalogue.fa'),'--genome-size',str(genome['genome_bp']),
                       '--k','21','--kmer-backend','rust','--no-progress','--outdir',str(folder/'output')]
            measured = run_process(command,folder/'stdout.log',folder/'stderr.log',timeout,
                                   {**os.environ,'PYTHONPATH':str(source)},source)
            execution = dict(seed=generation['seed'],condition_id=condition['condition_id'],command=command,
                             reads_sha256=sampling['files']['reads.fa'], catalogue_sha256=genome['files']['catalogue.fa'],**measured)
            (folder/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
            executions.append(execution)
            with (outdir/'run.log').open('a') as log:
                log.write(f's{generation["seed"]}/{condition["condition_id"]}\texit={measured["exit_code"]}\n')
            if measured['exit_code'] != 0 or measured['timed_out']:
                continue
            # The source-occupancy oracle is evaluation-only; no truth or error
            # probability is passed in the software command above.
            oracle = dict(sampled_repeat_bp={f:sampling['sampled_repeat_bp'].get(f,0) for f in factors},
                          actual_base_coverage=sampling['actual_source_coverage'])
            scores = score_copy_number(folder/'output/copy_number.tsv',truth,oracle)
            for score in scores:
                factor=factors[score['family_id']]
                metrics.append(dict(seed=generation['seed'],condition_id=condition['condition_id'],
                    coverage=condition['coverage'],error_model=condition['label'],
                    substitution_rate=condition['substitution_rate'],insertion_rate=condition['insertion_rate'],deletion_rate=condition['deletion_rate'],
                    unit_substitution_rate=float(factor['unit_substitution_rate']),period=int(factor['period']),
                    requested_gc=float(factor['requested_gc']), founder_gc=float(factor['founder_gc']),
                    array_scope='megabase' if int(factor['repeat_bp'])>=1_000_000 else 'factorial',
                    source_coverage=sampling['actual_source_coverage'],observed_coverage=sampling['actual_observed_coverage'],
                    read_count=sampling['read_count'],total_bases=sampling['total_bases'],**score))
            write_table(outdir/'copy_number_metrics.tsv',metrics,list(metrics[0]))
            summaries = aggregate(metrics,['coverage','error_model','unit_substitution_rate','array_scope'],'copy')
            write_table(outdir/'copy_number_summary.tsv',summaries,list(summaries[0]))
    validation = dict(complete=all(r['exit_code']==0 and not r['timed_out'] for r in executions),
                      executions=len(executions),successful=sum(r['exit_code']==0 and not r['timed_out'] for r in executions),
                      independent_simulated_genomes=len(datasets),copy_number_family_conditions=len(metrics),
                      scientific_acceptance='not_assumed_from_successful_execution')
    (outdir/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    if not validation['complete']:
        raise RuntimeError('Failed quantification stage; preserve executions and do not replace missing metrics by zeros')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets',type=Path,nargs='+',required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=900)
    args=parser.parse_args()
    run(args.datasets,args.outdir,args.timeout)
