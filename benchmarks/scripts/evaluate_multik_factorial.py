"""Paired development replay of fixed multi-k model and two single-k baselines."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import statistics
import sys

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from tandemx.discover.mvp import read_fasta
from tandemx.quantify.mvp import read_monomer_fasta
from tandemx.quantify.multik import DEFAULT_K_VALUES, estimate_multik


def summarize(metrics: list[dict]) -> list[dict]:
    fields=('method','coverage','error_model','unit_substitution_rate','array_scope')
    groups=defaultdict(list)
    for row in metrics:
        groups[tuple(row[k] for k in fields)].append(row)
    summaries=[]
    for key,group in sorted(groups.items()):
        valid=[r for r in group if r['estimate'] is not None]
        summaries.append(dict(zip(fields,key), family_conditions=len(group), estimates_available=len(valid),
            mean_signed_relative_error=statistics.mean(r['signed_relative_error'] for r in valid) if valid else None,
            mean_absolute_relative_error=statistics.mean(r['absolute_relative_error'] for r in valid) if valid else None,
            median_absolute_relative_error=statistics.median(r['absolute_relative_error'] for r in valid) if valid else None,
            mean_estimator_minus_oracle=statistics.mean(r['estimator_minus_sampling_oracle'] for r in valid) if valid else None))
    return summaries


def worker(baseline: Path, outdir: Path) -> None:
    environment=json.loads((baseline/'environment.json').read_text())
    valid=json.loads((baseline/'validation.json').read_text())
    if not valid['complete'] or valid['executions']!=valid['successful']:
        raise ValueError('Require a completed conditional baseline')
    base_rows=read_table(baseline/'copy_number_metrics.tsv')
    by_key={(int(r['seed']),r['condition_id'],r['family_id']):r for r in base_rows}
    if len(by_key)!=len(base_rows):
        raise ValueError('Duplicate baseline family conditions')
    outputs,points,metrics,inputs=[],[],[],[]
    scored_keys=set()
    for directory,receipt_hash in environment['datasets'].items():
        directory=Path(directory)
        if digest_file(directory/'generation_receipt.json')!=receipt_hash:
            raise ValueError('Changed generation receipt')
        generation=json.loads((directory/'generation_receipt.json').read_text())
        if generation['split']!='development' or generation['heldout_used'] or not generation['complete']:
            raise ValueError('Only completed development inputs permitted')
        genomic=directory/'genome'
        genome=json.loads((genomic/'manifest.json').read_text())
        if digest_file(genomic/'catalogue.fa')!=genome['files']['catalogue.fa']:
            raise ValueError('Changed catalogue')
        catalogue=list(read_monomer_fasta(genomic/'catalogue.fa'))
        for condition in generation['conditions_completed']:
            folder=directory/'reads'/condition['condition_id']
            if digest_file(folder/'manifest.json')!=condition['manifest_sha256']:
                raise ValueError('Changed condition manifest')
            sampling=json.loads((folder/'manifest.json').read_text())
            if digest_file(folder/'reads.fa')!=sampling['files']['reads.fa']:
                raise ValueError('Changed observed reads')
            # This function has no truth/error/coordinate inputs. Fixed k values
            # were chosen before running this experiment, not fit to its scores.
            result=estimate_multik((r.sequence for r in read_fasta(folder/'reads.fa')), catalogue,
                                   genome['genome_bp'], DEFAULT_K_VALUES, backend='rust')
            if result.read_count!=sampling['read_count'] or result.total_bases!=sampling['total_bases']:
                raise ValueError('Observed totals differ from generator receipt')
            context=dict(seed=generation['seed'],condition_id=condition['condition_id'])
            inputs.append(dict(**context,reads_sha256=sampling['files']['reads.fa'],
                               catalogue_sha256=genome['files']['catalogue.fa'],
                               read_count=result.read_count,total_bases=result.total_bases))
            outputs.extend(dict(**context,**r) for r in result.estimates)
            points.extend(dict(**context,**r) for r in result.per_k)
            for estimate in result.estimates:
                key=(generation['seed'],condition['condition_id'],estimate['family_id'])
                if key in scored_keys or key not in by_key:
                    raise ValueError('Repeated or unknown family-condition key')
                scored_keys.add(key)
                base=by_key[key]
                truth=float(base['truth_copies']); oracle=float(base['sampling_oracle_copy_estimate'])
                strata={k:base[k] for k in ('coverage','error_model','unit_substitution_rate','period',
                                          'requested_gc','founder_gc','array_scope')}
                for method,value in [('native_median_k21',float(base['estimate'])),
                                     ('mean_k21_exposure',estimate['uncorrected_mean_k21']),
                                     ('multik_loglinear',estimate['extrapolated_copy_number'])]:
                    metrics.append(dict(**context,**strata,family_id=estimate['family_id'],method=method,
                        truth_copies=truth,sampling_oracle_copy_estimate=oracle,estimate=value,
                        signed_relative_error=(value-truth)/truth if value is not None else None,
                        absolute_relative_error=abs(value-truth)/truth if value is not None else None,
                        estimator_minus_sampling_oracle=(value-oracle)/truth if value is not None else None,
                        status=estimate['status'] if method=='multik_loglinear' else 'single_k_baseline',
                        warning='development_only_no_confidence_interval_or_external_superiority'))
            for name,rows in [('estimates.tsv',outputs),('per_k.tsv',points),('metrics.tsv',metrics),('inputs.tsv',inputs),
                              ('summary.tsv',summarize(metrics))]:
                write_table(outdir/name,rows,list(rows[0]))
            with (outdir/'worker.log').open('a') as log:
                log.write(f's{generation["seed"]}/{condition["condition_id"]}\t{result.read_count} reads\n')
    if scored_keys!=set(by_key):
        raise ValueError('Replay omitted baseline family conditions')
    (outdir/'validation.json').write_text(json.dumps(dict(complete=True,input_conditions=len(inputs),
        family_conditions=len(outputs),paired_method_rows=len(metrics),k_values=list(DEFAULT_K_VALUES),
        heldout_used=False,scientific_acceptance='not_assumed_from_execution_success'),indent=2)+'\n')


def run(baseline: Path, outdir: Path) -> None:
    baseline,outdir=baseline.resolve(),outdir.resolve()
    outdir.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[2]
    snapshot=outdir/'source_snapshot'
    provenance=source_manifest(root,snapshot)
    target=snapshot/Path(__file__).relative_to(root)
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(__file__,target)
    hashes={name:digest_file(baseline/name) for name in ('environment.json','validation.json','copy_number_metrics.tsv')}
    provenance.update(script_sha256=digest_file(target),baseline=str(baseline),baseline_sha256=hashes,
        k_values=list(DEFAULT_K_VALUES),scope='conditional_known_founder_development_point_estimates_only',
        resource_note='diagnostic_during_other_jobs_not_final_resource_ranking')
    (outdir/'environment.json').write_text(json.dumps(provenance,indent=2)+'\n')
    command=[sys.executable,str(target),'--baseline',str(baseline),'--outdir',str(outdir),'--worker']
    measured=run_process(command,outdir/'stdout.log',outdir/'stderr.log',1800,
                         {**os.environ,'PYTHONPATH':str(snapshot)},snapshot)
    (outdir/'execution.json').write_text(json.dumps(dict(command=command,**measured),indent=2)+'\n')
    if measured['exit_code']!=0 or measured['timed_out']:
        raise RuntimeError('Multi-k replay failed; retain logs and missing values')
    if any(digest_file(baseline/name)!=value for name,value in hashes.items()):
        raise RuntimeError('Baseline changed during replay; reject this run')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    (worker if args.worker else run)(args.baseline,args.outdir)
