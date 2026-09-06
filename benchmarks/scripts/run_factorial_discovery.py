"""Run three discovery methods on completed factorial observed-read conditions."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import random
import shutil
import sys
import time

from benchmarks.abundance.discovery_truth import load_truth, score_predictions
from benchmarks.challenge.adapters import build_command, parse_arrays, read_fasta
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.challenge.sequence_metrics import score_threshold_recovery


def run(datasets: list[Path], conditions: list[str], outdir: Path, trf: Path, tidehunter: Path,
        timeout: float=900) -> None:
    if (not datasets or not conditions or len(set(conditions))!=len(conditions)
            or len({p.resolve() for p in datasets})!=len(datasets)
            or not math.isfinite(timeout) or timeout<=0):
        raise ValueError('Require unique development datasets/conditions and positive timeout')
    tools=dict(tandemx=sys.executable,trf=str(trf.resolve()),tidehunter=str(tidehunter.resolve()))
    if any(not Path(p).is_file() or not os.access(p,os.X_OK) for p in tools.values()):
        raise ValueError('Missing comparator executable')
    outdir=outdir.resolve();outdir.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[2];snapshot=outdir/'source_snapshot'
    provenance=source_manifest(root,snapshot)
    hashes={}
    for path in (Path(__file__),root/'benchmarks/abundance/discovery_truth.py'):
        target=snapshot/path.relative_to(root);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target);hashes[str(path.relative_to(root))]=digest_file(target)
    provenance.update(runner_source_sha256=hashes,tool_paths=tools,
        tool_hashes={t:digest_file(Path(p)) for t,p in tools.items()},conditions=conditions,
        datasets={str(p.resolve()):digest_file(p/'generation_receipt.json') for p in datasets},
        parameters=dict(min_period=30,max_period=1000,min_span=100,min_iou=.5,sequence_threshold=.9,
                        threads=1,repetitions=1,timeout=timeout),
        truth_policy='all_observed_planted_bases;array_endpoint_excludes_reads_with_ineligible_truth_fragments',
        sequence_policy='all_per_array_consensuses;all_planted_and_observed_eligible_denominators',
        resources='diagnostic_concurrent_jobs;direct_child_wait4_excludes_controller_and_evaluation')
    (outdir/'environment.json').write_text(json.dumps(provenance,indent=2)+'\n')
    summaries,inputs=[],[];seen_seeds=set()
    for dataset in datasets:
        dataset=dataset.resolve()
        seed=json.loads((dataset/'generation_receipt.json').read_text())['seed']
        if seed in seen_seeds:
            raise ValueError('Independent genome seeds must be unique')
        seen_seeds.add(seed)
        for condition in conditions:
            truth=load_truth(dataset,condition)
            inputs.append(truth.metadata)
            write_table(outdir/'inputs.tsv',inputs,list(inputs[0]))
            observed_founders={f:s for f,s in truth.founders.items() if f in {r.family_id for r in truth.eligible_segments}}
            order=list(tools);random.Random(f'{seed}:{condition}').shuffle(order)
            for tool in order:
                folder=outdir/'runs'/f's{seed}'/condition/tool;folder.mkdir(parents=True)
                command,output=build_command(tool,tools[tool],dataset/'reads'/condition/'reads.fa',folder,30,1000,100)
                if tool=='tandemx':
                    command[:1]=[sys.executable,'-m','tandemx.cli']
                    command+=['--discovery-method','elastic','--clustering-method','sequence','--cluster-identity','.95',
                              '--family-audit','related']
                measured=run_process(command,output if tool=='trf' else folder/'stdout.log',folder/'stderr.log',timeout,
                    {**os.environ,'PYTHONPATH':str(snapshot)} if tool=='tandemx' else None,
                    snapshot if tool=='tandemx' else folder)
                (folder/'execution.json').write_text(json.dumps(dict(command=command,**measured),indent=2)+'\n')
                row=dict(seed=seed,condition_id=condition,tool=tool,**measured,status='execution_failed',
                         evaluation_seconds=None,warning='development_truth;correlated_arrays_not_independent_replicates')
                if measured['exit_code']==0 and not measured['timed_out']:
                    began=time.perf_counter()
                    try:
                        arrays=parse_arrays(tool,output,30,1000,100)
                        metrics,details=score_predictions(arrays,truth)
                        sequences=(list(read_fasta(folder/'discover/candidate_monomers.fa').values()) if tool=='tandemx'
                                   else [r.sequence for r in arrays])
                        # Candidate consensus records represent the same interval scope.
                        if tool=='tandemx' and len(sequences)!=len(arrays):
                            raise ValueError('TandemX candidate consensus/interval scope differs')
                        for prefix,founders in [('all_planted',truth.founders),('observed_eligible',observed_founders)]:
                            family_metrics,family_details=score_threshold_recovery(sequences,founders,.9)
                            metrics.update({f'{prefix}_{k}':v for k,v in family_metrics.items()})
                            if family_details:
                                write_table(folder/f'{prefix}_family_recovery.tsv',family_details,list(family_details[0]))
                        write_table(folder/'normalized_arrays.tsv',[asdict(a) for a in arrays],
                                    ['read_id','start','end','period','sequence','family_id'])
                        if details:
                            write_table(folder/'eligible_array_details.tsv',details,list(details[0]))
                        json_metrics={k:None if isinstance(v,float) and not math.isfinite(v) else v for k,v in metrics.items()}
                        (folder/'metrics.json').write_text(json.dumps(json_metrics,indent=2,allow_nan=False)+'\n')
                        row.update(metrics,status='ok')
                    except Exception as exc:
                        row['status']='evaluation_failed';row['warning']+=';'+str(exc)
                    row['evaluation_seconds']=time.perf_counter()-began
                summaries.append(row)
                fields=list(dict.fromkeys(k for r in summaries for k in r))
                write_table(outdir/'summary.tsv',({k:r.get(k) for k in fields} for r in summaries),fields)
                with (outdir/'run.log').open('a') as log:
                    log.write(f's{seed}/{condition}/{tool}\t{row["status"]}\n')
    complete=all(r['status']=='ok' for r in summaries)
    (outdir/'validation.json').write_text(json.dumps(dict(complete=complete,executions=len(summaries),
        successful=sum(r['status']=='ok' for r in summaries),independent_source_genomes=len(seen_seeds),
        scientific_acceptance='not_inferred_from_success'),indent=2)+'\n')
    if not complete:
        raise RuntimeError('Failed discovery stage(s); retain missing metrics and inspect logs')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets',type=Path,nargs='+',required=True)
    parser.add_argument('--conditions',nargs='+',required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--trf',type=Path,required=True)
    parser.add_argument('--tidehunter',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=900)
    args=parser.parse_args()
    run(args.datasets,args.conditions,args.outdir,args.trf,args.tidehunter,args.timeout)
