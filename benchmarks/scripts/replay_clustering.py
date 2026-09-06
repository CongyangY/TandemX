"""Compare old/new clustering on identical serialized real candidate evidence."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import time

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table
from tandemx.discover.clustering import cluster_monomers
from tandemx.discover.mvp import CandidateRepeat


def run(previous_run: Path, outdir: Path) -> None:
    previous_run=previous_run.resolve();outdir=outdir.resolve()
    old=previous_run/'source_snapshot/tandemx/discover/clustering.py'
    root=Path(__file__).resolve().parents[2]
    old_environment=json.loads((previous_run/'environment.json').read_text())
    for relative in ('tandemx/discover/clustering.py','tandemx/discover/distance.py'):
        if digest_file(previous_run/'source_snapshot'/relative)!=old_environment['file_hashes'][relative]:
            raise ValueError('Changed original source snapshot')
    if digest_file(root/'tandemx/discover/distance.py')!=digest_file(previous_run/'source_snapshot/tandemx/discover/distance.py'):
        raise ValueError('Distance implementation changed; cannot isolate indexing')
    folder=previous_run/'tandemx/discover'
    sequences={name.split(';')[0].split('=',1)[1]:seq for name,seq in read_fasta(folder/'candidate_monomers.fa').items()}
    candidates=[]
    records=read_table(folder/'candidate_reads.tsv')
    if {r['candidate_id'] for r in records}!=sequences.keys() or len(records)!=len(sequences):
        raise ValueError('Candidate sequence/table IDs differ or repeat')
    for row in records:
        candidates.append(CandidateRepeat(row['read_id'],row['candidate_id'],sequences[row['candidate_id']],
            int(row['read_start']),int(row['read_end']),row['strand'],int(row['period_bp']),int(row['repeat_span_bp']),
            float(row['unit_count']),float(row['score']),row['low_complexity_flag']=='true',row['confidence'],row['warning']))
    outdir.mkdir(parents=True,exist_ok=False)
    provenance=source_manifest(root,outdir/'source_snapshot')
    target=outdir/'replay_clustering.py';shutil.copyfile(__file__,target)
    provenance.update(script_sha256=digest_file(target),previous_run=str(previous_run),
        baseline_clustering_sha256=digest_file(old),
        input_sha256={name:digest_file(folder/name) for name in ('candidate_reads.tsv','candidate_monomers.fa')},
        input_precision='serialized_rounded_float_values_not_the_original_live_candidate_objects',
        scope='old_and_new_algorithms_on_identical_serialized_candidates_not_full_pipeline_parity',
        resources='single_process_stage_timing_diagnostics;no_stage_peak_memory_measurement')
    (outdir/'environment.json').write_text(json.dumps(provenance,indent=2)+'\n')
    spec=importlib.util.spec_from_file_location('tandemx_baseline_clustering',old)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    outputs={};measured=[]
    # One diagnostic pair; isolated shuffled repeated timing remains future work.
    for name,function in [('baseline',module.cluster_monomers),('indexed',cluster_monomers)]:
        began=time.perf_counter();families,membership=function(candidates,1,.95,'rust')
        elapsed=time.perf_counter()-began
        payload=dict(families=[asdict(f) for f in families],membership=membership)
        serialized=json.dumps(payload,sort_keys=True,separators=(',',':')).encode()
        outputs[name]=payload
        measured.append(dict(method=name,runtime_seconds=elapsed,candidates=len(candidates),families=len(families),
                             output_sha256=hashlib.sha256(serialized).hexdigest()))
        (outdir/f'{name}.json').write_bytes(serialized)
        with (outdir/'run.log').open('a') as log:log.write(f'{name}\t{elapsed:.6f}\t{len(families)} families\n')
    validation=dict(complete=True,exact_family_and_membership_parity=outputs['baseline']==outputs['indexed'],
                    measurements=measured,scope=provenance['scope'])
    (outdir/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    if not validation['exact_family_and_membership_parity']:
        raise RuntimeError('Clustering parity failed; preserve both outputs')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-run',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args();run(args.previous_run,args.outdir)
