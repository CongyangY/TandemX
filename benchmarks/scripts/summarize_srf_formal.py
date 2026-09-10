import csv,json,statistics,hashlib,shutil
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description='Aggregate and archive frozen SRF comparator outputs without modifying raw runs')
parser.add_argument('--results', type=Path, required=True)
parser.add_argument('--archive', type=Path, required=True)
args=parser.parse_args()
root=args.results
completion=json.loads((root/'completion.json').read_text())
assert completion['completed_cells']==72 and completion['unrun_cells']==0
rows=list(csv.DictReader((root/'summary.tsv').open(),delimiter='\t'))
keys=['MARE','family_recovery','base_union_recall','base_union_precision','negative_read_predicted_bp','negative_read_attribution_fraction','unassigned_abundance_bp','native_catalogue_count','workflow_wall_seconds','cpu_user_seconds','cpu_system_seconds','peak_rss_mib','throughput_bp_per_second','shared_discovery_plus_mapping_wall_seconds','shared_discovery_plus_mapping_peak_rss_mib']
agg=[]
for condition in dict.fromkeys(r['condition'] for r in rows):
 for method in dict.fromkeys(r['method'] for r in rows):
  sub=[r for r in rows if r['condition']==condition and r['method']==method]
  assert len(sub)==3
  item={'condition':condition,'method':method,'simulation_seeds':3,'statuses':';'.join(f"{s}:{sum(r['status']==s for r in sub)}" for s in sorted({r['status'] for r in sub}))}
  for key in keys:
   values=[float(r[key]) for r in sub if r.get(key) not in (None,'','NA')]
   for suffix,fn in [('mean',statistics.mean),('min',min),('max',max)]:
    item[f'{key}_{suffix}']=fn(values) if values else 'NA'
  agg.append(item)
out=args.archive
out.mkdir(parents=True,exist_ok=False)
for name in ['summary.tsv','completion.json','environment.json','clean_kmc_build_receipt.json','earlier_comparator_build_provenance.json']:
 shutil.copyfile(root/name,out/name)
for p in (root/'runs').glob('*/*/*'):
 if p.name in {'result.json','family_abundance.tsv','correspondence.json','interval_metrics.json','native_abundance.json','workflow_receipt.json'}:
  dst=out/p.relative_to(root);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
for p in (root/'datasets').glob('*/*'):
 if p.name in {'manifest.json','truth_arrays.tsv','truth_reads.tsv','truth_monomers.fa'}:
  dst=out/p.relative_to(root);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
with (out/'condition_summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(agg[0]),delimiter='\t');w.writeheader();w.writerows(agg)
manifest={str(p.relative_to(out)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in out.rglob('*') if p.is_file()}
(out/'archive_manifest.json').write_text(json.dumps({'source':str(root),'protocol_commit':'7abafc3','files':manifest},indent=2)+'\n')
print(json.dumps(completion))
for r in agg:
 print(r['condition'],r['method'],'MARE',round(r['MARE_mean'],6),'rec',r['family_recovery_mean'],'neg',r['negative_read_predicted_bp_mean'],'unassigned',r['unassigned_abundance_bp_mean'],'wall',round(r['workflow_wall_seconds_mean'],3),'RSS',round(r['peak_rss_mib_mean'],2))
