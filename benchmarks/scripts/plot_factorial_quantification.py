"""Six-panel source-backed diagnostic of the fixed factorial multi-k replay."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, read_table, write_table


METHODS=('native_median_k21','mean_k21_exposure','multik_loglinear')
LABELS=('Median k21','Mean k21','Multi-k prototype')
COLORS=('#687785','#d68d32','#187d91')
ERRORS=('error_free','iid_low_indel','iid_high_indel')


def run(replay: Path, baseline: Path, outdir: Path) -> None:
    validation=json.loads((replay/'validation.json').read_text())
    if not validation['complete'] or validation['heldout_used']:
        raise ValueError('Require completed development replay')
    rows=read_table(replay/'metrics.tsv');fits=read_table(replay/'estimates.tsv')
    original=read_table(baseline/'copy_number_metrics.tsv')
    if len(rows)!=validation['paired_method_rows'] or len(fits)!=validation['family_conditions']:
        raise ValueError('Metric cardinality differs from receipt')
    outdir.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
                         'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(3,2,figsize=(12.4,12),layout='constrained')
    sources=[]
    def source(panel,statistic,x,y,**context):
        sources.append(dict(panel=panel,statistic=statistic,x=x,y=y,**context))
    def key(r): return (r['seed'],r['condition_id'],r['family_id'])
    by_method={m:{key(r):r for r in rows if r['method']==m} for m in METHODS}
    fit_by_key={key(r):r for r in fits}
    if any(set(d)!=set(fit_by_key) for d in by_method.values()):
        raise ValueError('Methods have unequal pairing keys')
    ax=axs[0,0]
    strata=[(d,e) for d in (0,.02) for e in ERRORS]
    for method,label,color in zip(METHODS,LABELS,COLORS):
        means=[]
        for j,(d,e) in enumerate(strata):
            selected=[r for r in by_method[method].values() if r['coverage']=='20' and r['error_model']==e
                      and r['array_scope']=='factorial' and float(r['unit_substitution_rate'])==d]
            if any(r['estimate']=='NA' for r in selected):
                raise ValueError('Unexpected missing 20x estimate; revise displayed denominators')
            value=100*statistics.mean(float(r['signed_relative_error']) for r in selected)
            means.append(value);source('A','mean_signed_error_percent',j,value,method=method,coverage=20,error=e,divergence=d,n=len(selected))
            for seed in sorted({r['seed'] for r in selected}):
                v=100*statistics.mean(float(r['signed_relative_error']) for r in selected if r['seed']==seed)
                ax.scatter(j,v,color=color,s=12,alpha=.4)
                source('A','genome_mean_signed_error_percent',j,v,method=method,seed=seed,error=e,divergence=d)
        ax.plot(range(6),means,'o-',ms=4,lw=1.5,color=color,label=label)
    ax.axhline(0,color='black',lw=.7,ls='--')
    ax.set(xticks=range(6),xticklabels=['0 / clean','0 / low','0 / high','2 / clean','2 / low','2 / high'],
           xlabel='Unit divergence (%) / read error model',ylabel='Mean signed copy-number error (%)',
           title='A  Bias across error and divergence at 20×')
    ax.tick_params(axis='x',labelrotation=25);ax.legend(fontsize=8,loc='lower left')

    ax=axs[0,1];counts=[]
    for method,label,color in zip(METHODS,LABELS,COLORS):
        means=[]
        for coverage in (1,5,20):
            common=[k for k,r in by_method[method].items() if r['coverage']==str(coverage)
                    and r['error_model']=='iid_high_indel' and r['array_scope']=='factorial'
                    and float(r['unit_substitution_rate'])==.02
                    and all(by_method[m][k]['estimate']!='NA' for m in METHODS)]
            if method==METHODS[0]:counts.append(len(common))
            value=100*statistics.mean(float(by_method[method][k]['absolute_relative_error']) for k in common)
            means.append(value);source('B','complete_pair_mean_absolute_error_percent',coverage,value,method=method,n=len(common))
        ax.plot((1,5,20),means,'o-',ms=4,color=color,label=label)
    ax.set(xscale='log',xticks=(1,5,20),xticklabels=[f'{c}×\n{n} pairs' for c,n in zip((1,5,20),counts)],
           ylabel='Mean absolute relative error (%)',xlabel='Source coverage; complete pairs only',
           title='B  Higher-error reads with 2% unit divergence')

    ax=axs[1,0];outcomes=[]
    for coverage in (1,5,20):
        count=Counter()
        for k,r in by_method[METHODS[0]].items():
            if r['coverage']!=str(coverage):continue
            new=by_method[METHODS[2]][k]
            if new['estimate']=='NA':count['No fit']+=1;continue
            delta=float(new['absolute_relative_error'])-float(r['absolute_relative_error'])
            count['Lower error' if delta < -1e-12 else 'Higher error' if delta>1e-12 else 'Equal']+=1
        outcomes.append(count)
    left=np.zeros(3)
    for category,color in [('Lower error','#187d91'),('Higher error','#bb4d51'),('Equal','#e0bb57'),('No fit','#c6cbd0')]:
        values=np.array([r[category] for r in outcomes]);totals=np.array([sum(r.values()) for r in outcomes])
        fraction=100*values/totals
        ax.barh(range(3),fraction,left=left,color=color,label=category)
        for j,(v,f) in enumerate(zip(values,fraction)):
            source('C','paired_outcome_count',j,int(v),coverage=(1,5,20)[j],status=category,n=int(totals[j]))
            if f>5:ax.text(left[j]+f/2,j,str(v),ha='center',va='center',fontsize=9,color='white' if category!='No fit' else '#30363b')
        left+=fraction
    ax.set(yticks=range(3),yticklabels=['1×','5×','20×'],xlim=(0,100),xlabel='All paired family conditions (%)',
           title='C  Multi-k versus native median: all outcomes')
    ax.legend(fontsize=8,ncols=2,loc='upper center',bbox_to_anchor=(.5,-.18))

    ax=axs[1,1];matrix=np.zeros((3,3))
    for i,coverage in enumerate((1,5,20)):
        for j,error in enumerate(ERRORS):
            chosen=[r for r in original if r['coverage']==str(coverage) and r['error_model']==error
                    and float(r['unit_substitution_rate'])==.02 and r['array_scope']=='factorial']
            value=100*sum(r['interval_contains_truth']=='True' for r in chosen)/len(chosen)
            matrix[i,j]=value;source('D','native_spread_truth_coverage_percent',j,value,coverage=coverage,error=error,n=len(chosen))
    im=ax.pcolormesh(np.arange(4)-.5,np.arange(4)-.5,matrix,vmin=0,vmax=100,cmap='cividis',rasterized=False)
    ax.invert_yaxis()
    for i in range(3):
        for j in range(3):ax.text(j,i,f'{matrix[i,j]:.1f}%',ha='center',va='center',color='white')
    ax.set(xticks=range(3),xticklabels=['Clean','Low error','High error'],yticks=range(3),yticklabels=['1×','5×','20×'],
           title='D  Native k-mer spread truth coverage',xlabel='2% unit divergence; spread is not a sampling CI')
    colorbar=fig.colorbar(im,ax=ax,fraction=.05,label='Observed coverage (%)')
    colorbar.solids.set_rasterized(False)

    ax=axs[2,0]
    for method,label,color in zip(METHODS,LABELS,COLORS):
        means=[]
        for coverage in (1,5,20):
            selected=[r for r in by_method[method].values() if r['coverage']==str(coverage)
                      and r['error_model']=='iid_high_indel' and r['array_scope']=='megabase']
            if any(r['estimate']=='NA' for r in selected):raise ValueError('Missing megabase estimate')
            values=[float(r['estimate'])/float(r['truth_copies']) for r in selected]
            means.append(statistics.mean(values))
            for r,v in zip(selected,values):
                ax.scatter(coverage,v,s=17,color=color,alpha=.5)
                source('E','estimate_to_truth_ratio',coverage,v,method=method,seed=r['seed'],family_id=r['family_id'])
        ax.plot((1,5,20),means,'o-',color=color,ms=4,label=label)
    ax.axhline(1,color='black',lw=.7,ls='--')
    ax.set(xscale='log',xticks=(1,5,20),xticklabels=['1×','5×','20×'],xlabel='Source coverage',ylabel='Estimated / planted copies',
           title='E  1.026-Mb arrays, high read errors')

    ax=axs[2,1]
    for divergence,color,label in [(0,'#687785','0%'),(.01,'#d68d32','1%: megabase'),(.02,'#187d91','2%')]:
        xs,ys=[],[]
        for k,r in by_method[METHODS[2]].items():
            if r['coverage']!='20' or float(r['unit_substitution_rate'])!=divergence or r['estimate']=='NA':continue
            x=float(fit_by_key[k]['max_absolute_log_residual']);y=100*abs(float(r['estimator_minus_sampling_oracle']))
            xs.append(x);ys.append(y)
            source('F','absolute_error_vs_source_oracle_percent',x,y,seed=r['seed'],condition_id=r['condition_id'],
                   family_id=r['family_id'],divergence=divergence)
        ax.scatter(xs,ys,s=12,color=color,alpha=.5,label=label)
    ax.set(xscale='symlog',xlabel='Maximum absolute log-fit residual',ylabel='Absolute estimator − oracle error (%)',
           title='F  Fit residual does not establish calibration')
    ax.set_xscale('symlog',linthresh=1e-4);ax.legend(title='Unit divergence',fontsize=8,title_fontsize=8)
    fig.suptitle('Conditional copy-number evaluation on three 10-Mb source genomes',fontsize=15)
    fig.supxlabel('Development data · 55 families/genome · 27 paired read conditions · points are not independent plants',fontsize=10)
    for extension in ('pdf','svg','png'):
        fig.savefig(outdir/f'factorial_quantification.{extension}',dpi=180)
    plt.close(fig)
    fields=list(dict.fromkeys(k for r in sources for k in r))
    write_table(outdir/'panel_source.tsv',({k:r.get(k) for k in fields} for r in sources),fields)
    input_paths=[replay/name for name in ('metrics.tsv','estimates.tsv','validation.json')]+[baseline/'copy_number_metrics.tsv']
    provenance=dict(script_sha256=digest_file(Path(__file__)),input_sha256={str(p):digest_file(p) for p in input_paths},
        output_sha256={p.name:digest_file(p) for p in outdir.iterdir() if p.is_file()},
        interpretation='development_diagnostic_not_final_paper_or_external_superiority')
    (outdir/'figure_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args();run(args.replay,args.baseline,args.outdir)
