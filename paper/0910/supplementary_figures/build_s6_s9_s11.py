#!/usr/bin/env python3
"""Render the final source-linked v5 supplementary figures S6, S9 and S11."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
EVID=ROOT/"paper"/"evidence"
OUT=Path(__file__).resolve().parent
COL={"tandemx":"#1679ab","srf_k151":"#e17c05","srf_k101":"#5b8f5a","competitive_mapping":"#7659a6"}
LAB={"tandemx":"TandemX","srf_k151":"SRF k=151","srf_k101":"SRF k=101","competitive_mapping":"competitive mapping*"}

def read(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter="\t"))
def num(x):return np.nan if x in {"", "NA", "N/A", None} else float(x)
def fig(n=2):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.linewidth":.7})
    f,a=plt.subplots(1,n,figsize=(15 if n==3 else 10.4,4.75));f.subplots_adjust(left=.06,right=.985,top=.84,bottom=.27,wspace=.36);return f,a
def caption(ax,letter,title,text):
    ax.set_title(f"{letter}  {title}",loc="left",fontsize=10,fontweight="bold",pad=8)
    ax.text(0,-.30,text,transform=ax.transAxes,ha="left",va="top",fontsize=7,color="#444")
def save(f,out,name):
    out.mkdir(parents=True,exist_ok=True)
    for ext in ("png","pdf","svg"):f.savefig(out/f"{name}.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
    plt.close(f)
def provenance(out, entries):
    with (out/"panel_source.tsv").open("w",newline="") as h:
        fields=["panel","source_file","selection","metric","evidence_boundary"]
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(entries)
def readme(out,title):
    (out/"README.md").write_text(f"# {title}\n\nRendered locally by `../build_s6_s9_s11.py` from the frozen sources in `panel_source.tsv`. PNG was visually reviewed; PDF is single-page; SVG preserves editable vector marks and text.\n")

def build_s6():
    source=EVID/"srf_formal_unified_v1"/"condition_summary.tsv"; data=read(source);out=OUT/"figure_s6_formal_srf_endpoint"
    methods=["tandemx","srf_k151","srf_k101","competitive_mapping"];conditions=[]
    for r in data:
        if r['condition'] not in conditions:conditions.append(r['condition'])
    lookup={(r['condition'],r['method']):r for r in data}; f,axes=fig(3);x=np.arange(len(conditions));w=.19
    for j,m in enumerate(methods):
        axes[0].bar(x+(j-1.5)*w,[num(lookup[(c,m)]['MARE_mean']) for c in conditions],w,color=COL[m],label=LAB[m])
    axes[0].set(xticks=x,xticklabels=conditions,ylabel="Mean absolute relative error",xlabel="Specified simulation condition")
    axes[0].tick_params(axis='x',rotation=35,labelsize=6);axes[0].grid(axis='y',color="#ddd",lw=.6);axes[0].legend(frameon=False,fontsize=6)
    caption(axes[0],"A","Abundance endpoint","Formal frozen record; each value is the mean across three simulation seeds.\n*Competitive mapping uses the TandemX catalogue and is not an independent discovery endpoint.")
    for j,m in enumerate(methods):
        axes[1].plot(x,[num(lookup[(c,m)]['base_union_recall_mean']) for c in conditions],"o-",color=COL[m],label=LAB[m])
    axes[1].set(xticks=x,xticklabels=conditions,ylim=(-.03,1.05),ylabel="Mean base-union recall",xlabel="Specified simulation condition")
    axes[1].tick_params(axis='x',rotation=35,labelsize=6);axes[1].grid(axis='y',color="#ddd",lw=.6);axes[1].legend(frameon=False,fontsize=6)
    caption(axes[1],"B","Assigned-base endpoint","All methods are shown at their formal specified endpoints; this panel does not rank runtime or memory.")
    for j,m in enumerate(methods):
        vals=[num(lookup[(c,m)]['family_recovery_mean']) for c in conditions]
        axes[2].scatter(x+(j-1.5)*.12,vals,s=25,color=COL[m],label=LAB[m])
    axes[2].set(xticks=x,xticklabels=conditions,ylim=(.7,1.05),ylabel="Mean family recovery",xlabel="Specified simulation condition")
    axes[2].tick_params(axis='x',rotation=35,labelsize=6);axes[2].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[2],"C","Catalogue endpoint","Competitive-mapping family recovery is unavailable because it reuses TandemX discovery; unavailable is not plotted as zero.")
    save(f,out,"figure_s6_formal_srf_endpoint")
    provenance(out,[
      {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"all formal conditions and methods","metric":"MARE_mean","evidence_boundary":"mean of three simulation seeds; mapping reuses TandemX catalogue"},
      {"panel":"B","source_file":str(source.relative_to(ROOT)),"selection":"all formal conditions and methods","metric":"base_union_recall_mean","evidence_boundary":"mean of three simulation seeds; no runtime/memory ranking"},
      {"panel":"C","source_file":str(source.relative_to(ROOT)),"selection":"all formal conditions and methods","metric":"family_recovery_mean","evidence_boundary":"competitive-mapping value unavailable by design and not zero"},
    ]);readme(out,"Supplementary Figure S6 — formal SRF endpoint")

def build_s9():
    pth=EVID/"macadamia_jansenii_donor_matched_collapse_v1"/"results"/"evaluation_primary_total_bases"/"family_metrics.tsv"
    sth=EVID/"macadamia_jansenii_donor_matched_collapse_v1"/"results"/"evaluation_reported_depth_sensitivity"/"family_metrics.tsv"
    ath=EVID/"macadamia_jansenii_donor_matched_collapse_v1"/"results"/"old_new_alignment_context"/"family_alignment_context.tsv"
    p=[r for r in read(pth) if r['eligibility']=='eligible'];s={r['family_id']:r for r in read(sth)};a={r['family_id']:r for r in read(ath)};out=OUT/"figure_s9_macadamia_reference_context";f,axes=fig(3)
    order=sorted(p,key=lambda r:num(r['new_read_ratio']));colors=["#c44e52" if r['reference_state']=='reference_collapse' else "#4c78a8" for r in order];x=np.arange(len(order))
    axes[0].scatter(x,[num(r['new_read_ratio']) for r in order],c=colors,s=28);axes[0].axhline(1,color="#999",lw=.8,ls="--")
    axes[0].set(yscale="log",xticks=x,xticklabels=[r['family_id'].replace('TXF','') for r in order],ylabel="New-assembly / read abundance ratio",xlabel="Eligible family")
    axes[0].tick_params(axis='x',rotation=65,labelsize=5);axes[0].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[0],"A","Primary reference comparison","Read-to-assembly ratios are abundance diagnostics.\nThe newer assembly is a reference proxy, not copy-number truth.")
    y1=[num(r['new_read_ratio']) for r in order];y2=[num(s[r['family_id']]['new_read_ratio']) for r in order]
    axes[1].scatter(y1,y2,c=colors,s=28);lo=min(y1+y2)*.8;hi=max(y1+y2)*1.2;axes[1].plot([lo,hi],[lo,hi],"--",color="#888",lw=.8)
    axes[1].set(xscale="log",yscale="log",xlim=(lo,hi),ylim=(lo,hi),xlabel="Primary ratio",ylabel="Reported-depth sensitivity ratio")
    axes[1].grid(color="#ddd",lw=.6);caption(axes[1],"B","Specified depth sensitivity","Same eligible families under two frozen read-depth inputs.\nRatios are not converted into physical sequence quantities.")
    aa=[r for r in order if r['family_id'] in a];x2=np.arange(len(aa))
    axes[2].scatter(x2,[num(a[r['family_id']]['any_alignment_fraction']) for r in aa],s=20,color="#4c78a8",label="any alignment")
    axes[2].scatter(x2,[num(a[r['family_id']]['primary_alignment_fraction']) for r in aa],s=20,color="#72b7b2",marker="s",label="primary alignment")
    axes[2].set(ylim=(-.03,1.05),xticks=x2,xticklabels=[r['family_id'].replace('TXF','') for r in aa],ylabel="New-assembly bases aligned",xlabel="Eligible family")
    axes[2].tick_params(axis='x',rotation=65,labelsize=5);axes[2].legend(frameon=False,fontsize=6);axes[2].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[2],"C","Old-to-new alignment context","Alignment supports reference comparison context.\nIt does not establish an absolute biological copy number.")
    save(f,out,"figure_s9_macadamia_reference_context")
    provenance(out,[
      {"panel":"A","source_file":str(pth.relative_to(ROOT)),"selection":"eligibility=eligible","metric":"new_read_ratio;reference_state","evidence_boundary":"new assembly is a donor-matched reference proxy"},
      {"panel":"B","source_file":f"{pth.relative_to(ROOT)};{sth.relative_to(ROOT)}","selection":"same eligible family IDs","metric":"new_read_ratio","evidence_boundary":"read-depth sensitivity only"},
      {"panel":"C","source_file":str(ath.relative_to(ROOT)),"selection":"eligible family IDs with alignment context","metric":"any_alignment_fraction;primary_alignment_fraction","evidence_boundary":"alignment context, not copy-number truth"},
    ]);readme(out,"Supplementary Figure S9 — Macadamia donor-matched reference context")

def build_s11():
    qc=ROOT/"paper"/"0910"/"source_data"/"figure6_txf000695_qc.tsv";orth=EVID/"orthogonal_abundance_validation_v1"/"figures_v2"/"panel_source.tsv";q=read(qc);o=read(orth);out=OUT/"figure_s11_txf000695_qc"
    f,axes=fig(3);ks=[r['k'] for r in q];x=np.arange(2)
    axes[0].bar(x-.16,[num(r['observed_fraction']) for r in q],.30,color="#4c78a8",label="observed ≥2")
    axes[0].bar(x+.16,[1-num(r['observed_fraction']) for r in q],.30,color="#c9c9c9",label="below 2")
    axes[0].set(ylim=(0,1),xticks=x,xticklabels=[f"k={k}" for k in ks],ylabel="Diagnostic k-mer fraction",xlabel="Frozen k-mer setting")
    axes[0].legend(frameon=False,fontsize=7);axes[0].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[0],"A","Diagnostic-k-mer observation","Observed fraction uses the frozen ≥2-count threshold.\nThis QC panel does not change the estimator or family eligibility.")
    axes[1].bar(x-.16,[num(r['median_depth']) for r in q],.30,color="#e17c05",label="median depth")
    axes[1].bar(x+.16,[num(r['depth_mad']) for r in q],.30,color="#72b7b2",label="depth MAD")
    axes[1].set(yscale="log",xticks=x,xticklabels=[f"k={k}" for k in ks],ylabel="Diagnostic-k-mer depth",xlabel="Frozen k-mer setting")
    axes[1].legend(frameon=False,fontsize=7);axes[1].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[1],"B","Median-boundary sensitivity","The two frozen k settings occupy different depth modes.\nThey are shown separately and are not combined into one estimate.")
    candidate=[r for r in o if r['panel']=='A' and r['species_key']=='macadamia' and r['family_id']=='TXF000695']
    labels=[];vals=[]
    for m in ['Newer assembly','Frozen HiFi','Illumina k=21','Illumina k=31']:
        r=next(r for r in candidate if r['platform_or_measure']==m);labels.append(m);vals.append(num(r['value_bp']))
    axes[2].bar(range(4),vals,color=["#4c78a8","#e17c05","#59a14f","#c44e52"])
    axes[2].set(yscale="log",xticks=range(4),xticklabels=["assembly","HiFi","Illumina\nk=21","Illumina\nk=31"],ylabel="Abundance estimate (bp; log scale)",xlabel="Evidence source")
    axes[2].grid(axis='y',color="#ddd",lw=.6)
    caption(axes[2],"C","Unresolved candidate context","The cross-k discrepancy leaves TXF000695 unresolved.\nThe final interpretation remains unchanged.")
    save(f,out,"figure_s11_txf000695_qc")
    provenance(out,[
      {"panel":"A","source_file":str(qc.relative_to(ROOT)),"selection":"TXF000695; k=21 and k=31","metric":"observed_fraction","evidence_boundary":"frozen threshold QC; no estimator change"},
      {"panel":"B","source_file":str(qc.relative_to(ROOT)),"selection":"TXF000695; k=21 and k=31","metric":"median_depth;depth_mad","evidence_boundary":"settings shown separately; not pooled"},
      {"panel":"C","source_file":str(orth.relative_to(ROOT)),"selection":"panel=A; macadamia TXF000695","metric":"value_bp;final_interpretation","evidence_boundary":"unresolved evidence context; no new Result"},
    ]);readme(out,"Supplementary Figure S11 — TXF000695 diagnostic QC")

if __name__=="__main__":build_s6();build_s9();build_s11()
