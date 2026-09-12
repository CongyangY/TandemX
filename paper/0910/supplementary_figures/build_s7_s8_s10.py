#!/usr/bin/env python3
"""Render evidence-bounded Supplementary Figures S7, S8 and S10 for v5."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[3]
EVID = ROOT / "paper" / "evidence"
OUT = Path(__file__).resolve().parent


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def num(x):
    return np.nan if x in {"", "NA", "N/A", None} else float(x)


def figure(ncols=2):
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":8, "axes.linewidth":.7})
    fig, axes = plt.subplots(1, ncols, figsize=(14.7 if ncols == 3 else 10.4, 4.75))
    fig.subplots_adjust(left=.06, right=.985, top=.84, bottom=.27, wspace=.37)
    return fig, axes


def save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png","pdf","svg"):
        fig.savefig(out / f"{name}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def caption(ax, letter, title, text):
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax.text(0,-.30,text,transform=ax.transAxes,ha="left",va="top",color="#444",fontsize=7,wrap=True)


def provenance(out, data):
    with (out/"panel_source.tsv").open("w",newline="") as h:
        fields=["panel","source_file","selection","metric","evidence_boundary"]
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(data)


def readme(out, title, script):
    (out/"README.md").write_text(f"# {title}\n\nRendered locally by `{script}` from the frozen sources listed in `panel_source.tsv`. PNG was visually reviewed; PDF is one page; SVG uses editable vector marks and text.\n")


def build_s7():
    qc = EVID/"multispecies_input_qc"/"figures_v2"/"panel_source.tsv"
    diag = EVID/"multispecies_real_diagnostics"/"figures_v2"/"panel_source.tsv"
    q, d = read(qc), read(diag); out=OUT/"figure_s7_multispecies_inputs"
    fig, axes=figure(3)
    a=[r for r in q if r['metric']=='total_bases_Gb']
    axes[0].bar(range(len(a)),[num(r['value']) for r in a],color="#377eb8")
    axes[0].set(xticks=range(len(a)),xticklabels=[r['reported_material'] for r in a],ylabel="Full-library bases (Gb)",xlabel="Reported material")
    axes[0].tick_params(axis="x",rotation=50,labelsize=6);axes[0].grid(axis="y",color="#ddd",lw=.6)
    caption(axes[0],"A","Input-library scale","Ten complete libraries from eight reported species; this panel reports input only.")
    by={}
    for r in q:
        if r['metric'] in {'median_read_length_kb','read_n50_kb'}: by.setdefault(r['reported_material'],{})[r['metric']]=num(r['value'])
    labs=list(by); xs=np.arange(len(labs)); axes[1].scatter(xs-.1,[by[x]['median_read_length_kb'] for x in labs],color="#377eb8",label="median")
    axes[1].scatter(xs+.1,[by[x]['read_n50_kb'] for x in labs],color="#7a5195",marker="s",label="N50")
    axes[1].set(xticks=xs,xticklabels=labs,ylabel="Read length (kb)",xlabel="Reported material")
    axes[1].tick_params(axis="x",rotation=50,labelsize=6); axes[1].legend(frameon=False,fontsize=7);axes[1].grid(axis="y",color="#ddd",lw=.6)
    caption(axes[1],"B","Input-read length","Per-library descriptive QC; materials are not biological replicates or a cross-species accuracy test.")
    t=[r for r in d if r['tool']=='tandemx']
    samples=sorted({r['material'] for r in t})
    for mat in samples:
        rr=sorted([r for r in t if r['material']==mat],key=lambda x:num(x['input_Gb']))
        axes[2].plot([num(r['input_Gb']) for r in rr],[num(r['observed_union_base_fraction']) for r in rr],"o-",label=mat)
    axes[2].set(xlabel="Random read subset (Gb)",ylabel="Observed called-base fraction")
    axes[2].grid(color="#ddd",lw=.6);axes[2].legend(frameon=False,fontsize=6,ncol=2)
    caption(axes[2],"C","Observed TandemX calls","Descriptive calls on matched read subsets; no external-tool resource ranking is included.")
    save(fig,out,"figure_s7_multispecies_inputs")
    provenance(out,[
      {"panel":"A","source_file":str(qc.relative_to(ROOT)),"selection":"metric=total_bases_Gb; all full libraries","metric":"total_bases_Gb","evidence_boundary":"input QC only"},
      {"panel":"B","source_file":str(qc.relative_to(ROOT)),"selection":"metric=median_read_length_kb or read_n50_kb","metric":"read-length summaries","evidence_boundary":"input QC only"},
      {"panel":"C","source_file":str(diag.relative_to(ROOT)),"selection":"tool=tandemx; matched random subsets","metric":"input_Gb;observed_union_base_fraction","evidence_boundary":"descriptive output; not biological replication"},
    ]); readme(out,"Supplementary Figure S7 — multispecies input and output context","../build_s7_s8_s10.py")


def build_s8():
    primary=EVID/"ey15_donor_matched_collapse_v1"/"results"/"evaluation_primary_total_bases"/"family_metrics.tsv"
    sens=EVID/"ey15_donor_matched_collapse_v1"/"results"/"evaluation_depth107_sensitivity"/"family_metrics.tsv"
    align=EVID/"ey15_donor_matched_collapse_v1"/"results"/"old_new_alignment_context"/"family_alignment_context.tsv"
    p=[r for r in read(primary) if r['eligibility']=='eligible']; s={r['family_id']:r for r in read(sens)}; al={r['family_id']:r for r in read(align)}
    out=OUT/"figure_s8_ey15_reference_context"; fig,axes=figure(3); order=sorted(p,key=lambda r:num(r['new_read_ratio']))
    colors=["#c44e52" if r['reference_state']=='reference_collapse' else "#4c78a8" for r in order]
    xs=np.arange(len(order)); axes[0].scatter(xs,[num(r['new_read_ratio']) for r in order],c=colors,s=33)
    axes[0].axhline(1,color="#999",lw=.8,ls="--"); axes[0].set(yscale="log",xticks=xs,xticklabels=[r['family_id'].replace('TXF','') for r in order],ylabel="New-assembly / read abundance ratio",xlabel="Eligible family")
    axes[0].tick_params(axis='x',rotation=65,labelsize=6);axes[0].grid(axis="y",color="#ddd",lw=.6)
    caption(axes[0],"A","Primary reference comparison","Read-to-assembly ratios are abundance diagnostics.\nThe newer assembly is a reference proxy, not copy-number truth.")
    y1=[num(r['new_read_ratio']) for r in order]; y2=[num(s[r['family_id']]['new_read_ratio']) for r in order]
    axes[1].scatter(y1,y2,c=colors,s=33); lo=min(y1+y2)*.8;hi=max(y1+y2)*1.2;axes[1].plot([lo,hi],[lo,hi],"--",color="#888",lw=.8)
    axes[1].set(xscale="log",yscale="log",xlim=(lo,hi),ylim=(lo,hi),xlabel="Primary ratio",ylabel="Depth-107 sensitivity ratio")
    axes[1].grid(color="#ddd",lw=.6);caption(axes[1],"B","Specified depth sensitivity","Same eligible families under the two frozen read-depth inputs.\nRatios are not converted into physical sequence quantities.")
    aa=[r for r in order if r['family_id'] in al]; x=np.arange(len(aa)); axes[2].bar(x-.17,[num(al[r['family_id']]['primary_alignment_fraction']) for r in aa],.32,color="#4c78a8",label="primary")
    axes[2].bar(x+.17,[num(al[r['family_id']]['same_chromosome_primary_mapq20_alignment_fraction']) for r in aa],.32,color="#72b7b2",label="same chr., MAPQ≥20")
    axes[2].set(ylim=(0,1.05),xticks=x,xticklabels=[r['family_id'].replace('TXF','') for r in aa],ylabel="New-assembly bases aligned",xlabel="Eligible family")
    axes[2].tick_params(axis='x',rotation=65,labelsize=6);axes[2].grid(axis="y",color="#ddd",lw=.6);axes[2].legend(frameon=False,fontsize=6)
    caption(axes[2],"C","Old-to-new alignment context","Alignment supports reference comparison context.\nIt does not establish an absolute biological copy number.")
    save(fig,out,"figure_s8_ey15_reference_context")
    provenance(out,[
      {"panel":"A","source_file":str(primary.relative_to(ROOT)),"selection":"eligibility=eligible","metric":"new_read_ratio;reference_state","evidence_boundary":"new assembly is a donor-matched reference proxy"},
      {"panel":"B","source_file":f"{primary.relative_to(ROOT)};{sens.relative_to(ROOT)}","selection":"same eligible family IDs","metric":"new_read_ratio","evidence_boundary":"read-depth sensitivity only"},
      {"panel":"C","source_file":str(align.relative_to(ROOT)),"selection":"eligible family IDs with alignment context","metric":"primary_alignment_fraction;same_chromosome_primary_mapq20_alignment_fraction","evidence_boundary":"alignment context, not copy-number truth"},
    ]);readme(out,"Supplementary Figure S8 — Ey15 donor-matched reference context","../build_s7_s8_s10.py")


def build_s10():
    source=EVID/"orthogonal_abundance_validation_v1"/"figures_v2"/"panel_source.tsv"; data=read(source);out=OUT/"figure_s10_orthogonal_candidates"
    vals=[r for r in data if r['panel']=='A']; families=[]
    for r in vals:
        key=(r['species_key'],r['family_id'])
        if key not in families:families.append(key)
    methods=["Newer assembly","Frozen HiFi","Illumina k=21","Illumina k=31"]
    fig,axes=figure(2); xs=np.arange(len(families)); width=.18
    for j,m in enumerate(methods):
        y=[]
        for sp,fam in families:
            rr=next((r for r in vals if r['species_key']==sp and r['family_id']==fam and r['platform_or_measure']==m),None)
            y.append(num(rr['value_bp']) if rr else np.nan)
        axes[0].bar(xs+(j-1.5)*width,y,width,label=m)
    axes[0].set(yscale="log",xticks=xs,xticklabels=[f"{sp[:3]}\n{fam.replace('TXF','')}" for sp,fam in families],ylabel="Abundance estimate (bp; log scale)",xlabel="Candidate family")
    axes[0].tick_params(axis="x",labelsize=6);axes[0].legend(frameon=False,fontsize=6,ncol=2);axes[0].grid(axis="y",color="#ddd",lw=.6)
    caption(axes[0],"A","Candidate abundance evidence","Methods estimate abundance, not physical copy number.\nONT directional ranges remain in the source record.")
    states={}
    for r in vals:
        if r['platform_or_measure']=='Frozen HiFi':states[(r['species_key'],r['family_id'])]=r['final_interpretation']
    for i,key in enumerate(families):
        state=states.get(key,'unresolved'); color="#c44e52" if state=='orthogonal_supports_residual_collapse' else "#bdbdbd"
        axes[1].add_patch(Rectangle((-.45,i-.42),.9,.84,facecolor=color,edgecolor="white")); axes[1].text(0,i,"supported" if color=="#c44e52" else "unresolved",ha="center",va="center",fontsize=8)
    axes[1].set(xlim=(-.5,.5),ylim=(len(families)-.5,-.5),xticks=[],yticks=range(len(families)),yticklabels=[f"{sp}: {fam}" for sp,fam in families],ylabel="Candidate family")
    caption(axes[1],"B","Evidence state after orthogonal review","Three candidates have direction support; three remain unresolved.\nThe evidence states do not measure physical sequence quantities.")
    save(fig,out,"figure_s10_orthogonal_candidates")
    provenance(out,[
      {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"panel=A; six candidate families; assembly/HiFi/Illumina measures","metric":"value_bp","evidence_boundary":"abundance estimates; ONT ranges directional only"},
      {"panel":"B","source_file":str(source.relative_to(ROOT)),"selection":"panel=A; Frozen HiFi evidence_state","metric":"evidence_state","evidence_boundary":"residual-collapse direction or unresolved; no physical missing-base inference"},
    ]);readme(out,"Supplementary Figure S10 — orthogonal abundance candidates","../build_s7_s8_s10.py")


if __name__=="__main__":
    build_s7();build_s8();build_s10()
