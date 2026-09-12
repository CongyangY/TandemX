#!/usr/bin/env python3
"""Render source-linked Supplementary Figures S2--S5 for manuscript v5.

Each panel is calculated directly from the frozen evidence TSV files listed in
the figure-local panel_source.tsv files.  These figures deliberately omit
runtime/RSS rankings, internal TandemX variants and retired exploratory replay.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[3]
EVID = ROOT / "paper" / "evidence"
OUT = Path(__file__).resolve().parent

COLORS = {"tandemx": "#1679ab", "tidehunter": "#dd8452", "trf": "#5b8f5a", "tidecluster": "#875ca5"}
LABELS = {"tandemx": "TandemX", "tidehunter": "TideHunter", "trf": "TRF", "tidecluster": "TideCluster"}


def rows(path: Path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def f(row, key):
    value = row.get(key, "")
    return np.nan if value in {"", "NA", "N/A"} else float(value)


def mean_or_nan(values):
    values = [x for x in values if not np.isnan(x)]
    return np.mean(values) if values else np.nan


def start(nrows=1, ncols=2):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7})
    fig, ax = plt.subplots(nrows, ncols, figsize=(10.1, 4.35 if nrows == 1 else 7.25), constrained_layout=True)
    return fig, np.ravel(ax)


def finish(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def panel_table(outdir, entries):
    with (outdir / "panel_source.tsv").open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=["panel", "source_file", "selection", "metric", "evidence_boundary"])
        writer.writeheader(); writer.writerows(entries)


def note(ax, letter, title, text):
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=10, pad=8)
    ax.text(0, -0.25, text, transform=ax.transAxes, ha="left", va="top", fontsize=7, color="#444444", wrap=True)


def build_s2():
    source = EVID / "cascade_gap_free_validation_v1" / "run" / "summary.tsv"
    data = rows(source)
    outdir = OUT / "figure_s2_production_validation"
    positive = [x for x in data if x["scenario"] not in {"at_rich_control", "low_complexity_control", "dispersed_control"}]
    scenarios = sorted({x["scenario"] for x in positive})
    fig, ax = start()
    for j, tool in enumerate(["tandemx", "tidehunter"]):
        vals = []
        for s in scenarios:
            r = next(x for x in positive if x["scenario"] == s and x["tool"] == tool)
            vals.append(f(r, "array_f1"))
        ax[0].scatter(np.arange(len(scenarios)) + (j-.5)*.20, vals, s=27, color=COLORS[tool], label=LABELS[tool], zorder=3)
    ax[0].set(ylim=(-.04, 1.08), ylabel="Array F1", xlabel="Specified simulated scenario")
    ax[0].set_xticks(range(len(scenarios))); ax[0].set_xticklabels([x.replace("_", "\n") for x in scenarios], rotation=0, fontsize=6)
    ax[0].grid(axis="y", color="#dddddd", lw=.6); ax[0].legend(frameon=False, ncol=2, loc="lower left")
    note(ax[0], "A", "Positive-array endpoint", "One frozen validation seed per condition; points are scenario outputs, not biological replicates.")
    controls = ["at_rich_control", "low_complexity_control", "dispersed_control"]
    for j, tool in enumerate(["tandemx", "tidehunter"]):
        vals = [f(next(x for x in data if x["scenario"] == s and x["tool"] == tool), "negative_read_call_rate") for s in controls]
        ax[1].bar(np.arange(3)+(j-.5)*.32, vals, width=.30, color=COLORS[tool], label=LABELS[tool])
    ax[1].set(ylabel="Negative-read call rate", xlabel="Negative control", ylim=(0, .55))
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels([x.replace("_control", "").replace("_", "\n") for x in controls])
    ax[1].grid(axis="y", color="#dddddd", lw=.6); ax[1].legend(frameon=False, ncol=2)
    note(ax[1], "B", "Specified negative controls", "Zero is an observed value; no external-tool runtime or memory ranking is shown.")
    finish(fig, outdir, "figure_s2_production_validation")
    panel_table(outdir, [
        {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"13 planted positive scenarios; TandemX and TideHunter","metric":"array_f1","evidence_boundary":"single frozen validation seed per scenario"},
        {"panel":"B","source_file":str(source.relative_to(ROOT)),"selection":"three specified negative controls; TandemX and TideHunter","metric":"negative_read_call_rate","evidence_boundary":"single frozen validation seed per control"},
    ])
    return outdir


def build_s3():
    source = EVID / "srf_formal_unified_v1" / "condition_summary.tsv"
    data = [x for x in rows(source) if x["method"] == "tandemx"]
    outdir = OUT / "figure_s3_production_quantification"
    conditions = [x["condition"] for x in data]
    fig, ax = start()
    mare = [f(x, "MARE_mean") for x in data]
    ax[0].bar(range(len(data)), mare, color="#1679ab")
    ax[0].set(xticks=range(len(data)), xticklabels=conditions, ylabel="Mean absolute relative error", xlabel="Specified simulation condition")
    ax[0].tick_params(axis="x", rotation=25); ax[0].grid(axis="y", color="#dddddd", lw=.6)
    note(ax[0], "A", "Read abundance endpoint", "TandemX production endpoint only; each bar summarizes three simulation seeds in the formal unified record.")
    rec = [f(x, "base_union_recall_mean") for x in data]
    prec = [f(x, "base_union_precision_mean") for x in data]
    xs = np.arange(len(data)); ax[1].plot(xs, rec, "o-", color="#1679ab", label="base-union recall")
    ax[1].plot(xs, prec, "s-", color="#425466", label="base-union precision")
    ax[1].set(xticks=xs, xticklabels=conditions, ylim=(.72, 1.01), ylabel="Mean endpoint value", xlabel="Specified simulation condition")
    ax[1].tick_params(axis="x", rotation=25); ax[1].grid(axis="y", color="#dddddd", lw=.6); ax[1].legend(frameon=False, fontsize=7)
    note(ax[1], "B", "Assigned-base endpoint", "This is simulation-based endpoint behavior and is not a population-level plant validation.")
    finish(fig, outdir, "figure_s3_production_quantification")
    panel_table(outdir, [
        {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"method=tandemx; all formal conditions","metric":"MARE_mean","evidence_boundary":"mean across three simulation seeds"},
        {"panel":"B","source_file":str(source.relative_to(ROOT)),"selection":"method=tandemx; all formal conditions","metric":"base_union_recall_mean;base_union_precision_mean","evidence_boundary":"mean across three simulation seeds"},
    ])
    return outdir


def build_s4():
    source = EVID / "tidecluster_factorial_validation_v3" / "results" / "summary.tsv"
    fates = EVID / "tidecluster_factorial_validation_v3" / "results" / "cell_fates.tsv"
    data = rows(source); outdir = OUT / "figure_s4_tidecluster_endpoint"
    good = [x for x in data if x["status"] == "ok"]
    labels = [f"{x['seed']}\n{x['setting'].replace('_sensitivity','')}" for x in good]
    fig, ax = start()
    xs = np.arange(len(good)); ax[0].bar(xs-.16, [f(x,"array_recall") for x in good], .30, color="#875ca5", label="array recall")
    ax[0].bar(xs+.16, [f(x,"array_precision") for x in good], .30, color="#4e6a88", label="array precision")
    ax[0].set(xticks=xs, xticklabels=labels, ylim=(.82,1.01), ylabel="Endpoint value", xlabel="Completed external-tool cells")
    ax[0].tick_params(axis="x", labelsize=6); ax[0].grid(axis="y", color="#dddddd", lw=.6); ax[0].legend(frameon=False, fontsize=7)
    note(ax[0], "A", "TideCluster completed cells", "Only cells with evaluable output are plotted; technical seeds do not represent biological replicates.")
    fate_rows = rows(fates); sets = ["default_primary", "matched_period_sensitivity"]; seeds=["6401","6402","6403"]
    status = {(x['seed'],x['setting']): x['status'] for x in fate_rows}
    matrix=np.array([[1 if status.get((seed,setting)) == 'ok' else 0 for setting in sets] for seed in seeds])
    for i in range(3):
        for j in range(2):
            ax[1].add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor="#875ca5" if matrix[i, j] else "#e6e6e6", edgecolor="white", linewidth=1))
    for i in range(3):
        for j in range(2): ax[1].text(j,i,"completed" if matrix[i,j] else "unavailable",ha="center",va="center",fontsize=7)
    ax[1].set(xlim=(-.5, 1.5), ylim=(2.5, -.5), xticks=[0,1],xticklabels=["default", "matched period"],yticks=[0,1,2],yticklabels=seeds,xlabel="Frozen setting",ylabel="Technical seed")
    note(ax[1], "B", "External-tool completion record", "Unavailable accuracy is not plotted as zero; failures are retained as a completion record, not a performance ranking.")
    finish(fig, outdir, "figure_s4_tidecluster_endpoint")
    panel_table(outdir, [
        {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"status=ok TideCluster cells","metric":"array_recall;array_precision","evidence_boundary":"technical seeds on planted factorial assembly"},
        {"panel":"B","source_file":str(fates.relative_to(ROOT)),"selection":"all seed-by-setting cells","metric":"status","evidence_boundary":"completion status only; unavailable accuracy is not zero"},
    ])
    return outdir


def build_s5():
    source = EVID / "cascade_native_screen_heldout_v1" / "run" / "summary.tsv"
    data = rows(source); outdir = OUT / "figure_s5_heldout_external_endpoint"
    pos = [x for x in data if x["scenario"] not in {"at_rich_control", "low_complexity_control", "dispersed_control"}]
    scen = sorted({x["scenario"] for x in pos}); methods=["tandemx", "tidehunter", "trf"]
    fig, ax = start()
    for j, tool in enumerate(methods):
        vals=[]
        for s in scen:
            rr=[x for x in pos if x["scenario"]==s and x["tool"]==tool]
            vals.append(mean_or_nan([f(x,"array_f1") for x in rr]))
        ax[0].scatter(np.arange(len(scen))+(j-1)*.18,vals,color=COLORS[tool],s=18,label=LABELS[tool],zorder=3)
    ax[0].set(ylim=(-.05,1.07),xticks=range(len(scen)),xticklabels=[s.replace('_','\n') for s in scen],ylabel="Mean array F1",xlabel="Held-out simulated scenario")
    ax[0].tick_params(axis="x",labelsize=6);ax[0].grid(axis="y",color="#dddddd",lw=.6);ax[0].legend(frameon=False,ncol=3,fontsize=7)
    note(ax[0], "A", "Held-out array endpoint", "Means summarize the three frozen held-out seeds per scenario; blank values remain unavailable.")
    controls=["at_rich_control","low_complexity_control","dispersed_control"]
    for j,tool in enumerate(methods):
        vals=[]
        for s in controls:
            rr=[x for x in data if x['scenario']==s and x['tool']==tool]
            vals.append(mean_or_nan([f(x,'negative_read_call_rate') for x in rr]))
        ax[1].bar(np.arange(3)+(j-1)*.24, vals,.22,color=COLORS[tool],label=LABELS[tool])
    ax[1].set(ylim=(0,.75),xticks=range(3),xticklabels=[s.replace('_control','').replace('_','\n') for s in controls],ylabel="Mean negative-read call rate",xlabel="Held-out negative control")
    ax[1].grid(axis="y",color="#dddddd",lw=.6);ax[1].legend(frameon=False,ncol=3,fontsize=7)
    note(ax[1], "B", "Held-out negative controls", "No runtime, memory, or internal-version comparison is presented in this supplementary endpoint record.")
    finish(fig, outdir, "figure_s5_heldout_external_endpoint")
    panel_table(outdir, [
        {"panel":"A","source_file":str(source.relative_to(ROOT)),"selection":"13 positive scenarios; TandemX, TideHunter and TRF; heldout seeds","metric":"array_f1","evidence_boundary":"mean across three held-out simulation seeds"},
        {"panel":"B","source_file":str(source.relative_to(ROOT)),"selection":"three negative controls; TandemX, TideHunter and TRF; heldout seeds","metric":"negative_read_call_rate","evidence_boundary":"mean across three held-out simulation seeds"},
    ])
    return outdir


def write_readme(outdir, title):
    (outdir / "README.md").write_text(f"# {title}\n\nRendered by `../build_s2_s5.py` from the frozen evidence files enumerated in `panel_source.tsv`. PNG is the visual-review asset; PDF is one page; SVG retains editable text and paths.\n")


if __name__ == "__main__":
    outputs = [(build_s2(), "Supplementary Figure S2 — specified production validation"), (build_s3(), "Supplementary Figure S3 — production quantification endpoint"), (build_s4(), "Supplementary Figure S4 — TideCluster endpoint and completion"), (build_s5(), "Supplementary Figure S5 — held-out external endpoint")]
    for folder, title in outputs: write_readme(folder, title)
