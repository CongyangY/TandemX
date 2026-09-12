#!/usr/bin/env python3
"""Render the three-panel Figure 2 v3.1 review preview from frozen TSVs."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).parent
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none"})
COL = {
    "tandemx": "#357bbd", "trf": "#df6b44", "tidehunter": "#3e9b6e",
    "srf_k151": "#8b6cc7", "srf_k101": "#c25cb4", "competitive_mapping": "#596a7c",
}
NAME = {
    "tandemx": "TandemX", "trf": "TRF", "tidehunter": "TideHunter",
    "srf_k151": "SRF k=151", "srf_k101": "SRF k=101",
    "competitive_mapping": "competitive mapping",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def panel_title(ax: plt.Axes, letter: str, title: str, subtitle: str) -> None:
    ax.text(-0.13, 1.10, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="bottom")
    ax.text(-0.02, 1.10, title, transform=ax.transAxes, fontsize=9.2,
            fontweight="bold", va="bottom")
    ax.text(-0.02, 1.035, subtitle, transform=ax.transAxes, fontsize=7,
            color="#607083", va="bottom")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d7dde5", linewidth=.55)
    ax.set_axisbelow(True)


def panel_a(ax: plt.Axes) -> None:
    rows = read_tsv(ROOT / "paper/evidence/cascade_native_screen_heldout_v1/run/summary.tsv")
    tools = ["tandemx", "trf", "tidehunter"]
    metrics = [
        ("sequence_family_recall", "family\nrecovery"),
        ("array_f1", "interval\nF1"),
        ("base_union_f1", "base-union\nF1"),
    ]
    scores: dict[str, list[float]] = {tool: [] for tool in tools}
    for field, _ in metrics:
        for tool in tools:
            values = [float(row[field]) for row in rows if row["tool"] == tool and row[field] not in ("", "NA")]
            scores[tool].append(float(np.mean(values)))
    x = np.arange(len(metrics)); width = .23
    for j, tool in enumerate(tools):
        bars = ax.bar(x + (j - 1) * width, scores[tool], width, color=COL[tool], label=NAME[tool])
        for bar, value in zip(bars, scores[tool]):
            if value < .9995:
                ax.text(bar.get_x() + bar.get_width() / 2, value + .016, f"{value:.3f}",
                        ha="center", va="bottom", fontsize=7)
    panel_title(ax, "A", "Read discovery and interval recovery", "13 held-out scenarios; 3 seeds each")
    ax.set_ylim(0, 1.08); ax.set_yticks([0, .25, .5, .75, 1.0]); ax.set_ylabel("mean score")
    ax.set_xticks(x, [label for _, label in metrics]); clean_axis(ax)
    ax.legend(frameon=False, ncol=3, fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.19))


def panel_b(ax: plt.Axes) -> None:
    rows = read_tsv(ROOT / "paper/evidence/srf_formal_unified_v1/condition_summary.tsv")
    by: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by[row["method"]][row["condition"]] = row
    discovery = ax.inset_axes([.00, .13, .32, .78])
    methods = ["tandemx", "srf_k151", "srf_k101"]
    recovered = [round(sum(float(by[m][condition]["family_recovery_mean"]) for condition in by[m]) * 9) for m in methods]
    bars = discovery.bar(range(3), recovered, color=[COL[m] for m in methods], width=.72)
    discovery.set_ylim(0, 58); discovery.set_yticks([0, 27, 54]); discovery.set_ylabel("recovered conditions", fontsize=7)
    discovery.set_xticks(range(3), ["Tandem\nX", "SRF\n151", "SRF\n101"], fontsize=7); discovery.set_title("Discovery", fontsize=8, pad=3)
    discovery.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, recovered):
        discovery.text(bar.get_x() + bar.get_width() / 2, value + 1, f"{value}/54", ha="center", fontsize=7)
    mare = ax.inset_axes([.52, .13, .48, .78])
    conditions = ["clean", "substitution", "indel", "divergence", "low_abundance", "shared_fragment"]
    for method in ["tandemx", "srf_k151", "srf_k101", "competitive_mapping"]:
        mare.plot(range(6), [float(by[method][condition]["MARE_mean"]) for condition in conditions],
                  marker="o", ms=3.4, lw=1.45, color=COL[method],
                  label="mapping (discovery N/A)" if method == "competitive_mapping" else NAME[method])
    mare.set_ylim(-.04, 1.08); mare.set_yticks([0, .5, 1]); mare.set_ylabel("positive-family MARE", fontsize=7)
    mare.set_xticks(range(6), ["clean", "subs.", "indel", "div.", "low", "shared"], rotation=28, ha="right", fontsize=7)
    mare.set_title("Abundance error", fontsize=8, pad=3); clean_axis(mare)
    mare.legend(frameon=False, fontsize=7, loc="upper left", handlelength=1.4)
    panel_title(ax, "B", "Discovery recovery and abundance error", "18 datasets; 3 seeds per condition")
    ax.set_axis_off()


def panel_c(ax: plt.Axes) -> None:
    rows = read_tsv(ROOT / "paper/evidence/srf_formal_unified_v1/condition_summary.tsv")
    shared = {row["method"]: row for row in rows if row["condition"] == "shared_fragment"}
    methods = ["tandemx", "srf_k151", "srf_k101", "competitive_mapping"]
    values = [float(shared[method]["negative_read_predicted_bp_mean"]) for method in methods]
    bars = ax.bar(range(4), values, color=[COL[m] for m in methods], width=.68)
    panel_title(ax, "C", "False attribution", "Shared-fragment background (n=3)")
    ax.set_ylabel("negative-read predicted bp")
    ax.set_xticks(range(4), ["TandemX", "SRF151", "SRF101", "mapping"], fontsize=7)
    ax.set_ylim(0, 36000); ax.set_yticks([0, 10000, 20000, 30000])
    clean_axis(ax)
    for bar, value in zip(bars, values):
        y = value + 1000 if value >= 1000 else value + 650
        ax.text(bar.get_x() + bar.get_width() / 2, y, f"{value:,.0f}", ha="center", va="bottom", fontsize=7,
                rotation=90 if value >= 10000 else 0)


def main() -> None:
    fig = plt.figure(figsize=(10.6, 4.35))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.35, 1.05], wspace=.60)
    panel_a(fig.add_subplot(grid[0, 0]))
    panel_b(fig.add_subplot(grid[0, 1]))
    panel_c(fig.add_subplot(grid[0, 2]))
    fig.subplots_adjust(left=.07, right=.985, top=.80, bottom=.22)
    fig.savefig(OUT / "figure2_v3_1_production_only_preview.png", dpi=300)


if __name__ == "__main__":
    main()
