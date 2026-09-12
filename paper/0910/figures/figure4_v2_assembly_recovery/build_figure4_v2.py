#!/usr/bin/env python3
"""Render Figure 4 from frozen Ey15-2 and Macadamia family metrics."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).parent
EY_METRICS = ROOT / "paper/evidence/ey15_donor_matched_collapse_v1/results/evaluation_primary_total_bases/family_metrics.tsv"
MAC_METRICS = ROOT / "paper/evidence/macadamia_jansenii_donor_matched_collapse_v1/results/evaluation_primary_total_bases/family_metrics.tsv"
EY_SUMMARY = ROOT / "paper/evidence/ey15_donor_matched_collapse_v1/results/evaluation_primary_total_bases/summary.json"
MAC_SUMMARY = ROOT / "paper/evidence/macadamia_jansenii_donor_matched_collapse_v1/results/evaluation_primary_total_bases/summary.json"

BLUE = "#2D76B7"
CHARCOAL = "#3B4653"
GREY = "#6B7885"
LIGHT = "#DCE4EC"
ORANGE = "#D57A1F"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.2, "svg.fonttype": "none", "axes.linewidth": 0.7})


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["eligibility"] == "eligible"]


def summary(path: Path) -> dict:
    return json.loads(path.read_text())


def axes_style(ax: plt.Axes, grid: str = "both") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid, color="#DCE2E7", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def stamp(ax: plt.Axes, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(-0.13, 1.065, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
    ax.text(-0.01, 1.07, title, transform=ax.transAxes, fontsize=10.1, fontweight="bold", va="bottom")
    if subtitle:
        ax.text(-0.01, 1.005, subtitle, transform=ax.transAxes, fontsize=7.0, color=GREY, va="bottom")


def paired_representation(ax: plt.Axes, rows: list[dict[str, str]], letter: str, title: str) -> None:
    old = np.array([float(r["old_assembly_bp"]) for r in rows])
    new = np.array([float(r["new_assembly_bp"]) for r in rows])
    recovered = np.array([r["reference_state"] == "reference_collapse" for r in rows])
    limit = max(old.max(), new.max()) * 1.6 + 1
    ax.scatter(new[~recovered] + 1, old[~recovered] + 1, marker="s", s=29, facecolor="white", edgecolor=CHARCOAL, linewidth=1.0, label="other eligible family", zorder=2)
    ax.scatter(new[recovered] + 1, old[recovered] + 1, marker="o", s=45, color=BLUE, edgecolor=CHARCOAL, linewidth=0.75, label="preferentially recovered", zorder=3)
    ax.plot([1, limit], [1, limit], color="#A6B2BF", lw=1.0, ls="--", label="equal representation")
    ax.plot([1, limit], [0.6, 0.6 * limit], color=ORANGE, lw=1.15, ls=":", label="old/new = 0.6")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.8, limit); ax.set_ylim(0.8, limit)
    ax.set_xlabel("newer-assembly localized bp (+1)"); ax.set_ylabel("historical-assembly localized bp (+1)")
    axes_style(ax); ax.legend(frameon=False, fontsize=6.5, loc="upper left", handletextpad=0.45)
    stamp(ax, letter, title, "Every point is one read-derived family in the frozen primary denominator")


def deficit_gain(ax: plt.Axes, rows: list[dict[str, str]], pearson: float, letter: str, title: str) -> None:
    gain = np.array([float(r["observed_gain_bp"]) for r in rows])
    deficit = np.array([float(r["predicted_missing_bp"]) for r in rows])
    recovered = np.array([r["reference_state"] == "reference_collapse" for r in rows])
    limit = max(gain.max(), deficit.max()) * 1.6 + 1
    ax.scatter(gain[~recovered] + 1, deficit[~recovered] + 1, marker="s", s=29, facecolor="white", edgecolor=CHARCOAL, linewidth=1.0, label="other eligible family", zorder=2)
    ax.scatter(gain[recovered] + 1, deficit[recovered] + 1, marker="o", s=45, color=BLUE, edgecolor=CHARCOAL, linewidth=0.75, label="preferentially recovered", zorder=3)
    ax.plot([1, limit], [1, limit], color=GREY, lw=1.0, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.8, limit); ax.set_ylim(0.8, limit)
    ax.set_xlabel("newer - historical assembly gain (bp +1)"); ax.set_ylabel("HiFi - historical abundance deficit (bp +1)")
    axes_style(ax); ax.legend(frameon=False, fontsize=6.5, loc="upper left", handletextpad=0.45)
    ax.text(0.035, 0.065, f"Pearson r = {pearson:.3f}\nMagnitude agreement is descriptive", transform=ax.transAxes, fontsize=6.9, color=GREY)
    stamp(ax, letter, title, "Deficit and assembly gain retain separate meanings")


def system_design(ax: plt.Axes, ey: dict, mac: dict) -> None:
    ax.set_axis_off(); stamp(ax, "A", "Two historical-to-newer assembly comparisons", "Primary denominator: >=15 kb localized in the newer assembly")
    records = [
        ("Ey15-2", "CLR-Canu", "HiFi-Hifiasm", ey["eligible_family_rows"], ey["confusion"]["TP"]),
        ("Macadamia", "CLR Falcon-Unzip", "HiFi IPA", mac["eligible_family_rows"], mac["confusion"]["TP"]),
    ]
    for y, (name, old, new, eligible, recovered) in zip([0.60, 0.22], records):
        ax.text(0.02, y + 0.09, name, transform=ax.transAxes, fontsize=10.4, fontweight="bold", color="#244B72")
        for x, text, face in [(0.29, old, "#FFF1DE"), (0.54, new, "#EAF4EB")]:
            ax.add_patch(plt.Rectangle((x, y), 0.18, 0.19, transform=ax.transAxes, facecolor=face, edgecolor="#AAB7C3", lw=0.8))
            ax.text(x + 0.09, y + 0.115, text, transform=ax.transAxes, ha="center", va="center", fontsize=7.0, fontweight="bold")
            ax.text(x + 0.09, y + 0.055, "assembly representation", transform=ax.transAxes, ha="center", va="center", fontsize=6.2, color=GREY)
        ax.annotate("", xy=(0.53, y + 0.095), xytext=(0.48, y + 0.095), xycoords=ax.transAxes, arrowprops={"arrowstyle": "->", "lw": 1.15, "color": GREY})
        ax.text(0.77, y + 0.115, f"{eligible} eligible", transform=ax.transAxes, fontsize=8.0, fontweight="bold")
        ax.text(0.77, y + 0.045, f"{recovered} recovered", transform=ax.transAxes, fontsize=7.2, color=BLUE)
    ax.text(0.02, 0.02, "Newer assemblies are reference proxies, not absolute biological copy-number truth.", transform=ax.transAxes, fontsize=6.9, color=GREY)


def exact_match(ax: plt.Axes, ey: dict, mac: dict) -> None:
    ax.set_axis_off(); stamp(ax, "F", "Read prioritization matches preferential recovery", "Complete match in each frozen primary family set")
    for y, (label, summary) in zip([0.68, 0.28], [("Ey15-2", ey), ("Macadamia", mac)]):
        total = summary["eligible_family_rows"]
        hits = summary["confusion"]["TP"]
        others = total - hits
        ax.text(0.02, y + 0.09, label, transform=ax.transAxes, fontsize=9.2, fontweight="bold", color="#244B72")
        ax.add_patch(plt.Rectangle((0.30, y), 0.19, 0.20, transform=ax.transAxes, facecolor="#E8F2FB", edgecolor=BLUE, lw=0.9))
        ax.text(0.395, y + 0.125, f"{hits}/{hits}", transform=ax.transAxes, ha="center", va="center", fontsize=13, color=BLUE, fontweight="bold")
        ax.text(0.395, y + 0.050, "read-prioritized\n& recovered", transform=ax.transAxes, ha="center", va="center", fontsize=6.3)
        ax.add_patch(plt.Rectangle((0.58, y), 0.19, 0.20, transform=ax.transAxes, facecolor="#F5F6F7", edgecolor="#AAB7C3", lw=0.8))
        ax.text(0.675, y + 0.125, f"{others}/{others}", transform=ax.transAxes, ha="center", va="center", fontsize=12, color=CHARCOAL, fontweight="bold")
        ax.text(0.675, y + 0.050, "not prioritized\n& not recovered", transform=ax.transAxes, ha="center", va="center", fontsize=6.3)
    ax.text(0.02, 0.01, "This is a retrospective reference-proxy result, not an independent population accuracy estimate.", transform=ax.transAxes, fontsize=6.8, color=GREY)


def main() -> None:
    ey, mac = table(EY_METRICS), table(MAC_METRICS)
    ey_summary, mac_summary = summary(EY_SUMMARY), summary(MAC_SUMMARY)
    if (len(ey), len(mac), ey_summary["confusion"], mac_summary["confusion"]) != (19, 43, {"FN": 0, "FP": 0, "TN": 11, "TP": 8}, {"FN": 0, "FP": 0, "TN": 41, "TP": 2}):
        raise ValueError("Unexpected frozen primary family tables")
    fig = plt.figure(figsize=(11.4, 10.25), facecolor="white")
    grid = fig.add_gridspec(3, 2, left=0.09, right=0.985, bottom=0.082, top=0.92, hspace=0.62, wspace=0.38)
    system_design(fig.add_subplot(grid[0, 0]), ey_summary, mac_summary)
    paired_representation(fig.add_subplot(grid[0, 1]), ey, "B", "Ey15-2: representation in historical and newer assemblies")
    deficit_gain(fig.add_subplot(grid[1, 0]), ey, ey_summary["missing_bp_pearson"], "C", "Ey15-2: read deficit versus assembly gain")
    paired_representation(fig.add_subplot(grid[1, 1]), mac, "D", "Macadamia: representation in historical and newer assemblies")
    deficit_gain(fig.add_subplot(grid[2, 0]), mac, mac_summary["missing_bp_pearson"], "E", "Macadamia: read deficit versus assembly gain")
    exact_match(fig.add_subplot(grid[2, 1]), ey_summary, mac_summary)
    fig.text(0.5, 0.982, "Figure 4 | Improved assemblies preferentially recover read-prioritized repeat families", ha="center", va="top", fontsize=13.7, fontweight="bold")
    fig.text(0.5, 0.018, "Continuous values are read-assembly abundance deficits and assembly gains, not physical missing-base truth. Full family tables and denominator sensitivities remain in the evidence archive.", ha="center", fontsize=7.1, color=GREY)
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("png", {"dpi": 350})):
        fig.savefig(OUT / f"figure4_v2_assembly_recovery.{suffix}", bbox_inches="tight", facecolor="white", **kwargs)


if __name__ == "__main__":
    main()
