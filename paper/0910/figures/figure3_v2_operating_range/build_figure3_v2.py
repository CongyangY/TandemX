#!/usr/bin/env python3
"""Render Figure 3 from complete-library QC and nested TandemX runs.

The figure describes operating range.  It deliberately does not compare tools,
rank species, or imply that the public libraries are biological replicates.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).parent
QC_SOURCE = ROOT / "paper/evidence/multispecies_input_qc/figures_v2/panel_source.tsv"
RUN_SOURCE = ROOT / "paper/evidence/multispecies_real_diagnostics/figures_v2/panel_source.tsv"

PALETTE = {
    "A. thaliana": "#6B8EAD",
    "A. sativa": "#D48A45",
    "G. soja": "#7C6AA5",
    "H. vulgare": "#5D8A66",
    "O. sativa": "#C07AA9",
    "S. cereale": "#767676",
    "T. aestivum": "#B8A41D",
    "Z. mays": "#3C9DAF",
}
MATERIAL_LABELS = {
    "Chinese Spring": "Wheat\nChinese Spring",
    "Morex": "Barley\nMorex",
    "Lo7": "Rye\nLo7",
    "Victoria": "Oat\nVictoria",
}
MATERIAL_COLORS = {
    "Chinese Spring": "#B8A41D",
    "Morex": "#5D8A66",
    "Lo7": "#767676",
    "Nipponbare": "#C07AA9",
    "Victoria": "#D48A45",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.3,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "axes.labelcolor": "#27313B",
        "xtick.color": "#27313B",
        "ytick.color": "#27313B",
    }
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concise_species(value: str) -> str:
    first, second = value.split()[:2]
    return f"{first[0]}. {second}"


def axes_style(ax: plt.Axes, grid: str = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid, color="#D9DEE3", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def panel_title(ax: plt.Axes, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(-0.13, 1.055, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
    ax.text(-0.01, 1.06, title, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom")
    if subtitle:
        ax.text(-0.01, 0.995, subtitle, transform=ax.transAxes, fontsize=7.2, color="#5D6975", va="bottom")


def selected_tandemx_rows(rows: list[dict[str, str]], materials: list[str]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for material in materials:
        candidates = [r for r in rows if r["material"] == material and r["tool"] == "tandemx"]
        if not candidates:
            raise ValueError(f"No TandemX rows for {material}")
        selected.append(min(candidates, key=lambda r: abs(float(r["input_Gb"]) - 0.012)))
    return selected


def main() -> None:
    qc = read_tsv(QC_SOURCE)
    runs = read_tsv(RUN_SOURCE)
    libraries = [row for row in qc if row["panel"] == "A"]
    if len(libraries) != 10 or len({r["reported_species"] for r in libraries}) != 8:
        raise ValueError("Expected ten libraries from eight reported species")
    accession_order = [r["run_accession"] for r in libraries]
    medians = {
        r["run_accession"]: float(r["value"])
        for r in qc
        if r["panel"] == "B" and r["metric"] == "median_read_length_kb"
    }
    n50s = {
        r["run_accession"]: float(r["value"])
        for r in qc
        if r["panel"] == "B" and r["metric"] == "read_n50_kb"
    }
    if set(accession_order) != set(medians) or set(accession_order) != set(n50s):
        raise ValueError("Read-length summaries do not match the complete-library cohort")

    labels = [f'{r["reported_material"]}\n({concise_species(r["reported_species"])})' for r in libraries]
    colors = [PALETTE[concise_species(r["reported_species"])] for r in libraries]
    y = np.arange(len(libraries))
    fig = plt.figure(figsize=(11.25, 8.35), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.09, right=0.985, bottom=0.105, top=0.905, hspace=0.58, wspace=0.34
    )

    # A. Input breadth, retaining every enrolled library.
    ax_a = fig.add_subplot(outer[0, 0])
    bases = [float(r["value"]) for r in libraries]
    bars = ax_a.barh(y, bases, height=0.73, color=colors, edgecolor="none", zorder=2)
    ax_a.set_yticks(y, labels, fontsize=7.6)
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 86)
    ax_a.set_xlabel("complete-library bases (Gb)")
    axes_style(ax_a)
    for bar, value in zip(bars, bases):
        ax_a.text(value + 1.0, bar.get_y() + bar.get_height() / 2, f"{value:.1f}", va="center", fontsize=7.4)
    panel_title(ax_a, "A", "Ten complete public HiFi libraries", "Eight reported species; all files passed integrity and FASTQ checks")

    # B. Pair of read-length positions avoids a redundant legend-heavy bar chart.
    ax_b = fig.add_subplot(outer[0, 1])
    for idx, (accession, color) in enumerate(zip(accession_order, colors)):
        ax_b.plot([medians[accession], n50s[accession]], [idx, idx], color="#B9C2CC", lw=1.15, zorder=1)
        ax_b.scatter(medians[accession], idx, s=31, color=color, marker="o", edgecolor="white", linewidth=0.45, zorder=3)
        ax_b.scatter(n50s[accession], idx, s=36, color=color, marker="D", edgecolor="white", linewidth=0.45, zorder=3)
    ax_b.set_yticks(y, labels, fontsize=7.6)
    ax_b.invert_yaxis()
    ax_b.set_xlim(12.5, 23.7)
    ax_b.set_xlabel("read length (kb)")
    axes_style(ax_b)
    ax_b.scatter([], [], color="#526273", marker="o", s=28, label="median")
    ax_b.scatter([], [], color="#526273", marker="D", s=32, label="read N50")
    ax_b.legend(frameon=False, fontsize=7.2, ncol=2, handletextpad=0.45, loc="lower right")
    panel_title(ax_b, "B", "Read-length range", "Circles, median; diamonds, read N50")

    # C. Separate aligned scales so output count and support fraction are never read on one axis.
    materials = ["Chinese Spring", "Morex", "Lo7", "Victoria"]
    chosen = selected_tandemx_rows(runs, materials)
    inner_c = outer[1, 0].subgridspec(1, 2, wspace=0.12)
    ax_c1 = fig.add_subplot(inner_c[0, 0])
    ax_c2 = fig.add_subplot(inner_c[0, 1], sharey=ax_c1)
    cy = np.arange(len(chosen))
    c_labels = [MATERIAL_LABELS[r["material"]] for r in chosen]
    c_colors = [MATERIAL_COLORS[r["material"]] for r in chosen]
    density = [float(r["observed_call_density_per_Mbp"]) for r in chosen]
    fraction = [100 * float(r["observed_positive_read_fraction"]) for r in chosen]
    ax_c1.hlines(cy, 0, density, color="#CBD4DC", linewidth=1.2, zorder=1)
    ax_c1.scatter(density, cy, s=54, color=c_colors, edgecolor="white", linewidth=0.65, zorder=3)
    ax_c1.set_yticks(cy, c_labels, fontsize=7.8)
    ax_c1.invert_yaxis()
    ax_c1.set_xlim(0, 82)
    ax_c1.set_xlabel("calls per Mb")
    axes_style(ax_c1)
    ax_c1.grid(axis="y", visible=False)
    ax_c2.hlines(cy, 0, fraction, color="#CBD4DC", linewidth=1.2, zorder=1)
    ax_c2.scatter(fraction, cy, s=54, color=c_colors, edgecolor="white", linewidth=0.65, zorder=3)
    ax_c2.set_xlim(0, 80)
    ax_c2.set_xlabel("positive reads (%)")
    axes_style(ax_c2)
    ax_c2.grid(axis="y", visible=False)
    ax_c2.tick_params(axis="y", left=False, labelleft=False)
    panel_title(ax_c1, "C", "Observed repeat-call landscape", "Matched 11.4–12.7 Mb nested inputs; TandemX outputs")

    # D. Nested inputs are explicitly connected only within material.
    ax_d = fig.add_subplot(outer[1, 1])
    for material in ["Chinese Spring", "Morex", "Lo7", "Nipponbare", "Victoria"]:
        rows = sorted(
            [r for r in runs if r["material"] == material and r["tool"] == "tandemx"],
            key=lambda r: float(r["input_Gb"]),
        )
        values = [100 * float(r["observed_union_base_fraction"]) for r in rows]
        inputs = [float(r["input_Gb"]) for r in rows]
        ax_d.plot(inputs, values, marker="o", ms=4.3, lw=1.45, color=MATERIAL_COLORS[material], label=material)
    ax_d.set_xscale("log")
    ax_d.set_xlim(0.009, 1.45)
    ax_d.set_ylim(2.35, 9.7)
    ax_d.set_xlabel("nested input (Gb)")
    ax_d.set_ylabel("called-base fraction (%)")
    axes_style(ax_d, grid="both")
    ax_d.legend(frameon=False, fontsize=6.8, ncol=2, loc="upper right", handlelength=1.8, columnspacing=1.3)
    panel_title(ax_d, "D", "Called repeat-base fraction across nested inputs", "Within-library nesting; includes the 1.17-Gb barley run")

    fig.text(0.5, 0.992, "Figure 3 | TandemX operates across diverse plant long-read datasets", ha="center", va="top", fontsize=14.5, fontweight="bold")
    fig.text(
        0.5,
        0.028,
        "Complete libraries establish the input range. Nested runs are descriptive operating data, not cross-species accuracy estimates or biological replicates.",
        ha="center",
        va="center",
        fontsize=7.7,
        color="#5D6975",
    )
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("png", {"dpi": 350})):
        fig.savefig(OUT / f"figure3_v2_operating_range.{suffix}", bbox_inches="tight", facecolor="white", **kwargs)


if __name__ == "__main__":
    main()
