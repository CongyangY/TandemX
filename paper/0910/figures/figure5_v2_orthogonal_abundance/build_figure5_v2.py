#!/usr/bin/env python3
"""Render Figure 5 from the frozen orthogonal-abundance evidence tables."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).parent
ORTHO = ROOT / "paper/evidence/orthogonal_abundance_validation_v1/figures_v2/panel_source.tsv"
QC = ROOT / "paper/0910/source_data/figure6_txf000695_qc.tsv"
BLUE = "#2D76B7"
ORANGE = "#D57A1F"
CHARCOAL = "#39424E"
GREY = "#6D7A87"
LIGHT = "#D8E0E8"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.1, "svg.fonttype": "none", "axes.linewidth": 0.7})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float | None:
    return float(value) if value.strip() else None


def axes_style(ax: plt.Axes, grid: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid, color="#DCE2E7", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def stamp(ax: plt.Axes, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(-0.13, 1.075, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
    ax.text(-0.01, 1.08, title, transform=ax.transAxes, fontsize=10.2, fontweight="bold", va="bottom")
    if subtitle:
        ax.text(-0.01, 1.012, subtitle, transform=ax.transAxes, fontsize=7.1, color=GREY, va="bottom")


def measure_rows(rows: list[dict[str, str]], family: str) -> list[dict[str, str]]:
    return [r for r in rows if r["panel"] == "A" and r["family_id"] == family]


def get_measure(rows: list[dict[str, str]], name: str) -> dict[str, str]:
    matches = [r for r in rows if r["platform_or_measure"] == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one {name!r}; found {len(matches)}")
    return matches[0]


def abundance_panel(ax: plt.Axes, rows: list[dict[str, str]], families: list[str], title: str, letter: str, ont: bool = False) -> None:
    labels = {"TXF000002": "TXF000002", "TXF000154": "TXF000154", "TXF000496": "TXF000496"}
    x = np.arange(len(families))
    positions = {"Newer assembly": -0.27, "Frozen HiFi": -0.09, "Illumina k=21": 0.09, "Illumina k=31": 0.27}
    styles = {
        "Newer assembly": ("s", "white", CHARCOAL),
        "Frozen HiFi": ("D", CHARCOAL, CHARCOAL),
        "Illumina k=21": ("o", BLUE, BLUE),
        "Illumina k=31": ("o", "white", BLUE),
    }
    for idx, family in enumerate(families):
        family_rows = measure_rows(rows, family)
        for name, offset in positions.items():
            row = get_measure(family_rows, name)
            marker, fill, edge = styles[name]
            ax.scatter(idx + offset, float(row["value_bp"]), marker=marker, s=46, facecolor=fill, edgecolor=edge, linewidth=1.15, zorder=3)
        if ont:
            row = get_measure(family_rows, "ONT direction range")
            low, high = float(row["value_bp_min"]), float(row["value_bp_max"])
            ax.vlines(idx + 0.43, low, high, color=ORANGE, lw=2.3, zorder=2)
            ax.hlines([low, high], idx + 0.38, idx + 0.48, color=ORANGE, lw=1.2, zorder=2)
    ax.set_yscale("log")
    ax.set_xticks(x, [labels[f] for f in families])
    ax.set_ylabel("estimated abundance or localized bp")
    ax.set_ylim(8e3, 3e6)
    axes_style(ax)
    stamp(ax, letter, title, "Assembly square; HiFi diamond; Illumina circles" + ("; ONT range" if ont else ""))
    if ont:
        ax.text(0.98, 0.04, "ONT is direction-level evidence", ha="right", transform=ax.transAxes, fontsize=6.6, color=GREY)


def main() -> None:
    rows = read_tsv(ORTHO)
    qc = read_tsv(QC)
    states = {r["family_id"]: r["final_interpretation"] for r in rows if r["panel"] == "C" and r["platform_or_measure"] == "Final"}
    expected = {
        "TXF000002": "orthogonal_supports_residual_collapse",
        "TXF000154": "orthogonal_supports_residual_collapse",
        "TXF001517": "unresolved",
        "TXF000496": "orthogonal_supports_residual_collapse",
        "TXF000563": "unresolved",
        "TXF000695": "unresolved",
    }
    if states != expected or len(qc) != 2:
        raise ValueError("Unexpected frozen orthogonal evidence state")

    fig = plt.figure(figsize=(12.0, 8.65), facecolor="white")
    grid = fig.add_gridspec(3, 3, left=0.075, right=0.985, bottom=0.095, top=0.90, hspace=0.78, wspace=0.48, height_ratios=[0.78, 1.08, 1.0])

    # A: the preselected evidence design.
    ax_a = fig.add_subplot(grid[0, :])
    ax_a.set_axis_off()
    stamp(ax_a, "A", "Preselected candidates were assessed with orthogonal reads", "Six newer-assembly families were fixed before independent-read inspection")
    groups = [("Ey15-2", ["TXF000002", "TXF000154", "TXF001517"], "Illumina k=21 and k=31"), ("Macadamia", ["TXF000496", "TXF000563", "TXF000695"], "Illumina k=21 and k=31; ONT direction")]
    for y, (system, families, evidence) in zip([0.57, 0.18], groups):
        ax_a.text(0.02, y + 0.08, system, color="#244B72", fontweight="bold", fontsize=10, transform=ax_a.transAxes)
        for j, family in enumerate(families):
            x = 0.18 + j * 0.13
            outcome = states[family]
            face = "#E9F2FB" if outcome.startswith("orthogonal_supports") else "#F4F4F4"
            ax_a.add_patch(plt.Rectangle((x, y), 0.11, 0.20, transform=ax_a.transAxes, facecolor=face, edgecolor="#AAB6C1", lw=0.8))
            ax_a.text(x + 0.055, y + 0.10, family, ha="center", va="center", fontsize=7.1, transform=ax_a.transAxes)
        ax_a.annotate("", xy=(0.64, y + 0.10), xytext=(0.59, y + 0.10), xycoords=ax_a.transAxes, arrowprops={"arrowstyle": "->", "lw": 1.15, "color": GREY})
        ax_a.text(0.67, y + 0.125, evidence, fontsize=8.0, transform=ax_a.transAxes)
        ax_a.text(0.67, y + 0.035, "assembly/read ratio < 0.6 at both k values", fontsize=7.0, color=GREY, transform=ax_a.transAxes)
    ax_a.text(0.965, 0.45, "3 supported\n3 unresolved", ha="right", va="center", fontsize=11, fontweight="bold", color=BLUE, transform=ax_a.transAxes)

    abundance_panel(fig.add_subplot(grid[1, 0]), rows, ["TXF000002", "TXF000154"], "Ey15-2 supported families", "B")
    abundance_panel(fig.add_subplot(grid[1, 1]), rows, ["TXF000496"], "Macadamia supported family", "C", ont=True)

    # D: ratios preserve the final classifier boundary without translating it to physical missing bases.
    ax_d = fig.add_subplot(grid[1, 2])
    unresolved = ["TXF001517", "TXF000563", "TXF000695"]
    x = np.arange(3)
    for i, family in enumerate(unresolved):
        family_rows = [r for r in rows if r["panel"] == "B" and r["family_id"] == family]
        for label, marker, fill, edge, offset in [("k=21", "o", BLUE, BLUE, -0.12), ("k=31", "o", "white", BLUE, 0.12)]:
            row = get_measure(family_rows, f"Illumina {label}")
            ax_d.scatter(i + offset, float(row["assembly_read_ratio"]), marker=marker, s=45, facecolor=fill, edgecolor=edge, linewidth=1.15, zorder=3)
        ont_rows = [r for r in family_rows if r["platform_or_measure"] == "ONT direction range"]
        if ont_rows:
            low, high = float(ont_rows[0]["assembly_read_ratio_min"]), float(ont_rows[0]["assembly_read_ratio_max"])
            ax_d.vlines(i + 0.30, low, high, color=ORANGE, lw=2)
            ax_d.hlines([low, high], i + 0.25, i + 0.35, color=ORANGE, lw=1)
    ax_d.axhline(0.6, color=CHARCOAL, lw=1.0, ls="--", zorder=1)
    ax_d.set_yscale("log"); ax_d.set_ylim(0.01, 30); ax_d.set_xticks(x, unresolved); ax_d.set_ylabel("newer assembly / read estimate")
    axes_style(ax_d)
    stamp(ax_d, "D", "Three candidates remain unresolved", "Illumina k=21 filled; k=31 open; orange, ONT range")
    ax_d.text(0.02, 0.44, "frozen threshold = 0.6", transform=ax_d.transAxes, fontsize=6.8, color=CHARCOAL)

    # E: exact TXF000695 depth-bin counts make instability inspectable.
    ax_e = fig.add_subplot(grid[2, 0:2])
    bins = ["<2", "2–99", "100–999", "1k–9,999", "50k–99,999"]
    k21 = [int(qc[0][key]) for key in ["depth_below_2_count", "depth_below_100_count", "depth_below_1000_count", "depth_below_10000_count", "depth_50000_to_99999_count"]]
    k31 = [int(qc[1][key]) for key in ["depth_below_2_count", "depth_below_100_count", "depth_below_1000_count", "depth_below_10000_count", "depth_50000_to_99999_count"]]
    # cumulative values become non-overlapping histogram bins.
    k21 = [k21[0], k21[1] - k21[0], k21[2] - k21[1], k21[3] - k21[2], k21[4]]
    k31 = [k31[0], k31[1] - k31[0], k31[2] - k31[1], k31[3] - k31[2], k31[4]]
    bx = np.arange(len(bins)); width = 0.36
    ax_e.bar(bx - width / 2, k21, width, color=BLUE, label="k=21")
    ax_e.bar(bx + width / 2, k31, width, facecolor="white", edgecolor=BLUE, linewidth=1.2, label="k=31")
    ax_e.set_xticks(bx, bins); ax_e.set_ylabel("diagnostic words")
    axes_style(ax_e)
    ax_e.legend(frameon=False, fontsize=7.2, ncol=2, loc="upper left")
    stamp(ax_e, "E", "TXF000695 has a cross-k median-boundary shift", "467 diagnostic words per k; low and high depth components are retained")
    ax_e.text(0.985, 0.90, "median depth: 77,091 → 101\nabundance: 1.10 Mb → 1.65 kb", transform=ax_e.transAxes, ha="right", va="top", fontsize=7.1, color=GREY)

    # F: categorical disposition, not a numerical population estimate.
    ax_f = fig.add_subplot(grid[2, 2])
    all_families = list(expected)
    for idx, family in enumerate(all_families):
        supported = states[family].startswith("orthogonal_supports")
        ax_f.scatter(0, 5 - idx, s=110, marker="o" if supported else "D", color=BLUE if supported else "white", edgecolor=BLUE if supported else CHARCOAL, linewidth=1.2, zorder=3)
        ax_f.text(0.16, 5 - idx, family, va="center", fontsize=7.8)
        ax_f.text(1.62, 5 - idx, "supported" if supported else "unresolved", va="center", ha="right", fontsize=7.4, color=BLUE if supported else CHARCOAL)
    ax_f.set_xlim(-0.24, 1.78); ax_f.set_ylim(-0.55, 5.6); ax_f.set_xticks([]); ax_f.set_yticks([])
    ax_f.spines[:].set_visible(False)
    stamp(ax_f, "F", "Evidence-state boundary", "Circle, cross-platform support; diamond, unresolved")
    ax_f.text(0.0, -0.12, "A family-level disposition, not a population accuracy estimate.", transform=ax_f.transAxes, fontsize=6.7, color=GREY)

    fig.text(0.5, 0.986, "Figure 5 | Independent reads support three residual family-specific deficits and define three unresolved cases", ha="center", va="top", fontsize=13.5, fontweight="bold")
    fig.text(0.5, 0.022, "Read–assembly abundance ratios and ranges show direction-level evidence. They do not measure physical missing bases or establish absolute repeat copy number.", ha="center", fontsize=7.3, color=GREY)
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("png", {"dpi": 350})):
        fig.savefig(OUT / f"figure5_v2_orthogonal_abundance.{suffix}", bbox_inches="tight", facecolor="white", **kwargs)


if __name__ == "__main__":
    main()
