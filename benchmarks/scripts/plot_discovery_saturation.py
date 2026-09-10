"""Plot the formal TandemX discovery-saturation validation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmarks.challenge.schema import read_table, write_table


BLUE = "#2C6E9B"
GOLD = "#C58A2A"
GREEN = "#3A7D5D"
RED = "#A94A45"
GRAY = "#6B7280"


def grouped(rows: list[dict], field: str) -> tuple[list[float], list[float], list[float], list[float]]:
    by_depth: dict[float, list[float]] = {}
    for row in rows:
        by_depth.setdefault(float(row["coverage"]), []).append(float(row[field]))
    depths = sorted(by_depth)
    return (
        depths,
        [statistics.mean(by_depth[d]) for d in depths],
        [min(by_depth[d]) for d in depths],
        [max(by_depth[d]) for d in depths],
    )


def draw_band(axis, rows, field, label, color, marker="o"):
    x, mean, low, high = grouped(rows, field)
    axis.plot(x, mean, marker=marker, color=color, linewidth=1.8, label=label)
    axis.fill_between(x, low, high, color=color, alpha=0.15, linewidth=0)
    return x, mean


def run(evaluation: Path, outdir: Path) -> None:
    metrics = read_table(evaluation / "depth_metrics.tsv")
    transitions = read_table(evaluation / "transition_metrics.tsv")
    decision = json.loads((evaluation / "saturation_decision.json").read_text())
    outdir.mkdir(parents=True, exist_ok=False)

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.6), constrained_layout=True)
    axis = axes[0, 0]
    draw_band(axis, metrics, "unique_candidate_count", "Unique candidates", GRAY, "s")
    draw_band(axis, metrics, "operational_family_count", "Operational families", BLUE)
    axis.set_ylabel("Count")
    axis.set_title("A  Discovery yield")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[0, 1]
    transition_rows = [
        {**row, "coverage": row["higher_depth"]} for row in transitions
    ]
    draw_band(axis, transition_rows, "new_families_per_added_x", "New families per added 1x", GOLD)
    axis.set_ylabel("Families per added 1x", color=GOLD)
    twin = axis.twinx()
    draw_band(twin, transition_rows, "new_family_fraction", "New-family fraction", BLUE, "s")
    twin.axhline(0.05, color=BLUE, linestyle="--", linewidth=1)
    twin.set_ylabel("New-family fraction", color=BLUE)
    axis.set_title("B  Marginal discovery")

    axis = axes[1, 0]
    draw_band(axis, transition_rows, "family_jaccard", "Adjacent-depth Jaccard", GREEN)
    axis.axhline(0.95, color=RED, linestyle="--", linewidth=1, label="Saturation threshold")
    axis.set_ylim(0, 1.02)
    axis.set_ylabel("Family-level Jaccard")
    axis.set_title("C  Catalogue stability")
    axis.legend(frameon=False, fontsize=8)

    axis = axes[1, 1]
    draw_band(axis, metrics, "high_abundance_recall", "High", RED)
    draw_band(axis, metrics, "medium_abundance_recall", "Medium", GOLD, "s")
    draw_band(axis, metrics, "low_abundance_recall", "Low", BLUE, "^")
    axis.set_ylim(0, 1.02)
    axis.set_ylabel("Planted-family recall")
    axis.set_title("D  Recall by abundance tier")
    axis.legend(frameon=False, fontsize=8)

    for axis in axes.flat:
        axis.set_xscale("log", base=2)
        axis.set_xticks([0.5, 1, 2, 5, 10, 20, 30])
        axis.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "30"])
        axis.set_xlabel("Nominal read depth (x)")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.6)

    if decision["saturation_reached"]:
        depth = float(decision["saturation_depth"])
        for axis in axes.flat:
            axis.axvline(depth, color="#111827", linestyle=":", linewidth=1)
    fig.suptitle("TandemX discovery saturation across three validation genomes", fontsize=12)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"discovery_saturation.{suffix}", dpi=300)
    plt.close(fig)

    panel_rows = []
    for row in metrics:
        for panel, field in [
            ("A", "unique_candidate_count"),
            ("A", "operational_family_count"),
            ("D", "high_abundance_recall"),
            ("D", "medium_abundance_recall"),
            ("D", "low_abundance_recall"),
        ]:
            panel_rows.append(
                {
                    "panel": panel,
                    "seed": row["seed"],
                    "coverage": row["coverage"],
                    "metric": field,
                    "value": row[field],
                }
            )
    for row in transitions:
        for panel, field in [
            ("B", "new_families_per_added_x"),
            ("B", "new_family_fraction"),
            ("C", "family_jaccard"),
        ]:
            panel_rows.append(
                {
                    "panel": panel,
                    "seed": row["seed"],
                    "coverage": row["higher_depth"],
                    "metric": field,
                    "value": row[field],
                }
            )
    write_table(outdir / "panel_source.tsv", panel_rows, list(panel_rows[0]))
    legend = (
        "**Discovery saturation across three validation genomes.** A, unique "
        "candidate monomers and operational families. B, new operational families "
        "per added 1x and the fraction new relative to the higher-depth catalogue. "
        "C, one-to-one cyclic-sequence family Jaccard between adjacent nested depths. "
        "D, planted-family recall for low (20-copy), medium (80-copy) and high "
        "(at least 200-copy) families. Lines show three-seed means and ribbons show "
        "ranges. Dashed lines mark the pre-specified 5% and 0.95 transition criteria. "
        f"Primary saturation reached: {decision['saturation_reached']}; depth: "
        f"{decision['saturation_depth']}."
    )
    (outdir / "figure_legend.md").write_text(legend + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    run(args.evaluation, args.outdir)
