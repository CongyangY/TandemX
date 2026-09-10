#!/usr/bin/env python3
"""Plot the frozen three-panel orthogonal abundance validation figure."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from benchmarks.challenge.schema import digest_file


BLUE = "#2A6FBB"
BLUE_LIGHT = "#CFE0F2"
ORANGE = "#D97706"
ORANGE_LIGHT = "#F5D8B5"
CHARCOAL = "#30343B"
GREY = "#9299A1"
LIGHT_GREY = "#E4E7EA"
GRID = "#D9DDE2"

SOURCE_FIELDS = (
    "panel",
    "species_key",
    "family_id",
    "platform_or_measure",
    "value_bp",
    "value_bp_min",
    "value_bp_max",
    "assembly_read_ratio",
    "assembly_read_ratio_min",
    "assembly_read_ratio_max",
    "evidence_state",
    "final_interpretation",
    "warning",
)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def style_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(CHARCOAL)
    axis.spines["bottom"].set_color(CHARCOAL)
    axis.tick_params(colors=CHARCOAL, labelsize=7.5)
    axis.grid(True, axis="y", color=GRID, linewidth=0.55, alpha=0.75, zorder=0)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.05,
        label,
        transform=axis.transAxes,
        fontsize=11,
        fontweight="bold",
        color=CHARCOAL,
        va="top",
    )


def x_label(row: dict[str, str]) -> str:
    species = "Ey15" if row["species_key"] == "ey15" else "Mac"
    return f"{species}\n{row['family_id']}"


def matrix_label(row: dict[str, str]) -> str:
    species = "Ey15-2" if row["species_key"] == "ey15" else "Macadamia"
    return f"{species}\n{row['family_id']}"


def ont_candidate_rows(path: Path) -> dict[str, dict[str, str]]:
    return {
        row["family_id"]: row
        for row in load_tsv(path)
        if row["preorthogonal_candidate"].lower() == "true"
    }


def evidence_code(state: str) -> str:
    if "supports_residual" in state:
        return "support"
    if "supports_quantification_bias" in state or "does_not_support_residual" in state:
        return "bias_direction"
    if state == "not_available":
        return "not_available"
    return "unresolved"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ey15-final", required=True, type=Path)
    parser.add_argument("--macadamia-final", required=True, type=Path)
    parser.add_argument("--macadamia-ont", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    rows = load_tsv(args.ey15_final) + load_tsv(args.macadamia_final)
    ont_rows = ont_candidate_rows(args.macadamia_ont)
    if len(rows) != 6 or len({row["family_id"] for row in rows}) != 6:
        raise ValueError("expected exactly six unique frozen candidate families")
    final_counts: dict[str, int] = {}
    for row in rows:
        state = row["final_interpretation"]
        final_counts[state] = final_counts.get(state, 0) + 1
    if final_counts != {
        "orthogonal_supports_residual_collapse": 3,
        "unresolved": 3,
    }:
        raise ValueError(f"unexpected final interpretations: {final_counts}")

    labels = [x_label(row) for row in rows]
    matrix_labels = [matrix_label(row) for row in rows]
    xvalues = list(range(len(rows)))
    source_rows: list[dict[str, Any]] = []

    matplotlib.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titleweight": "bold",
            "axes.titlesize": 9,
            "axes.labelsize": 8,
        }
    )
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(12.2, 4.15),
        gridspec_kw={"width_ratios": [1.22, 1.14, 1.35]},
    )

    # A: absolute abundance.
    axis = axes[0]
    abundance_specs = (
        ("Newer assembly", "newer_assembly_bp", -0.24, "s", CHARCOAL, "white"),
        ("Frozen HiFi", "frozen_hifi_estimated_bp", -0.08, "D", CHARCOAL, CHARCOAL),
        ("Illumina k=21", "k21_orthogonal_estimated_bp", 0.08, "o", BLUE, BLUE),
        ("Illumina k=31", "k31_orthogonal_estimated_bp", 0.24, "o", BLUE, "white"),
    )
    for name, field, offset, marker, edge, face in abundance_specs:
        values = [float(row[field]) for row in rows]
        if not all(math.isfinite(value) and value > 0 for value in values):
            raise ValueError(f"non-positive or non-finite {field}")
        axis.scatter(
            [x + offset for x in xvalues],
            values,
            marker=marker,
            s=34,
            facecolors=face,
            edgecolors=edge,
            linewidths=1.0,
            zorder=4,
            label=name,
        )
        for row, value in zip(rows, values):
            source_rows.append(
                {
                    "panel": "A",
                    "species_key": row["species_key"],
                    "family_id": row["family_id"],
                    "platform_or_measure": name,
                    "value_bp": value,
                    "final_interpretation": row["final_interpretation"],
                    "warning": "abundance_estimate_not_physical_copy_truth",
                }
            )
    for x, row in zip(xvalues, rows):
        if row["species_key"] != "macadamia":
            continue
        ont = ont_rows[row["family_id"]]
        low = float(ont["calibrated_ONT_bp_at_616000000"])
        high = float(ont["calibrated_ONT_bp_at_780000000"])
        axis.vlines(x + 0.38, low, high, color=ORANGE, linewidth=2.0, zorder=3)
        axis.scatter(
            [x + 0.38],
            [math.sqrt(low * high)],
            marker="_",
            s=55,
            color=ORANGE,
            linewidths=1.5,
            zorder=4,
        )
        source_rows.append(
            {
                "panel": "A",
                "species_key": row["species_key"],
                "family_id": row["family_id"],
                "platform_or_measure": "ONT direction range",
                "value_bp_min": low,
                "value_bp_max": high,
                "evidence_state": ont["ONT_direction_result"],
                "final_interpretation": row["final_interpretation"],
                "warning": "direction_only_genome_size_sensitivity_not_copy_truth",
            }
        )
    axis.set_yscale("log")
    axis.set_xticks(xvalues, labels, rotation=27, ha="right")
    axis.set_ylabel("Estimated abundance or localized bp (log scale)")
    axis.set_title("Read and assembly abundance")
    axis.legend(
        handles=[
            Line2D([], [], marker="s", markerfacecolor="white", markeredgecolor=CHARCOAL, linestyle="", label="Newer assembly"),
            Line2D([], [], marker="D", color=CHARCOAL, linestyle="", label="Frozen HiFi"),
            Line2D([], [], marker="o", color=BLUE, linestyle="", label="Illumina k=21"),
            Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE, linestyle="", label="Illumina k=31"),
            Line2D([], [], marker="_", color=ORANGE, linestyle="", label="ONT range"),
        ],
        frameon=False,
        fontsize=6.5,
        ncol=2,
        loc="lower left",
    )
    style_axis(axis)
    panel_label(axis, "A")

    # B: assembly/read ratio and frozen threshold.
    axis = axes[1]
    ratio_specs = (
        ("Illumina k=21", "k21_orthogonal_estimated_bp", -0.10, "o", BLUE, BLUE),
        ("Illumina k=31", "k31_orthogonal_estimated_bp", 0.10, "o", BLUE, "white"),
    )
    for name, field, offset, marker, edge, face in ratio_specs:
        ratios = [float(row["newer_assembly_bp"]) / float(row[field]) for row in rows]
        axis.scatter(
            [x + offset for x in xvalues],
            ratios,
            marker=marker,
            s=35,
            facecolors=face,
            edgecolors=edge,
            linewidths=1.0,
            zorder=4,
            label=name,
        )
        for row, ratio in zip(rows, ratios):
            source_rows.append(
                {
                    "panel": "B",
                    "species_key": row["species_key"],
                    "family_id": row["family_id"],
                    "platform_or_measure": name,
                    "assembly_read_ratio": ratio,
                    "final_interpretation": row["final_interpretation"],
                    "warning": "ratio_below_threshold_supports_direction_not_physical_missing_bp",
                }
            )
    for x, row in zip(xvalues, rows):
        if row["species_key"] != "macadamia":
            continue
        ont = ont_rows[row["family_id"]]
        low = min(
            float(ont["newer_assembly_ONT_ratio_at_616000000"]),
            float(ont["newer_assembly_ONT_ratio_at_780000000"]),
        )
        high = max(
            float(ont["newer_assembly_ONT_ratio_at_616000000"]),
            float(ont["newer_assembly_ONT_ratio_at_780000000"]),
        )
        axis.vlines(x + 0.27, low, high, color=ORANGE, linewidth=2.0, zorder=3)
        axis.scatter([x + 0.27], [math.sqrt(low * high)], marker="_", s=55, color=ORANGE, linewidths=1.5, zorder=4)
        source_rows.append(
            {
                "panel": "B",
                "species_key": row["species_key"],
                "family_id": row["family_id"],
                "platform_or_measure": "ONT direction range",
                "assembly_read_ratio_min": low,
                "assembly_read_ratio_max": high,
                "evidence_state": ont["ONT_direction_result"],
                "final_interpretation": row["final_interpretation"],
                "warning": "direction_only_genome_size_sensitivity_not_copy_truth",
            }
        )
    axis.axhline(0.6, color=CHARCOAL, linestyle="--", linewidth=1.0, zorder=2)
    axis.set_yscale("log")
    axis.set_xticks(xvalues, labels, rotation=27, ha="right")
    axis.set_ylabel("Newer assembly / read abundance (log scale)")
    axis.set_title("Direction and threshold stability")
    axis.legend(
        handles=[
            Line2D([], [], marker="o", color=BLUE, linestyle="", label="Illumina k=21"),
            Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE, linestyle="", label="Illumina k=31"),
            Line2D([], [], color=CHARCOAL, linestyle="--", linewidth=1.0, label="Frozen threshold 0.6"),
        ],
        frameon=False,
        fontsize=6.6,
        loc="upper left",
    )
    style_axis(axis)
    panel_label(axis, "B")

    # C: platform evidence states and final interpretation.
    axis = axes[2]
    columns = ("HiFi deficit", "Illumina k=21", "Illumina k=31", "ONT", "Final")
    state_style = {
        "support": (BLUE, CHARCOAL, "o", "S"),
        "bias_direction": (ORANGE_LIGHT, ORANGE, "s", "B"),
        "unresolved": (LIGHT_GREY, CHARCOAL, "D", "U"),
        "not_available": ("white", GREY, "x", "--"),
    }
    for row_index, row in enumerate(rows):
        raw_states = [
            "orthogonal_supports_residual_collapse_at_this_k",
            row["k21_result"],
            row["k31_result"],
            row["ont_direction_result"],
            row["final_interpretation"],
        ]
        for column_index, raw_state in enumerate(raw_states):
            code = evidence_code(raw_state)
            face, edge, marker, text_value = state_style[code]
            scatter_options: dict[str, Any] = {
                "marker": marker,
                "s": 90 if marker != "x" else 50,
                "linewidths": 1.0,
                "zorder": 3,
            }
            if marker == "x":
                scatter_options["color"] = edge
            else:
                scatter_options["facecolors"] = face
                scatter_options["edgecolors"] = edge
            axis.scatter([column_index], [row_index], **scatter_options)
            if marker != "x":
                axis.text(column_index, row_index, text_value, ha="center", va="center", fontsize=5.8, color=CHARCOAL, zorder=4)
            source_rows.append(
                {
                    "panel": "C",
                    "species_key": row["species_key"],
                    "family_id": row["family_id"],
                    "platform_or_measure": columns[column_index],
                    "evidence_state": raw_state,
                    "final_interpretation": row["final_interpretation"],
                    "warning": "binary_state_separate_from_magnitude",
                }
            )
    axis.axvline(3.5, color=GRID, linewidth=1.0)
    axis.set_xlim(-0.55, 4.55)
    axis.set_ylim(len(rows) - 0.5, -0.5)
    axis.set_xticks(range(len(columns)), columns, rotation=25, ha="right")
    axis.set_yticks(range(len(rows)), matrix_labels)
    axis.set_title("Evidence concordance and final state")
    axis.set_ylabel("Frozen candidate")
    axis.grid(True, color=GRID, linewidth=0.55, alpha=0.75, zorder=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.tick_params(colors=CHARCOAL, labelsize=7.2, length=0)
    panel_label(axis, "C")

    figure.legend(
        handles=[
            Line2D([], [], marker="o", markerfacecolor=BLUE, markeredgecolor=CHARCOAL, linestyle="", label="S  residual support"),
            Line2D([], [], marker="s", markerfacecolor=ORANGE_LIGHT, markeredgecolor=ORANGE, linestyle="", label="B  assembly/bias direction"),
            Line2D([], [], marker="D", markerfacecolor=LIGHT_GREY, markeredgecolor=CHARCOAL, linestyle="", label="U  unresolved"),
            Line2D([], [], marker="x", color=GREY, linestyle="", label="--  unavailable"),
        ],
        frameon=False,
        fontsize=6.8,
        loc="lower center",
        bbox_to_anchor=(0.72, 0.075),
        ncol=4,
    )

    figure.text(
        0.06,
        0.015,
        "Six candidates were fixed before orthogonal inspection. Ratios and ranges are read--assembly evidence, not physical missing-bp truth; 3 supported and 3 unresolved is not a population accuracy estimate.",
        fontsize=7.0,
        color=CHARCOAL,
        ha="left",
    )
    figure.subplots_adjust(left=0.07, right=0.99, bottom=0.31, top=0.88, wspace=0.34)

    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = args.outdir / "orthogonal_abundance_validation"
    svg = stem.with_suffix(".svg")
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    figure.savefig(svg, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=240, bbox_inches="tight")
    plt.close(figure)

    panel_source = args.outdir / "panel_source.tsv"
    with panel_source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SOURCE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(source_rows)

    legend = args.outdir / "figure_legend.md"
    legend.write_text(
        "**Orthogonal validation of read-based abundance and residual assembly "
        "under-representation.** A, newer-assembly localized bases, frozen HiFi "
        "abundance, independent Illumina k=21/k=31 abundance and Macadamia "
        "direction-only ONT ranges for six predeclared candidates. B, newer-"
        "assembly/read ratios with the frozen 0.6 threshold; ONT ranges span "
        "the 616--780-Mb genome-size sensitivity. C, platform-level direction "
        "and final frozen interpretation. Ey15 TXF000002/TXF000154 and "
        "Macadamia TXF000496 support residual under-representation; the other "
        "three families are unresolved. The selected 3/6 split is not an "
        "accuracy rate, and continuous values are read--assembly abundance "
        "deficits rather than physical missing-base truth.\n",
        encoding="utf-8",
    )

    svg_text = svg.read_text(encoding="utf-8")
    receipt = {
        "schema_version": 1,
        "status": "complete_pending_visual_inspection",
        "panel_count": 3,
        "candidate_count": len(rows),
        "final_interpretation_counts": final_counts,
        "panel_source_rows": len(source_rows),
        "svg_text_element_count": svg_text.count("<text"),
        "svg_raster_image_element_count": svg_text.count("<image"),
        "sources": {
            "plotting_script": {
                "path": str(Path(__file__)),
                "sha256": digest_file(Path(__file__)),
            },
            "ey15_final": {"path": str(args.ey15_final), "sha256": digest_file(args.ey15_final)},
            "macadamia_final": {"path": str(args.macadamia_final), "sha256": digest_file(args.macadamia_final)},
            "macadamia_ont": {"path": str(args.macadamia_ont), "sha256": digest_file(args.macadamia_ont)},
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": digest_file(path)}
            for path in (svg, pdf, png, panel_source, legend)
        },
        "boundary": "binary_interpretation_separate_from_continuous_magnitude;read_assembly_deficit_not_physical_missing_bp_truth",
    }
    (args.outdir / "receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
