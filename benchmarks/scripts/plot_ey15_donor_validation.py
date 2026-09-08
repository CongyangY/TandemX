#!/usr/bin/env python3
"""Plot the six-panel Ey15 donor-matched validation figure."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from benchmarks.challenge.schema import digest_file


BLUE = "#2A6FBB"
BLUE_LIGHT = "#CFE0F2"
ORANGE = "#D97706"
ORANGE_LIGHT = "#F5D8B5"
CHARCOAL = "#30343B"
GREY = "#9299A1"
LIGHT_GREY = "#E4E7EA"
GRID = "#D9DDE2"
PANELS = "ABCDEF"
SYMLOG_LINTHRESH = 100.0
SYMLOG_LOWER_LIMIT = -10.0
SOURCE_FIELDS = (
    "panel",
    "record_type",
    "family_id",
    "group",
    "metric",
    "threshold_bp",
    "normalization",
    "state",
    "x_value",
    "y_value",
    "numerator",
    "denominator",
    "value",
    "warning",
)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric_rows(summary: dict[str, Any], normalization: str) -> list[dict[str, Any]]:
    specs = (
        ("Sensitivity", "sensitivity", "sensitivity_wilson95"),
        ("Precision", "precision", "precision_wilson95"),
        ("False-positive rate", "false_positive_rate", "false_positive_rate_wilson95"),
    )
    rows = []
    for label, value_key, interval_key in specs:
        interval = summary[interval_key]
        rows.append(
            {
                "normalization": normalization,
                "metric": label,
                "value": summary[value_key],
                "lower": interval[0] if interval is not None else None,
                "upper": interval[1] if interval is not None else None,
            }
        )
    return rows


def annotation_class_recall(
    rows: list[dict[str, str]], annotation_class: str
) -> float:
    matches = [row for row in rows if row["annotation_class"] == annotation_class]
    if len(matches) != 1:
        raise ValueError(
            f"expected one annotation-class row for {annotation_class!r}; "
            f"found {len(matches)}"
        )
    return float(matches[0]["annotation_bp_recall"])


def family_style(state: str) -> dict[str, Any]:
    if state == "reference_collapse":
        return {"facecolors": BLUE, "edgecolors": CHARCOAL, "marker": "o"}
    return {"facecolors": "white", "edgecolors": CHARCOAL, "marker": "s"}


def scatter_by_state(
    axis: plt.Axes,
    rows: list[dict[str, str]],
    x_field: str,
    y_field: str,
) -> None:
    for state, label in (
        ("reference_collapse", "Reference collapse"),
        ("reference_retained", "Reference retained"),
    ):
        selected = [row for row in rows if row["reference_state"] == state]
        if not selected:
            continue
        axis.scatter(
            [float(row[x_field]) for row in selected],
            [float(row[y_field]) for row in selected],
            s=38,
            linewidths=0.9,
            label=f"{label} (n={len(selected)})",
            zorder=3,
            **family_style(state),
        )


def line_extent(rows: list[dict[str, str]], fields: tuple[str, ...]) -> float:
    return max(1.0, max(float(row[field]) for row in rows for field in fields))


def style_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(CHARCOAL)
    axis.spines["bottom"].set_color(CHARCOAL)
    axis.tick_params(colors=CHARCOAL, labelsize=8)
    axis.grid(True, color=GRID, linewidth=0.55, alpha=0.7, zorder=0)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.13,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=12,
        fontweight="bold",
        color=CHARCOAL,
        va="top",
    )


def add_blossom(figure: plt.Figure) -> None:
    centre_x, centre_y = 0.977, 0.982
    for index in range(5):
        angle = 2 * math.pi * index / 5
        figure.add_artist(
            Circle(
                (centre_x + 0.009 * math.cos(angle), centre_y + 0.007 * math.sin(angle)),
                0.004,
                transform=figure.transFigure,
                facecolor=BLUE_LIGHT if index % 2 == 0 else ORANGE_LIGHT,
                edgecolor=CHARCOAL,
                linewidth=0.35,
            )
        )
    figure.add_artist(
        Circle(
            (centre_x, centre_y),
            0.003,
            transform=figure.transFigure,
            facecolor=CHARCOAL,
            edgecolor="none",
        )
    )


def plot(source: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {outdir}")
    paths = {
        "primary_rows": source / "evaluation_primary_total_bases/family_metrics.tsv",
        "primary_summary": source / "evaluation_primary_total_bases/summary.json",
        "primary_verification": source
        / "evaluation_primary_total_bases/independent_verification.json",
        "depth107_rows": source
        / "evaluation_depth107_sensitivity/family_metrics.tsv",
        "depth107_summary": source / "evaluation_depth107_sensitivity/summary.json",
        "depth107_verification": source
        / "evaluation_depth107_sensitivity/independent_verification.json",
        "annotation_rows": source
        / "reference_annotation_context_posthoc_v2/reference_annotation_rows.tsv",
        "annotation_summary": source
        / "reference_annotation_context_posthoc_v2/summary.json",
        "alignment_rows": source
        / "old_new_alignment_context/family_alignment_context.tsv",
        "alignment_summary": source / "old_new_alignment_context/summary.json",
        "alignment_verification": source
        / "old_new_alignment_context/independent_verification.json",
        "author_class_summary": source
        / "author_annotation_audit/annotation_class_summary.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"figure inputs are incomplete: {missing}")
    for key in ("primary_verification", "depth107_verification", "alignment_verification"):
        verification = load_json(paths[key])
        if verification.get("verification_passed") is not True:
            raise ValueError(f"independent verification failed: {paths[key]}")

    primary_all = load_tsv(paths["primary_rows"])
    primary = [row for row in primary_all if row["eligibility"] == "eligible"]
    depth107_all = load_tsv(paths["depth107_rows"])
    depth107 = [row for row in depth107_all if row["eligibility"] == "eligible"]
    primary_summary = load_json(paths["primary_summary"])
    depth107_summary = load_json(paths["depth107_summary"])
    annotation_rows = load_tsv(paths["annotation_rows"])
    annotation_summary = load_json(paths["annotation_summary"])
    alignment_rows = load_tsv(paths["alignment_rows"])
    alignment_summary = load_json(paths["alignment_summary"])
    author_class_rows = load_tsv(paths["author_class_summary"])
    centromere_annotation_recall = annotation_class_recall(
        author_class_rows, "centromere"
    )
    if len(primary) != 19 or {row["family_id"] for row in primary} != {
        row["family_id"] for row in depth107
    }:
        raise ValueError("expected the same frozen 19-family primary denominator")
    if {row["family_id"] for row in primary} != {
        row["family_id"] for row in annotation_rows
    } or {row["family_id"] for row in primary} != {
        row["family_id"] for row in alignment_rows
    }:
        raise ValueError("orthogonal context does not cover the complete primary denominator")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "text.color": CHARCOAL,
            "axes.labelcolor": CHARCOAL,
            "axes.titlecolor": CHARCOAL,
        }
    )
    figure, axes = plt.subplots(3, 2, figsize=(8.3, 10.4), constrained_layout=False)
    figure.subplots_adjust(left=0.10, right=0.97, top=0.925, bottom=0.075, hspace=0.47, wspace=0.31)
    figure.suptitle("Ey15-2 donor-matched assembly-collapse validation", x=0.10, ha="left", fontsize=15, fontweight="bold")
    figure.text(
        0.10,
        0.945,
        "Complete HiFi reads; CLR-Canu old assembly; HiFi-Hifiasm reference proxy; frozen 15-kb primary denominator (n=19)",
        fontsize=8.5,
        color=CHARCOAL,
        ha="left",
    )
    add_blossom(figure)
    source_rows: list[dict[str, Any]] = []

    # A: reference contrast.
    axis = axes[0, 0]
    scatter_by_state(axis, primary, "new_assembly_bp", "old_assembly_bp")
    maximum = line_extent(primary, ("new_assembly_bp", "old_assembly_bp")) * 1.25
    axis.plot([0, maximum], [0, maximum], color=GREY, linewidth=1.0, linestyle="--", label="Equality")
    axis.plot([0, maximum], [0, 0.6 * maximum], color=ORANGE, linewidth=1.1, linestyle=":", label="Collapse threshold")
    axis.set_xscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_yscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_xlim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_ylim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_xlabel("New-assembly localized bp")
    axis.set_ylabel("Old-assembly localized bp")
    axis.set_title("Assembly reference states")
    axis.legend(frameon=False, fontsize=6.7, loc="upper left")
    style_axis(axis)
    panel_label(axis, "A")
    for row in primary:
        source_rows.append(
            {
                "panel": "A",
                "record_type": "family_point",
                "family_id": row["family_id"],
                "state": row["reference_state"],
                "x_value": row["new_assembly_bp"],
                "y_value": row["old_assembly_bp"],
                "denominator": 19,
            }
        )

    # B: read prediction contrast.
    axis = axes[0, 1]
    scatter_by_state(axis, primary, "read_estimated_bp", "old_assembly_bp")
    maximum = line_extent(primary, ("read_estimated_bp", "old_assembly_bp")) * 1.25
    axis.plot([0, maximum], [0, maximum], color=GREY, linewidth=1.0, linestyle="--")
    axis.plot([0, maximum], [0, 0.6 * maximum], color=ORANGE, linewidth=1.1, linestyle=":")
    axis.set_xscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_yscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_xlim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_ylim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_xlabel("Read-estimated repeat bp")
    axis.set_ylabel("Old-assembly localized bp")
    axis.set_title("Read-based collapse predictions")
    style_axis(axis)
    panel_label(axis, "B")
    for row in primary:
        source_rows.append(
            {
                "panel": "B",
                "record_type": "family_point",
                "family_id": row["family_id"],
                "state": row["prediction_state"],
                "x_value": row["read_estimated_bp"],
                "y_value": row["old_assembly_bp"],
                "denominator": 19,
            }
        )

    # C: amount agreement.
    axis = axes[1, 0]
    scatter_by_state(axis, primary, "observed_gain_bp", "predicted_missing_bp")
    maximum = line_extent(primary, ("observed_gain_bp", "predicted_missing_bp")) * 1.25
    axis.plot([0, maximum], [0, maximum], color=CHARCOAL, linewidth=1.0, linestyle="--")
    axis.set_xscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_yscale("symlog", linthresh=SYMLOG_LINTHRESH)
    axis.set_xlim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_ylim(SYMLOG_LOWER_LIMIT, maximum)
    axis.set_xlabel("Observed old-to-new gain (bp)")
    axis.set_ylabel("Predicted missing bp")
    axis.set_title("Missing-sequence magnitude")
    axis.text(
        0.03,
        0.96,
        f"Pearson r = {primary_summary['missing_bp_pearson']:.3f}\nSpearman ρ = {primary_summary['missing_bp_spearman']:.3f}",
        transform=axis.transAxes,
        va="top",
        fontsize=7.4,
    )
    style_axis(axis)
    panel_label(axis, "C")
    for row in primary:
        source_rows.append(
            {
                "panel": "C",
                "record_type": "family_point",
                "family_id": row["family_id"],
                "state": row["reference_state"],
                "x_value": row["observed_gain_bp"],
                "y_value": row["predicted_missing_bp"],
                "denominator": 19,
            }
        )

    # D: frozen size-threshold sensitivity.
    axis = axes[1, 1]
    threshold_rows = primary_summary["min_new_bp_sensitivity"]
    thresholds = [row["min_new_bp"] for row in threshold_rows]
    bottoms = [0] * len(thresholds)
    outcomes = (
        ("TP", BLUE, ""),
        ("TN", GREY, ""),
        ("FN", ORANGE_LIGHT, "///"),
        ("FP", ORANGE, "\\\\"),
    )
    for outcome, color, hatch in outcomes:
        values = [row[outcome] for row in threshold_rows]
        axis.bar(
            range(len(thresholds)),
            values,
            bottom=bottoms,
            color=color,
            edgecolor=CHARCOAL,
            linewidth=0.6,
            hatch=hatch,
            label=outcome,
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        for row, value in zip(threshold_rows, values):
            source_rows.append(
                {
                    "panel": "D",
                    "record_type": "confusion_count",
                    "metric": outcome,
                    "threshold_bp": row["min_new_bp"],
                    "numerator": value,
                    "denominator": row["families"],
                    "value": value,
                    "warning": "predeclared_threshold_sensitivity",
                }
            )
    for index, row in enumerate(threshold_rows):
        axis.text(index, row["families"] + 1.0, f"n={row['families']}", ha="center", fontsize=7.2)
    axis.set_xticks(range(len(thresholds)), [f"{value / 1000:g} kb" for value in thresholds])
    axis.set_ylabel("Families")
    axis.set_title("Predeclared denominator sensitivity")
    axis.legend(frameon=False, fontsize=7, ncol=4, loc="upper right")
    style_axis(axis)
    panel_label(axis, "D")

    # E: normalization sensitivity with Wilson intervals.
    axis = axes[2, 0]
    normalizations = (
        ("Total bases / 143.12 Mb", metric_rows(primary_summary, "total_bases"), BLUE, "o", -0.09),
        ("Explicit 107×", metric_rows(depth107_summary, "depth107"), ORANGE, "s", 0.09),
    )
    metric_labels = [row["metric"] for row in normalizations[0][1]]
    for label, rows, color, marker, offset in normalizations:
        xvalues = [index + offset for index in range(len(rows))]
        values = [row["value"] for row in rows]
        lower = [value - row["lower"] for value, row in zip(values, rows)]
        upper = [row["upper"] - value for value, row in zip(values, rows)]
        axis.errorbar(
            xvalues,
            values,
            yerr=[lower, upper],
            fmt=marker,
            color=color,
            markerfacecolor="white" if marker == "s" else color,
            markeredgecolor=CHARCOAL,
            capsize=3,
            linewidth=1.1,
            markersize=5,
            label=label,
        )
        for row in rows:
            source_rows.append(
                {
                    "panel": "E",
                    "record_type": "metric_with_wilson95",
                    "metric": row["metric"],
                    "normalization": row["normalization"],
                    "value": row["value"],
                    "x_value": row["lower"],
                    "y_value": row["upper"],
                    "denominator": 19,
                    "warning": "same_reads_not_biological_replicates",
                }
            )
    axis.set_xticks(range(len(metric_labels)), metric_labels, rotation=12, ha="right")
    axis.set_ylim(-0.05, 1.08)
    axis.set_ylabel("Estimate and Wilson 95% interval")
    axis.set_title("Depth-normalization sensitivity")
    axis.legend(frameon=False, fontsize=7, loc="lower left")
    style_axis(axis)
    panel_label(axis, "E")

    # F: orthogonal context.
    axis = axes[2, 1]
    context = {
        "Author annotation\noverlap": {
            "reference_collapse": annotation_summary["contingency"]["collapse_annotated"]
            / annotation_summary["reference_collapse_families"],
            "other_reference_state": annotation_summary["contingency"]["other_annotated"]
            / annotation_summary["other_reference_state_families"],
        },
        "Same-chromosome primary\nalignment coverage": {
            state: alignment_summary["groups"][state][
                "same_chromosome_primary_alignment_fraction_median"
            ]
            for state in ("reference_collapse", "other_reference_state")
        },
    }
    xvalues = range(len(context))
    width = 0.34
    for offset, state, label, color, hatch in (
        (-width / 2, "reference_collapse", "Reference collapse (n=8)", BLUE, ""),
        (width / 2, "other_reference_state", "Other state (n=11)", LIGHT_GREY, "//"),
    ):
        values = [row[state] for row in context.values()]
        axis.bar(
            [x + offset for x in xvalues],
            values,
            width,
            color=color,
            edgecolor=CHARCOAL,
            linewidth=0.7,
            hatch=hatch,
            label=label,
        )
        for metric, value in zip(context, values):
            source_rows.append(
                {
                    "panel": "F",
                    "record_type": "orthogonal_context_fraction",
                    "group": state,
                    "metric": metric.replace("\n", " "),
                    "value": value,
                    "denominator": 8 if state == "reference_collapse" else 11,
                    "warning": (
                        "annotation_is_posthoc;alignment_is_context_not_independent_copy_truth"
                    ),
                }
            )
    axis.set_xticks(list(xvalues), list(context))
    axis.set_ylim(0, 1.08)
    axis.set_ylabel("Fraction or median coverage")
    axis.set_title("Orthogonal assembly context")
    axis.text(
        0.02,
        0.97,
        f"Annotation Fisher P={annotation_summary['fisher_exact_greater_p']:.4g} (post hoc)\nAuthor centromere recall={centromere_annotation_recall:.3%}",
        transform=axis.transAxes,
        va="top",
        fontsize=7.0,
    )
    axis.legend(frameon=False, fontsize=6.7, loc="lower right")
    style_axis(axis)
    panel_label(axis, "F")

    figure.text(
        0.10,
        0.025,
        "High-quality reference proxy, not absolute truth. New assembly shares HiFi evidence; same reads underlie both normalization choices. Full family fates and failures remain in source tables.",
        fontsize=7.2,
        color=CHARCOAL,
        ha="left",
    )
    outdir.mkdir(parents=True)
    svg = outdir / "ey15_donor_validation.svg"
    pdf = outdir / "ey15_donor_validation.pdf"
    png = outdir / "ey15_donor_validation.png"
    figure.savefig(svg, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=240, bbox_inches="tight")
    plt.close(figure)
    panel_source = outdir / "panel_source.tsv"
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
    svg_text = svg.read_text(encoding="utf-8")
    receipt = {
        "schema_version": 1,
        "complete": True,
        "panel_count": 6,
        "primary_family_count": len(primary),
        "panel_source_rows": len(source_rows),
        "svg_text_element_count": svg_text.count("<text"),
        "svg_raster_image_element_count": svg_text.count("<image"),
        "source_files": {
            name: {"path": str(path), "sha256": digest_file(path)}
            for name, path in paths.items()
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": digest_file(path)}
            for path in (svg, pdf, png, panel_source)
        },
        "warning": (
            "single_donor_reference_proxy;normalization_sensitivity_not_replication;"
            "annotation_association_posthoc;alignment_context_not_copy_truth"
        ),
    }
    (outdir / "figure_provenance.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.source, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
