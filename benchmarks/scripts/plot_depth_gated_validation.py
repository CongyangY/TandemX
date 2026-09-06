"""Render six-panel evidence for frozen depth-gated held-out validation."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, read_table, write_table
from benchmarks.scripts.archive_classifier_validation import (
    confusion,
    validate_paired_classifier_rows,
)


BASELINE = "single_k21"
SELECTED = "depth_gated_blend_v3"
LABELS = {BASELINE: "Single k=21", SELECTED: "Depth-gated v3"}
COLORS = {BASELINE: "#6b6b6b", SELECTED: "#4477aa"}
DELTA_COLORS = {
    "sensitivity": "#4477aa",
    "false_positive_rate": "#d97732",
    "precision": "#aa4465",
}


def _portable_input_path(path: Path, outdir: Path) -> str:
    """Return an input path that remains valid after moving the evidence tree."""
    return Path(os.path.relpath(path.resolve(), start=outdir.resolve())).as_posix()


def _same_confusion(observed: dict, expected: dict) -> None:
    for name in ("available", "unavailable", "TP", "FN", "FP", "TN"):
        if observed[name] != expected.get(name):
            raise ValueError(f"Classifier confusion differs: {name}")
    for name in ("sensitivity", "false_positive_rate", "precision"):
        if not math.isclose(
            float(observed[name]), float(expected.get(name)), rel_tol=0, abs_tol=1e-15
        ):
            raise ValueError(f"Classifier metric differs: {name}")


def load_inputs(
    development: Path, baseline: Path, heldout: Path
) -> tuple[list[dict], list[dict], dict, dict]:
    development_validation = json.loads((development / "validation.json").read_text())
    heldout_validation = json.loads((heldout / "validation.json").read_text())
    heldout_environment = json.loads((heldout / "environment.json").read_text())
    baseline_validation = json.loads((baseline / "validation.json").read_text())
    rows = read_table(heldout / "comparison_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
    })
    localizations = read_table(baseline / "localization_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "assembly_fraction",
        "family_id", "true_assembly_bp", "predicted_assembly_bp", "base_recall",
        "base_precision", "fragments", "truth_fragments",
    })
    if (
        development_validation.get("complete") is not True
        or development_validation.get("acceptance", {}).get("passed") is not True
        or development_validation.get("selected_candidate") != SELECTED
        or heldout_validation.get("complete") is not True
        or heldout_validation.get("acceptance", {}).get("passed") is not True
        or heldout_validation.get("selected_candidate") != SELECTED
        or heldout_environment.get("classifier_model") != SELECTED
        or heldout_environment.get("no_heldout_fit_or_selection") is not True
        or baseline_validation.get("complete") is not True
        or baseline_validation.get("localization_family_rows") != len(localizations)
        or digest_file(heldout / "comparison_metrics.tsv")
            != heldout_validation.get("comparison_metrics_sha256")
    ):
        raise ValueError("Classifier figure requires the complete passed v3 held-out chain")
    grouped = validate_paired_classifier_rows(
        rows, BASELINE, SELECTED,
        int(heldout_validation.get("paired_family_conditions", 0)),
    )
    _same_confusion(confusion(grouped[BASELINE]), heldout_validation["baseline_confusion"])
    _same_confusion(confusion(grouped[SELECTED]), heldout_validation["selected_confusion"])
    seeds = {str(seed) for seed in heldout_environment.get("heldout_seeds", [])}
    localization_keys = {
        (
            row["seed"], row["unit_substitution_rate"], row["array_fragments"],
            row["assembly_fraction"], row["family_id"],
        )
        for row in localizations
    }
    if (
        not seeds
        or {row["seed"] for row in rows} != seeds
        or {row["seed"] for row in localizations} != seeds
        or len(localization_keys) != len(localizations)
    ):
        raise ValueError("Held-out localization and classifier seeds or rows differ")
    return rows, localizations, development_validation, heldout_validation


def _source(result: list[dict], panel: str, metric: str, value: float, **context) -> None:
    row = {
        "panel": panel,
        "metric": metric,
        "dataset": "NA",
        "method": "NA",
        "seed": "NA",
        "coverage": "NA",
        "unit_substitution_rate": "NA",
        "array_fragments": "NA",
        "value": value,
        "numerator": "NA",
        "denominator": "NA",
    }
    row.update(context)
    result.append(row)


def _subset_metrics(rows: list[dict], method: str, **filters: str) -> dict:
    subset = [
        row for row in rows
        if row["method"] == method
        and all(row[name] == value for name, value in filters.items())
    ]
    if not subset:
        raise ValueError(f"No classifier rows for {filters}")
    return confusion(subset)


def plot(development: Path, baseline: Path, heldout: Path, outdir: Path) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    rows, localizations, development_validation, heldout_validation = load_inputs(
        development, baseline, heldout
    )
    outdir.mkdir(parents=True)
    panel_source: list[dict] = []
    methods = (BASELINE, SELECTED)
    overall = {method: _subset_metrics(rows, method) for method in methods}

    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "DejaVu Sans",
        "font.size": 8.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.8))
    fig.subplots_adjust(left=.06, right=.985, top=.91, bottom=.10, wspace=.25, hspace=.34)
    fig.suptitle(
        "Fresh held-out validation of the frozen depth-gated abundance classifier",
        x=.06, ha="left", fontsize=13, fontweight="bold",
    )

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A   Frozen evaluation sequence", loc="left", fontweight="bold")
    design = [
        ("Failure audit", "CONSUMED", "5701-5703; added FP only at nominal 1x", "#d97732"),
        ("v3 development", "PASSED", "single k at estimated depth <2; blend otherwise", "#4477aa"),
        ("Code + rule", "LOCKED", "62892a6; Ubuntu + macOS CI", "#4477aa"),
        ("Fresh held-out", "PASSED", "5801-5803; no fitting or selection", "#228833"),
    ]
    for index, (name, state, detail, color) in enumerate(design):
        y = .85 - index * .21
        ax.add_patch(plt.Rectangle((.02, y-.065), .955, .135, transform=ax.transAxes,
                                   facecolor="#f7f8fa", edgecolor="#d9dde2", lw=.8))
        ax.add_patch(plt.Rectangle((.02, y-.065), .012, .135, transform=ax.transAxes,
                                   facecolor=color, edgecolor=color))
        ax.text(.055, y+.018, name, transform=ax.transAxes, fontweight="bold", va="center")
        ax.text(.055, y-.027, detail, transform=ax.transAxes, fontsize=7.25,
                va="center", color="#555555")
        ax.text(.93, y, state, transform=ax.transAxes, ha="right", va="center",
                color=color, fontweight="bold", fontsize=7.6)
        _source(panel_source, "A", name.lower().replace(" ", "_"),
                1.0 if state in {"PASSED", "LOCKED"} else 0.0, dataset=state.lower())
    low = int(heldout_validation["low_depth_rows"])
    pairs = int(heldout_validation["paired_family_conditions"])
    ax.text(.02, .02, f"{pairs:,} paired conditions; {low:,} below depth 2; known catalogue",
            transform=ax.transAxes, fontsize=7.2, color="#555555")
    _source(panel_source, "A", "paired_family_conditions", float(pairs), denominator=pairs)
    _source(panel_source, "A", "low_depth_rows", float(low), denominator=pairs)

    metric_names = ("sensitivity", "false_positive_rate", "precision")
    ax = axes[0, 1]
    x = np.arange(3)
    width = .34
    delta_sets = (
        ("Six-genome development", development_validation["acceptance"], "#88aadd"),
        ("Fresh held-out", heldout_validation["acceptance"], "#228833"),
    )
    for offset, (dataset, acceptance, color) in zip((-.17, .17), delta_sets):
        values = [acceptance[f"full_{name}_delta"] for name in metric_names]
        bars = ax.bar(x + offset, values, width, color=color, label=dataset)
        ax.bar_label(bars, [f"{value:+.4f}" for value in values], padding=2, fontsize=6.8)
        for metric, value in zip(metric_names, values):
            _source(panel_source, "B", f"{metric}_delta", value, dataset=dataset)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylabel("Delta versus single k=21")
    ax.set_ylim(-.01, .075)
    ax.set_title("B   Development and held-out deltas", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.2)
    ax.grid(axis="y", alpha=.15)

    ax = axes[0, 2]
    outcomes = ("TP", "FN", "FP", "TN")
    x = np.arange(len(outcomes))
    for offset, method in zip((-.17, .17), methods):
        values = [overall[method][outcome] for outcome in outcomes]
        bars = ax.bar(x + offset, values, width, color=COLORS[method], label=LABELS[method])
        ax.bar_label(bars, [f"{value:,}" for value in values], padding=2, fontsize=7)
        for outcome, value in zip(outcomes, values):
            _source(panel_source, "C", outcome, float(value), method=method)
    ax.set_xticks(x, outcomes)
    ax.set_ylabel("Family-condition count")
    ax.set_ylim(0, 1080)
    ax.set_title("C   Held-out confusion counts", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.2)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 0]
    seeds = sorted({row["seed"] for row in rows}, key=int)
    x = np.arange(len(seeds))
    width = .23
    for offset, metric in zip((-width, 0, width), metric_names):
        values = []
        for seed in seeds:
            baseline_metric = _subset_metrics(rows, BASELINE, seed=seed)[metric]
            selected_metric = _subset_metrics(rows, SELECTED, seed=seed)[metric]
            value = selected_metric - baseline_metric
            values.append(value)
            _source(panel_source, "D", f"{metric}_delta", value, seed=seed)
        bars = ax.bar(x + offset, values, width, color=DELTA_COLORS[metric],
                      label=metric.replace("_", " ").title())
        ax.bar_label(bars, [f"{value:+.4f}" for value in values], padding=2,
                     fontsize=6.4, rotation=90)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(x, seeds)
    ax.set_ylabel("Depth-gated v3 - baseline")
    ax.set_ylim(-.008, .098)
    ax.set_title("D   Metric deltas by held-out seed", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 1]
    coverages = sorted({row["coverage"] for row in rows}, key=float)
    x = np.arange(len(coverages))
    for offset, metric in zip((-width, 0, width), metric_names):
        values = []
        for coverage in coverages:
            baseline_metric = _subset_metrics(rows, BASELINE, coverage=coverage)[metric]
            selected_metric = _subset_metrics(rows, SELECTED, coverage=coverage)[metric]
            value = selected_metric - baseline_metric
            values.append(value)
            _source(panel_source, "E", f"{metric}_delta", value, coverage=coverage)
        bars = ax.bar(x + offset, values, width, color=DELTA_COLORS[metric],
                      label=metric.replace("_", " ").title())
        ax.bar_label(bars, [f"{value:+.4f}" for value in values], padding=2,
                     fontsize=6.5, rotation=90)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(x, [f"{float(value):g}x" for value in coverages])
    ax.set_ylabel("Depth-gated v3 - baseline")
    ax.set_ylim(-.008, .115)
    ax.set_title("E   Metric deltas by nominal coverage", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 2]
    full = [row for row in localizations if float(row["assembly_fraction"]) == 1]
    positive = [row for row in localizations if float(row["assembly_fraction"]) > 0]
    absent = [row for row in localizations if float(row["assembly_fraction"]) == 0]
    divergences = sorted({row["unit_substitution_rate"] for row in full}, key=float)
    fragments = sorted({row["array_fragments"] for row in full}, key=int)
    x = np.arange(len(divergences))
    width = .34
    for offset, fragment_count, color in zip((-.17, .17), fragments, ("#4477aa", "#d97732")):
        values = []
        for divergence in divergences:
            subset = [
                row for row in full
                if row["unit_substitution_rate"] == divergence
                and row["array_fragments"] == fragment_count
            ]
            recalls = [float(row["base_recall"]) for row in subset]
            value = sum(recalls) / len(recalls)
            values.append(value)
            _source(panel_source, "F", "full_assembly_mean_base_recall", value,
                    unit_substitution_rate=divergence, array_fragments=fragment_count,
                    numerator=sum(recalls), denominator=len(recalls))
        bars = ax.bar(x + offset, values, width, color=color,
                      label=f"{fragment_count} segment{'s' if fragment_count != '1' else ''}")
        ax.bar_label(bars, [f"{value:.3f}" for value in values], padding=2, fontsize=6.8)
    ax.axhline(.95, color="#555555", lw=1, linestyle=":",
               label="Aggregate target (reference)")
    ax.set_xticks(x, [f"{float(value)*100:g}%" for value in divergences])
    ax.set_ylim(.90, 1.012)
    ax.set_ylabel("Mean base recall")
    ax.set_title("F   Full-assembly localization strata", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, loc="lower left")
    ax.grid(axis="y", alpha=.15)
    precisions = [float(row["base_precision"]) for row in positive]
    mean_precision = sum(precisions) / len(precisions)
    absent_fp = sum(float(row["predicted_assembly_bp"]) > 0 for row in absent)
    _source(panel_source, "F", "positive_assembly_mean_base_precision", mean_precision,
            numerator=sum(precisions), denominator=len(precisions))
    _source(panel_source, "F", "absent_family_false_positive_rows", float(absent_fp),
            numerator=absent_fp, denominator=len(absent))
    ax.text(.99, .02, f"Precision {mean_precision:.5f}; absent FP {absent_fp}/{len(absent)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.9,
            color="#555555")

    write_table(outdir / "panel_source.tsv", panel_source, list(panel_source[0]))
    stem = outdir / "depth_gated_validation"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
    provenance = {
        "complete": True,
        "panels": ["A", "B", "C", "D", "E", "F"],
        "selected_candidate": SELECTED,
        "development_gate_passed": True,
        "heldout_gate_passed": True,
        "no_heldout_fit_or_selection": True,
        "source_rows": len(panel_source),
        "inputs": {
            "development_validation": {
                "path": _portable_input_path(development / "validation.json", outdir),
                "sha256": digest_file(development / "validation.json"),
            },
            "heldout_validation": {
                "path": _portable_input_path(heldout / "validation.json", outdir),
                "sha256": digest_file(heldout / "validation.json"),
            },
            "heldout_comparison_metrics": {
                "path": _portable_input_path(
                    heldout / "comparison_metrics.tsv", outdir
                ),
                "sha256": digest_file(heldout / "comparison_metrics.tsv"),
            },
            "baseline_localization_metrics": {
                "path": _portable_input_path(
                    baseline / "localization_metrics.tsv", outdir
                ),
                "sha256": digest_file(baseline / "localization_metrics.tsv"),
            },
        },
        "outputs": {
            path.name: {"sha256": digest_file(path), "bytes": path.stat().st_size}
            for path in (
                outdir / "panel_source.tsv",
                stem.with_suffix(".svg"),
                stem.with_suffix(".pdf"),
                stem.with_suffix(".png"),
            )
        },
        "warning": (
            "known_catalogue_IID_substitution_simulation;frozen_heldout_gate_passed;"
            "no_heldout_refitting;not_biological_validation"
        ),
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--heldout", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.development, args.baseline, args.heldout, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
