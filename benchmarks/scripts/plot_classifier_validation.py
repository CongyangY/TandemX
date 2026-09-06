"""Render six-panel evidence for the frozen abundance-classifier validation."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
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
SELECTED = "blend_a0.5_t0.5"
LABELS = {BASELINE: "Single k=21", SELECTED: "Frozen blend"}
COLORS = {BASELINE: "#6b6b6b", SELECTED: "#4477aa"}
DELTA_COLORS = {
    "sensitivity": "#4477aa",
    "false_positive_rate": "#d97732",
    "precision": "#aa4465",
}


def _same_confusion(observed: dict, expected: dict) -> None:
    for name in ("available", "unavailable", "TP", "FN", "FP", "TN"):
        if observed[name] != expected.get(name):
            raise ValueError(f"Classifier confusion differs: {name}")
    for name in ("sensitivity", "false_positive_rate", "precision"):
        if not math.isclose(float(observed[name]), float(expected.get(name)), abs_tol=1e-15):
            raise ValueError(f"Classifier metric differs: {name}")


def load_inputs(
    development: Path,
    baseline: Path,
    heldout: Path,
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
        or heldout_validation.get("acceptance", {}).get("passed") is not False
        or heldout_validation.get("selected_candidate") != SELECTED
        or heldout_environment.get("classifier_model") != SELECTED
        or baseline_validation.get("complete") is not True
        or baseline_validation.get("localization_family_rows") != len(localizations)
        or digest_file(heldout / "comparison_metrics.tsv")
            != heldout_validation.get("comparison_metrics_sha256")
    ):
        raise ValueError("Classifier figure requires the complete failed held-out chain")
    grouped = validate_paired_classifier_rows(
        rows, BASELINE, SELECTED,
        int(heldout_validation.get("paired_family_conditions", 0)),
    )
    _same_confusion(confusion(grouped[BASELINE]), heldout_validation["baseline_confusion"])
    _same_confusion(confusion(grouped[SELECTED]), heldout_validation["selected_confusion"])
    seeds = {str(seed) for seed in heldout_environment.get("heldout_seeds", [])}
    localization_keys = {
        (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
         row["assembly_fraction"], row["family_id"])
        for row in localizations
    }
    if (
        not seeds or {row["seed"] for row in rows} != seeds
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
        if row["method"] == method and all(row[name] == value for name, value in filters.items())
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
    fig.subplots_adjust(left=.06, right=.985, top=.91, bottom=.11, wspace=.24, hspace=.32)
    fig.suptitle(
        "Development and held-out evaluation of the frozen abundance classifier",
        x=.06, ha="left", fontsize=13, fontweight="bold",
    )

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A   Frozen evaluation design", loc="left", fontweight="bold")
    design = [
        ("Development v1", "FAILED", "pooled selection unstable"),
        ("Development v2", "PASSED", "seed-robust alpha=0.5, t=0.5"),
        ("Code + model", "FROZEN", "commit 300e48d; CI passed"),
        ("Held-out", "FAILED", "seeds 5701-5703; no refitting"),
    ]
    y_values = (.86, .65, .44, .23)
    for (name, state, detail), y in zip(design, y_values):
        face = "#eef4f8" if state in {"PASSED", "FROZEN"} else "#fff1e8"
        edge = "#4477aa" if state in {"PASSED", "FROZEN"} else "#d97732"
        ax.text(.03, y, name, transform=ax.transAxes, fontweight="bold", va="center")
        ax.text(.46, y, state, transform=ax.transAxes, ha="center", va="center",
                color=edge, fontweight="bold",
                bbox={"boxstyle": "round,pad=.28", "fc": face, "ec": edge, "lw": 1})
        ax.text(.61, y, detail, transform=ax.transAxes, va="center", fontsize=7.6)
        _source(panel_source, "A", name.lower().replace(" ", "_"),
                1.0 if state in {"PASSED", "FROZEN"} else 0.0,
                dataset=state.lower())
    ax.annotate("", xy=(.34, .13), xytext=(.34, .92), xycoords="axes fraction",
                arrowprops={"arrowstyle": "->", "color": "#9a9a9a", "lw": 1})
    ax.text(.03, .05, "2,430 paired family conditions; known catalogue; synthetic IID substitutions",
            transform=ax.transAxes, color="#555555", fontsize=7.4)

    ax = axes[0, 1]
    metric_names = ("sensitivity", "false_positive_rate", "precision")
    x = np.arange(3)
    width = .34
    dev_accept = development_validation["acceptance"]
    held_accept = heldout_validation["acceptance"]
    deltas = {
        "Development v2": [
            dev_accept["full_sensitivity_delta"],
            dev_accept["full_false_positive_rate_delta"],
            dev_accept["full_precision_delta"],
        ],
        "Held-out": [
            held_accept["full_sensitivity_delta"],
            held_accept["full_false_positive_rate_delta"],
            held_accept["full_precision_delta"],
        ],
    }
    for offset, (dataset, values), color in zip(
        (-width / 2, width / 2), deltas.items(), ("#4477aa", "#d97732")
    ):
        bars = ax.bar(x + offset, values, width, color=color, label=dataset)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, value + (.004 if value >= 0 else -.006),
                    f"{value:+.3f}", ha="center", va="bottom" if value >= 0 else "top",
                    fontsize=7)
        for metric, value in zip(metric_names, values):
            _source(panel_source, "B", f"{metric}_delta", value, dataset=dataset)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylabel("Delta versus single k=21")
    ax.set_ylim(-.035, .135)
    ax.set_title("B   Development-to-validation transfer", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", alpha=.15)

    ax = axes[0, 2]
    x = np.arange(3)
    for offset, method in zip((-width / 2, width / 2), methods):
        values = [overall[method][name] for name in metric_names]
        bars = ax.bar(x + offset, values, width, color=COLORS[method], label=LABELS[method])
        ax.bar_label(bars, [f"{value:.3f}" for value in values], padding=2, fontsize=7)
        for metric, value in zip(metric_names, values):
            _source(panel_source, "C", metric, float(value), method=method)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylim(0, 1.12)
    ax.set_title("C   Held-out aggregate metrics", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 0]
    seeds = sorted({row["seed"] for row in rows}, key=int)
    x = np.arange(len(seeds))
    width = .23
    for offset, metric in zip((-width, 0, width), metric_names):
        values = []
        for seed in seeds:
            base_metrics = _subset_metrics(rows, BASELINE, seed=seed)
            selected_metrics = _subset_metrics(rows, SELECTED, seed=seed)
            value = selected_metrics[metric] - base_metrics[metric]
            values.append(value)
            _source(panel_source, "D", f"{metric}_delta", value, seed=seed)
        bars = ax.bar(x + offset, values, width, color=DELTA_COLORS[metric],
                      label=metric.replace("_", " ").title())
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, value + (.004 if value >= 0 else -.005),
                    f"{value:+.3f}", ha="center", va="bottom" if value >= 0 else "top",
                    fontsize=6.5, rotation=90)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(x, seeds)
    ax.set_ylabel("Selected - baseline")
    ax.set_ylim(-.045, .13)
    ax.set_title("D   Metric deltas by held-out seed", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.8, ncol=1, loc="upper left")
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 1]
    coverages = sorted({row["coverage"] for row in rows}, key=float)
    for metric, marker in (("sensitivity", "o"), ("false_positive_rate", "s")):
        values = []
        for coverage in coverages:
            base_metrics = _subset_metrics(rows, BASELINE, coverage=coverage)
            selected_metrics = _subset_metrics(rows, SELECTED, coverage=coverage)
            value = selected_metrics[metric] - base_metrics[metric]
            values.append(value)
            _source(panel_source, "E", f"{metric}_delta", value, coverage=coverage)
        ax.plot(np.arange(len(coverages)), values, marker=marker, lw=2,
                color=DELTA_COLORS[metric], label=metric.replace("_", " ").title())
        for xpos, value in enumerate(values):
            ax.text(xpos, value + .006, f"{value:+.3f}", ha="center", fontsize=7)
    ax.axhline(0, color="#333333", lw=.8)
    ax.set_xticks(np.arange(len(coverages)), [f"{value}x" for value in coverages])
    ax.set_ylabel("Selected - baseline")
    ax.set_ylim(-.015, .12)
    ax.set_title("E   Coverage-specific trade-off", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 2]
    full = [row for row in localizations if float(row["assembly_fraction"]) == 1]
    positive = [row for row in localizations if float(row["assembly_fraction"]) > 0]
    absent = [row for row in localizations if float(row["assembly_fraction"]) == 0]
    divergences = sorted({row["unit_substitution_rate"] for row in full}, key=float)
    fragment_counts = sorted({row["array_fragments"] for row in full}, key=int)
    for fragments, color, marker, linestyle in zip(
        fragment_counts, ("#4477aa", "#d97732"), ("o", "s"), ("-", "--")
    ):
        values = []
        for divergence in divergences:
            subset = [
                row for row in full
                if row["unit_substitution_rate"] == divergence
                and row["array_fragments"] == fragments
            ]
            value = sum(float(row["base_recall"]) for row in subset) / len(subset)
            values.append(value)
            _source(panel_source, "F", "full_assembly_mean_base_recall", value,
                    unit_substitution_rate=divergence, array_fragments=fragments,
                    numerator=sum(float(row["overlap_bp"]) for row in subset),
                    denominator=sum(float(row["true_assembly_bp"]) for row in subset))
        ax.plot(np.arange(len(divergences)), values, marker=marker, linestyle=linestyle,
                color=color, lw=2, label=f"{fragments} segment{'s' if fragments != '1' else ''}")
        for xpos, value in enumerate(values):
            ax.text(xpos, value + .004, f"{value:.3f}", ha="center", fontsize=6.8)
    ax.axhline(.95, color="#555555", lw=1, linestyle=":", label="Predeclared mean gate")
    ax.set_xticks(np.arange(len(divergences)), [f"{float(value)*100:g}%" for value in divergences])
    ax.set_ylim(.90, 1.01)
    ax.set_ylabel("Mean base recall")
    ax.set_title("F   Held-out full-assembly localization", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.grid(axis="y", alpha=.15)
    mean_precision = sum(float(row["base_precision"]) for row in positive) / len(positive)
    absent_fp = sum(float(row["predicted_assembly_bp"]) > 0 for row in absent)
    _source(panel_source, "F", "positive_assembly_mean_base_precision", mean_precision,
            numerator=len(positive), denominator=len(positive))
    _source(panel_source, "F", "absent_family_false_positive_rows", float(absent_fp),
            numerator=absent_fp, denominator=len(absent))
    ax.text(.99, .02, f"Positive precision {mean_precision:.5f}; absent FP {absent_fp}/{len(absent)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#555555")

    fields = list(panel_source[0])
    write_table(outdir / "panel_source.tsv", panel_source, fields)
    stem = outdir / "classifier_validation"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
    provenance = {
        "complete": True,
        "panels": ["A", "B", "C", "D", "E", "F"],
        "selected_candidate": SELECTED,
        "development_gate_passed": True,
        "heldout_gate_passed": False,
        "source_rows": len(panel_source),
        "inputs": {
            "development_validation": {
                "path": str((development / "validation.json").resolve()),
                "sha256": digest_file(development / "validation.json"),
            },
            "heldout_validation": {
                "path": str((heldout / "validation.json").resolve()),
                "sha256": digest_file(heldout / "validation.json"),
            },
            "heldout_comparison_metrics": {
                "path": str((heldout / "comparison_metrics.tsv").resolve()),
                "sha256": digest_file(heldout / "comparison_metrics.tsv"),
            },
            "baseline_localization_metrics": {
                "path": str((baseline / "localization_metrics.tsv").resolve()),
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
            "known_catalogue_IID_substitution_simulation;heldout_classifier_gate_failed;"
            "no_posthoc_refitting;not_biological_validation"
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
