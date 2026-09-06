"""Render six-panel evidence for divergent and interrupted assembly validation."""
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


BASELINE = "single_k21"
FROZEN = "multik_fallback_depth_rule"
METHODS = (BASELINE, FROZEN)
METHOD_LABELS = {BASELINE: "Single k=21", FROZEN: "Frozen multi-k rule"}
COLORS = {BASELINE: "#6b6b6b", FROZEN: "#b23a48"}
MARKERS = {1: "o", 3: "s"}


def confusion(rows: list[dict]) -> dict[str, float | int]:
    counts = Counter(row["outcome"] for row in rows)
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    return {
        "TP": tp, "FN": fn, "FP": fp, "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else float("nan"),
        "false_positive_rate": fp / (fp + tn) if fp + tn else float("nan"),
        "precision": tp / (tp + fp) if tp + fp else float("nan"),
    }


def load_inputs(
    metrics_path: Path,
    validation_path: Path,
    localization_path: Path,
    baseline_validation_path: Path,
) -> tuple[list[dict], list[dict], dict]:
    rows = read_table(metrics_path, {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
    })
    validation = json.loads(validation_path.read_text())
    baseline_validation = json.loads(baseline_validation_path.read_text())
    if (
        validation.get("complete") is not True or validation.get("split") != "heldout"
        or baseline_validation.get("complete") is not True
    ):
        raise ValueError("Require complete paired held-out and baseline evidence")
    selected = [row for row in rows if row["method"] in METHODS]
    grouped = {method: [row for row in selected if row["method"] == method]
               for method in METHODS}
    keys = {
        method: {
            (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
             row["coverage"], row["substitution_rate"], row["assembly_fraction"],
             row["family_id"])
            for row in group
        }
        for method, group in grouped.items()
    }
    if (
        not selected or len(keys[BASELINE]) != len(grouped[BASELINE])
        or keys[BASELINE] != keys[FROZEN]
    ):
        raise ValueError("Baseline and frozen-rule domain-shift rows are not exactly paired")
    for method, group in grouped.items():
        observed = confusion(group)
        expected = validation.get("method_confusion", {}).get(method, {})
        if any(observed[name] != expected.get(name) for name in ("TP", "FN", "FP", "TN")):
            raise ValueError(f"Validation confusion counts differ for {method}")

    localization = read_table(localization_path, {
        "seed", "unit_substitution_rate", "array_fragments", "assembly_fraction",
        "family_id", "base_recall", "base_precision", "fragments", "truth_fragments",
    })
    localization_keys = {
        (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
         row["assembly_fraction"], row["family_id"])
        for row in localization
    }
    if (
        len(localization) != baseline_validation.get("localization_family_rows")
        or len(localization_keys) != len(localization)
    ):
        raise ValueError("Localization rows are missing or duplicated")
    return selected, localization, validation


def _source(result: list[dict], panel: str, metric: str, value: float, **context) -> None:
    row = {
        "panel": panel, "metric": metric, "method": "NA", "unit_substitution_rate": "NA",
        "array_fragments": "NA", "assembly_fraction": "NA", "value": value,
        "numerator": "NA", "denominator": "NA",
    }
    row.update(context)
    result.append(row)


def plot(
    metrics_path: Path,
    validation_path: Path,
    localization_path: Path,
    baseline_validation_path: Path,
    outdir: Path,
) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    rows, localization, validation = load_inputs(
        metrics_path, validation_path, localization_path, baseline_validation_path
    )
    outdir.mkdir(parents=True)
    grouped = {method: [row for row in rows if row["method"] == method] for method in METHODS}
    overall = {method: confusion(grouped[method]) for method in METHODS}
    divergences = sorted({float(row["unit_substitution_rate"]) for row in rows})
    fragment_counts = sorted({int(row["array_fragments"]) for row in rows})
    panel_source: list[dict] = []

    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans",
        "font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 3, figsize=(14.7, 8.8))
    fig.subplots_adjust(left=.055, right=.985, top=.91, bottom=.12, wspace=.18, hspace=.28)

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A   Predeclared domain-shift design", loc="left", fontweight="bold")
    design = [
        ("Fresh genomes", 3, "seeds 5401-5403"),
        ("Unit divergence", 3, "1%, 3%, 5% substitutions"),
        ("Array structures", 2, "one or three segments"),
        ("Paired conditions", len(grouped[BASELINE]), "all baseline and frozen rows"),
        ("Multi-k unavailable", validation["method_confusion"]["multik_loglinear"]["unavailable"],
         "retained as NA"),
    ]
    for index, (name, value, detail) in enumerate(design):
        y = 0.91 - index * 0.165
        ax.text(0.02, y, name, transform=ax.transAxes, fontweight="bold", va="center")
        ax.text(0.46, y, f"{value:,}", transform=ax.transAxes, color="#225588",
                fontsize=12.5, fontweight="bold", va="center")
        ax.text(0.62, y, detail, transform=ax.transAxes, va="center", fontsize=8)
        _source(panel_source, "A", name.lower().replace(" ", "_"), float(value),
                numerator=value, denominator=value)
    ax.text(0.02, 0.03, "Frozen exact-copy rule applied without domain-shift fitting",
            transform=ax.transAxes, color="#555555")

    ax = axes[0, 1]
    metrics = ("sensitivity", "false_positive_rate", "precision")
    x = np.arange(3)
    width = 0.34
    for offset, method in zip((-width/2, width/2), METHODS):
        values = [float(overall[method][metric]) for metric in metrics]
        bars = ax.bar(x+offset, values, width, color=COLORS[method], label=METHOD_LABELS[method])
        ax.bar_label(bars, [f"{value:.3f}" for value in values], padding=2, fontsize=7)
        for metric, value in zip(metrics, values):
            _source(panel_source, "B", metric, value, method=method)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylim(0, 1.12)
    ax.set_title("B   Overall domain-shift metrics", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=.15)

    ax = axes[0, 2]
    outcomes = ("TP", "FN", "FP", "TN")
    x = np.arange(4)
    for offset, method in zip((-width/2, width/2), METHODS):
        values = [int(overall[method][name]) for name in outcomes]
        bars = ax.bar(x+offset, values, width, color=COLORS[method])
        ax.bar_label(bars, fontsize=7, padding=2)
        for name, value in zip(outcomes, values):
            _source(panel_source, "C", name, float(value), method=method,
                    numerator=value, denominator=len(grouped[method]))
    ax.set_xticks(x, outcomes)
    ax.set_ylabel("Family conditions")
    ax.set_title("C   Confusion counts", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=.15)

    for panel, ax, metric, title in (
        ("D", axes[1, 0], "sensitivity", "D   Sensitivity by divergence"),
        ("E", axes[1, 1], "false_positive_rate", "E   False-positive rate by divergence"),
    ):
        for method in METHODS:
            for fragments in fragment_counts:
                values = []
                for divergence in divergences:
                    subset = [
                        row for row in grouped[method]
                        if float(row["unit_substitution_rate"]) == divergence
                        and int(row["array_fragments"]) == fragments
                    ]
                    value = float(confusion(subset)[metric])
                    values.append(value)
                    _source(panel_source, panel, metric, value, method=method,
                            unit_substitution_rate=f"{divergence:g}",
                            array_fragments=fragments)
                ax.plot(
                    np.array(divergences)*100, values, marker=MARKERS.get(fragments, "o"),
                    linestyle="-" if fragments == 1 else "--", color=COLORS[method], lw=2,
                    label=f"{METHOD_LABELS[method]}, {fragments} segment{'s' if fragments != 1 else ''}",
                )
        ax.set_xticks(np.array(divergences)*100, [f"{value*100:g}%" for value in divergences])
        ax.set_ylim(-.03, 1.05)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.grid(axis="y", alpha=.15)
    axes[1, 0].set_ylabel("Sensitivity")
    axes[1, 1].set_ylabel("False-positive rate")
    axes[1, 1].legend(frameon=False, fontsize=7, loc="lower right")

    ax = axes[1, 2]
    x = np.arange(len(divergences))
    width = .34
    for offset, fragments, color in zip((-width/2, width/2), fragment_counts, ("#4477aa", "#ccbb44")):
        values = []
        for divergence in divergences:
            subset = [
                row for row in localization
                if float(row["unit_substitution_rate"]) == divergence
                and int(row["array_fragments"]) == fragments
                and float(row["assembly_fraction"]) == 1.0
            ]
            numeric = [float(row["base_recall"]) for row in subset if row["base_recall"] != "NA"]
            value = sum(numeric)/len(numeric) if numeric else float("nan")
            values.append(value)
            _source(panel_source, "F", "mean_full_assembly_base_recall", value,
                    unit_substitution_rate=f"{divergence:g}", array_fragments=fragments,
                    assembly_fraction="1", numerator=len(numeric), denominator=len(subset))
        bars = ax.bar(x+offset, values, width, color=color,
                      label=f"{fragments} segment{'s' if fragments != 1 else ''}")
        ax.bar_label(bars, [f"{value:.3f}" if math.isfinite(value) else "NA" for value in values],
                     padding=2, fontsize=7)
    ax.set_xticks(x, [f"{value*100:g}%" for value in divergences])
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Mean base recall")
    ax.set_title("F   Full-assembly localization", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=.15)

    fig.suptitle("Frozen collapse rule fails to control false positives under array divergence",
                 fontsize=14, fontweight="bold")
    fig.text(
        .5, .025,
        "Known-catalogue, length-preserving substitutions and interrupted arrays; this is adverse simulation evidence, not biological truth.",
        ha="center", fontsize=8, color="#555555",
    )
    write_table(outdir/"panel_source.tsv", panel_source, list(panel_source[0]))
    outputs = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir/f"domain_shift_validation.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        outputs[path.name] = digest_file(path)
    plt.close(fig)
    provenance = {
        "complete": True,
        "methods": list(METHODS),
        "paired_conditions": len(grouped[BASELINE]),
        "unit_substitution_rates": divergences,
        "array_fragment_counts": fragment_counts,
        "input_sha256": {
            str(path.resolve()): digest_file(path)
            for path in (metrics_path, validation_path, localization_path, baseline_validation_path)
        },
        "script_sha256": digest_file(Path(__file__)),
        "panel_source_sha256": digest_file(outdir/"panel_source.tsv"),
        "outputs": outputs,
        "validation_method_confusion": validation["method_confusion"],
        "interpretation": (
            "predeclared domain-shift known-catalogue simulation; increased sensitivity accompanied by "
            "higher false-positive rate; not biological collapse truth"
        ),
    }
    (outdir/"figure_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--localization", required=True, type=Path)
    parser.add_argument("--baseline-validation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.metrics, args.validation, args.localization, args.baseline_validation, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
