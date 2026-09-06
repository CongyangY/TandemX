"""Render six-panel evidence for the divergence-aware localizer validation."""
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


DATASETS = ("development_v1", "development_v2", "heldout")
LABELS = {
    "development_v1": "Development v1",
    "development_v2": "Development v2",
    "heldout": "Held-out v1",
}
COLORS = {
    "development_v1": "#777777",
    "development_v2": "#4477aa",
    "heldout": "#b23a48",
}
BASELINE = "single_k21"
FROZEN = "multik_fallback_depth_rule"


def _load_localization(root: Path) -> tuple[list[dict], dict, dict, dict]:
    validation = json.loads((root / "validation.json").read_text())
    config = json.loads((root / "run_config.json").read_text())
    environment = json.loads((root / "environment.json").read_text())
    rows = read_table(root / "localization_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "assembly_fraction",
        "family_id", "true_assembly_bp", "predicted_assembly_bp", "base_recall",
        "base_precision", "fragments", "truth_fragments",
    })
    split = environment.get("split")
    seeds = {str(seed) for seed in config.get("seeds", {}).get(split, [])}
    keys = {
        (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
         row["assembly_fraction"], row["family_id"])
        for row in rows
    }
    if (
        validation.get("complete") is not True
        or validation.get("localization_family_rows") != len(rows)
        or validation.get("successful") != validation.get("executions")
        or not seeds or {row["seed"] for row in rows} != seeds
        or len(keys) != len(rows)
    ):
        raise ValueError(f"Incomplete or duplicated localization evidence: {root}")
    return rows, validation, config, environment


def _confusion(rows: list[dict]) -> dict[str, float | int]:
    counts = Counter(row["outcome"] for row in rows)
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    return {
        "TP": tp, "FN": fn, "FP": fp, "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else math.nan,
        "false_positive_rate": fp / (fp + tn) if fp + tn else math.nan,
        "precision": tp / (tp + fp) if tp + fp else math.nan,
    }


def load_inputs(
    development_v1: Path,
    development_v2: Path,
    heldout: Path,
    multik: Path,
) -> tuple[dict[str, list[dict]], dict, list[dict], dict]:
    roots = {
        "development_v1": development_v1,
        "development_v2": development_v2,
        "heldout": heldout,
    }
    loaded = {name: _load_localization(root) for name, root in roots.items()}
    localizations = {name: value[0] for name, value in loaded.items()}
    development_keys = []
    development_truth = []
    for name in ("development_v1", "development_v2"):
        rows = localizations[name]
        development_keys.append({
            (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
             row["assembly_fraction"], row["family_id"])
            for row in rows
        })
        development_truth.append({
            (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
             row["assembly_fraction"], row["family_id"]):
            (row["true_assembly_bp"], row["truth_fragments"])
            for row in rows
        })
    development_seeds = {row["seed"] for row in localizations["development_v1"]}
    heldout_seeds = {row["seed"] for row in localizations["heldout"]}
    if (
        development_keys[0] != development_keys[1]
        or development_truth[0] != development_truth[1]
        or development_seeds & heldout_seeds
    ):
        raise ValueError("Development versions are not paired or held-out seeds overlap")

    multik_validation = json.loads((multik / "validation.json").read_text())
    multik_environment = json.loads((multik / "environment.json").read_text())
    multik_rows = read_table(multik / "comparison_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
    })
    selected = [row for row in multik_rows if row["method"] in {BASELINE, FROZEN}]
    grouped = {method: [row for row in selected if row["method"] == method]
               for method in (BASELINE, FROZEN)}
    paired_keys = {
        method: {
            (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
             row["coverage"], row["substitution_rate"], row["assembly_fraction"],
             row["family_id"])
            for row in group
        }
        for method, group in grouped.items()
    }
    heldout_validation = loaded["heldout"][1]
    heldout_environment = loaded["heldout"][3]
    if (
        multik_validation.get("complete") is not True
        or multik_validation.get("split") != "heldout"
        or paired_keys[BASELINE] != paired_keys[FROZEN]
        or len(paired_keys[BASELINE]) != len(grouped[BASELINE])
        or {row["seed"] for row in selected} != heldout_seeds
        or multik_environment.get("previous_validation_sha256")
            != digest_file(heldout / "validation.json")
        or multik_environment.get("previous_environment_sha256")
            != digest_file(heldout / "environment.json")
        or multik_environment.get("previous_config_sha256")
            != digest_file(heldout / "run_config.json")
        or heldout_environment.get("split") != "heldout"
        or heldout_validation.get("comparison_family_rows") != len(grouped[BASELINE])
    ):
        raise ValueError("Held-out localization and multi-k evidence are not paired")
    for method, group in grouped.items():
        observed = _confusion(group)
        expected = multik_validation.get("method_confusion", {}).get(method, {})
        if any(observed[name] != expected.get(name) for name in ("TP", "FN", "FP", "TN")):
            raise ValueError(f"Confusion receipt differs for {method}")
    metadata = {
        name: {"validation": value[1], "config": value[2], "environment": value[3]}
        for name, value in loaded.items()
    }
    return localizations, metadata, selected, multik_validation


def _source(result: list[dict], panel: str, metric: str, value: float, **context) -> None:
    row = {
        "panel": panel, "metric": metric, "dataset": "NA", "method": "NA",
        "unit_substitution_rate": "NA", "array_fragments": "NA",
        "assembly_fraction": "NA", "value": value,
        "numerator": "NA", "denominator": "NA",
    }
    row.update(context)
    result.append(row)


def _mean(rows: list[dict], field: str) -> float:
    values = [float(row[field]) for row in rows if row[field] not in {"", "NA", "None"}]
    if not values:
        raise ValueError(f"No finite values for {field}")
    return sum(values) / len(values)


def plot(
    development_v1: Path,
    development_v2: Path,
    heldout: Path,
    multik: Path,
    outdir: Path,
) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    localizations, metadata, comparisons, multik_validation = load_inputs(
        development_v1, development_v2, heldout, multik
    )
    outdir.mkdir(parents=True)
    divergences = sorted({float(row["unit_substitution_rate"])
                          for row in localizations["heldout"]})
    fragment_counts = sorted({int(row["array_fragments"])
                              for row in localizations["heldout"]})
    fractions = sorted({float(row["assembly_fraction"])
                        for row in localizations["heldout"] if float(row["assembly_fraction"]) > 0})
    panel_source: list[dict] = []

    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans",
        "font.size": 8.2, "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.8))
    fig.subplots_adjust(left=.06, right=.985, top=.91, bottom=.115, wspace=.22, hspace=.30)

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A   Predeclared validation design", loc="left", fontweight="bold")
    design = (
        ("Failed development", 3, "seeds 5301-5303; v1"),
        ("Selected development", 3, "same seeds; bounded bridge"),
        ("Held-out genomes", 3, "fresh seeds 5501-5503"),
        ("Array scenarios", 6, "1/3/5% x 1/3 segments"),
        ("Held-out commands", metadata["heldout"]["validation"]["executions"],
         "locate + quantify + compare"),
    )
    for index, (name, value, detail) in enumerate(design):
        y = .91 - index * .165
        ax.text(.02, y, name, transform=ax.transAxes, fontweight="bold", va="center")
        ax.text(.49, y, f"{value:,}", transform=ax.transAxes, color="#225588",
                fontsize=12.3, fontweight="bold", va="center")
        ax.text(.68, y, detail, transform=ax.transAxes, va="center", fontsize=7.3)
        _source(panel_source, "A", name.lower().replace(" ", "_"), float(value),
                numerator=value, denominator=value)
    ax.text(.02, .025, "Known catalogue; substitutions only; no biological truth",
            transform=ax.transAxes, color="#555555", fontsize=7.7)

    ax = axes[0, 1]
    for dataset in DATASETS:
        for fragments in fragment_counts:
            values = []
            for divergence in divergences:
                subset = [row for row in localizations[dataset]
                          if float(row["assembly_fraction"]) == 1
                          and float(row["unit_substitution_rate"]) == divergence
                          and int(row["array_fragments"]) == fragments]
                value = _mean(subset, "base_recall")
                values.append(value)
                _source(panel_source, "B", "mean_full_assembly_base_recall", value,
                        dataset=dataset, unit_substitution_rate=divergence,
                        array_fragments=fragments, assembly_fraction=1,
                        numerator=len(subset), denominator=len(subset))
            ax.plot(np.array(divergences) * 100, values, color=COLORS[dataset],
                    linestyle="-" if fragments == 1 else "--",
                    marker="o" if fragments == 1 else "s", lw=1.8, ms=4.5,
                    label=f"{LABELS[dataset]}, {fragments} segment{'s' if fragments != 1 else ''}")
    ax.axhline(.95, color="#222222", lw=1, linestyle=":", label="0.95 gate")
    ax.set_xticks(np.array(divergences) * 100, [f"{value*100:g}%" for value in divergences])
    ax.set_ylim(.70, 1.015)
    ax.set_ylabel("Mean base recall")
    ax.set_title("B   Full-assembly localization recall", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=.15)
    ax.legend(frameon=False, fontsize=6.4, ncol=2, loc="lower left")

    ax = axes[0, 2]
    x = np.arange(len(divergences))
    width = .24
    for offset, dataset in zip((-width, 0, width), DATASETS):
        values = []
        for divergence in divergences:
            subset = [row for row in localizations[dataset]
                      if float(row["assembly_fraction"]) == 1
                      and float(row["unit_substitution_rate"]) == divergence]
            value = sum(abs(int(row["fragments"])-int(row["truth_fragments"]))
                        for row in subset) / len(subset)
            values.append(value)
            _source(panel_source, "C", "mean_absolute_fragment_count_error", value,
                    dataset=dataset, unit_substitution_rate=divergence,
                    assembly_fraction=1, numerator=len(subset), denominator=len(subset))
        bars = ax.bar(x+offset, values, width, color=COLORS[dataset], label=LABELS[dataset])
        ax.bar_label(bars, [f"{value:.2f}" for value in values], fontsize=6.3, padding=2)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(x, [f"{value*100:g}%" for value in divergences])
    ax.set_ylabel("Mean absolute count error (symlog)")
    ax.set_title("C   Predicted array fragmentation", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 0]
    heldout_rows = localizations["heldout"]
    for divergence, marker in zip(divergences, ("o", "s", "^")):
        values = []
        for fraction in fractions:
            subset = [row for row in heldout_rows
                      if float(row["unit_substitution_rate"]) == divergence
                      and float(row["assembly_fraction"]) == fraction]
            value = _mean(subset, "base_recall")
            values.append(value)
            _source(panel_source, "D", "mean_base_recall", value, dataset="heldout",
                    unit_substitution_rate=divergence, assembly_fraction=fraction,
                    numerator=len(subset), denominator=len(subset))
        ax.plot(np.array(fractions)*100, values, marker=marker, lw=1.8,
                label=f"{divergence*100:g}% divergence")
    ax.axhline(.95, color="#222222", lw=1, linestyle=":")
    ax.set_xticks(np.array(fractions)*100, [f"{value*100:g}%" for value in fractions])
    ax.set_ylim(.70, 1.015)
    ax.set_xlabel("Assembly repeat retention")
    ax.set_ylabel("Mean base recall")
    ax.set_title("D   Held-out recall across collapse", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 1]
    for divergence, marker in zip(divergences, ("o", "s", "^")):
        values = []
        for fraction in fractions:
            subset = [row for row in heldout_rows
                      if float(row["unit_substitution_rate"]) == divergence
                      and float(row["assembly_fraction"]) == fraction]
            errors = [abs(int(row["predicted_assembly_bp"])-int(row["true_assembly_bp"]))
                      / int(row["true_assembly_bp"]) for row in subset]
            value = sum(errors)/len(errors)
            values.append(value)
            _source(panel_source, "E", "mean_absolute_relative_assembly_bp_error", value,
                    dataset="heldout", unit_substitution_rate=divergence,
                    assembly_fraction=fraction, numerator=len(subset), denominator=len(subset))
        ax.plot(np.array(fractions)*100, values, marker=marker, lw=1.8,
                label=f"{divergence*100:g}% divergence")
    ax.set_xticks(np.array(fractions)*100, [f"{value*100:g}%" for value in fractions])
    ax.set_xlabel("Assembly repeat retention")
    ax.set_ylabel("Mean absolute relative error")
    ax.set_title("E   Held-out localized repeat-bp error", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="y", alpha=.15)

    ax = axes[1, 2]
    grouped = {method: [row for row in comparisons if row["method"] == method]
               for method in (BASELINE, FROZEN)}
    overall = {method: _confusion(rows) for method, rows in grouped.items()}
    metric_names = ("sensitivity", "false_positive_rate", "precision")
    x = np.arange(len(metric_names))
    width = .34
    for offset, method, label, color in (
        (-width/2, BASELINE, "Single k=21", "#777777"),
        (width/2, FROZEN, "Frozen multi-k", "#b23a48"),
    ):
        values = [float(overall[method][name]) for name in metric_names]
        bars = ax.bar(x+offset, values, width, color=color, label=label)
        ax.bar_label(bars, [f"{value:.3f}" for value in values], fontsize=6.8, padding=2)
        for name, value in zip(metric_names, values):
            if name == "sensitivity":
                numerator = overall[method]["TP"]
                denominator = overall[method]["TP"] + overall[method]["FN"]
            elif name == "false_positive_rate":
                numerator = overall[method]["FP"]
                denominator = overall[method]["FP"] + overall[method]["TN"]
            else:
                numerator = overall[method]["TP"]
                denominator = overall[method]["TP"] + overall[method]["FP"]
            _source(panel_source, "F", name, value, method=method,
                    numerator=numerator, denominator=denominator)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylim(0, 1.12)
    ax.set_title("F   Assembly/read classification", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="y", alpha=.15)

    fig.suptitle("Divergence-aware anchor bridging passes predeclared localization gates",
                 fontsize=14, fontweight="bold")
    fig.text(
        .5, .025,
        "Fresh seeds validate localization; the frozen multi-k classifier gains sensitivity while increasing false positives.",
        ha="center", fontsize=8, color="#555555",
    )
    write_table(outdir / "panel_source.tsv", panel_source, list(panel_source[0]))
    outputs = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir / f"localizer_validation.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        outputs[path.name] = digest_file(path)
    plt.close(fig)
    inputs = []
    for root in (development_v1, development_v2, heldout, multik):
        for name in ("validation.json", "environment.json", "run_config.json"):
            path = root / name
            inputs.append(path)
    inputs.extend([
        development_v1 / "localization_metrics.tsv",
        development_v2 / "localization_metrics.tsv",
        heldout / "localization_metrics.tsv",
        multik / "comparison_metrics.tsv",
    ])
    provenance = {
        "complete": True,
        "development_versions_exactly_paired": True,
        "development_and_heldout_seeds_disjoint": True,
        "heldout_localization_gates": {
            "full_assembly_mean_base_recall": _mean(
                [row for row in heldout_rows if float(row["assembly_fraction"]) == 1], "base_recall"
            ),
            "positive_assembly_mean_base_precision": _mean(
                [row for row in heldout_rows if float(row["assembly_fraction"]) > 0], "base_precision"
            ),
            "absent_family_false_positive_rate": (
                sum(int(row["predicted_assembly_bp"]) > 0 for row in heldout_rows
                    if float(row["assembly_fraction"]) == 0)
                / sum(float(row["assembly_fraction"]) == 0 for row in heldout_rows)
            ),
        },
        "method_confusion": multik_validation["method_confusion"],
        "input_sha256": {str(path.resolve()): digest_file(path) for path in inputs},
        "script_sha256": digest_file(Path(__file__)),
        "panel_source_sha256": digest_file(outdir / "panel_source.tsv"),
        "outputs": outputs,
        "interpretation": (
            "predeclared known-catalogue substitution simulation; localization gates passed on fresh seeds; "
            "the frozen multi-k classifier trades a sensitivity gain for a higher false-positive rate; "
            "not alignment identity or biological validation"
        ),
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-v1", required=True, type=Path)
    parser.add_argument("--development-v2", required=True, type=Path)
    parser.add_argument("--heldout", required=True, type=Path)
    parser.add_argument("--multik", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.development_v1, args.development_v2, args.heldout, args.multik, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
