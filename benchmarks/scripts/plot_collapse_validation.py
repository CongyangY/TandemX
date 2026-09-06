"""Render the paired held-out assembly-comparison validation as six panels."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, read_table, write_table


BASELINE = "single_k21"
FROZEN = "multik_fallback_depth_rule"
METHOD_LABELS = {BASELINE: "Single k=21", FROZEN: "Frozen multi-k rule"}
COLORS = {BASELINE: "#777777", FROZEN: "#228833"}


def confusion(rows: list[dict]) -> dict[str, float | int]:
    counts = Counter(row["outcome"] for row in rows)
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    return {
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else float("nan"),
        "false_positive_rate": fp / (fp + tn) if fp + tn else float("nan"),
        "precision": tp / (tp + fp) if tp + fp else float("nan"),
    }


def load_pairs(metrics_path: Path, validation_path: Path) -> tuple[list[dict], dict]:
    rows = read_table(metrics_path, {
        "seed", "coverage", "substitution_rate", "assembly_fraction",
        "family_id", "method", "outcome",
    })
    validation = json.loads(validation_path.read_text())
    if validation.get("complete") is not True or validation.get("split") != "heldout":
        raise ValueError("Require complete held-out validation")
    selected = [row for row in rows if row["method"] in {BASELINE, FROZEN}]
    by_method = {method: [row for row in selected if row["method"] == method]
                 for method in (BASELINE, FROZEN)}
    keys = {
        method: {(row["seed"], row["coverage"], row["substitution_rate"],
                  row["assembly_fraction"], row["family_id"]) for row in group}
        for method, group in by_method.items()
    }
    if not selected or len(keys[BASELINE]) != len(by_method[BASELINE]) or keys[BASELINE] != keys[FROZEN]:
        raise ValueError("Baseline and frozen-rule conditions are not exactly paired")
    for method, group in by_method.items():
        observed = confusion(group)
        expected = validation.get("method_confusion", {}).get(method, {})
        if any(observed[name] != expected.get(name) for name in ("TP", "FN", "FP", "TN")):
            raise ValueError(f"Validation confusion counts differ for {method}")
    return selected, validation


def _add_source(
    result: list[dict], panel: str, metric: str, value: float,
    *, method: str = "NA", seed: str = "NA", coverage: str = "NA",
    substitution_rate: str = "NA", assembly_fraction: str = "NA",
    numerator: int | str = "NA", denominator: int | str = "NA",
) -> None:
    result.append({
        "panel": panel, "metric": metric, "method": method, "seed": seed,
        "coverage": coverage, "substitution_rate": substitution_rate,
        "assembly_fraction": assembly_fraction, "value": value,
        "numerator": numerator, "denominator": denominator,
    })


def plot(metrics_path: Path, validation_path: Path, outdir: Path) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    rows, validation = load_pairs(metrics_path, validation_path)
    outdir.mkdir(parents=True)
    methods = (BASELINE, FROZEN)
    grouped = {method: [row for row in rows if row["method"] == method] for method in methods}
    overall = {method: confusion(grouped[method]) for method in methods}
    seeds = sorted({row["seed"] for row in rows}, key=int)
    panel_source: list[dict] = []

    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans",
        "font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.4), layout="constrained")

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A   Predeclared validation design", loc="left", fontweight="bold")
    design = [
        ("Fresh genomes", "3", "5201-5203"),
        ("Coverage", "3", "1x, 5x, 20x"),
        ("Read-error tiers", "3", "0, 0.1%, 1% substitutions"),
        ("Assembly fractions", "5", "100%, 75%, 50%, 25%, 0%"),
        ("Family conditions", "405", "243 positive + 162 control"),
    ]
    for index, (name, count, detail) in enumerate(design):
        y = 0.9 - index * 0.16
        ax.text(0.02, y, name, transform=ax.transAxes, fontweight="bold", va="center")
        ax.text(0.43, y, count, transform=ax.transAxes, color="#225588",
                fontsize=13, fontweight="bold", va="center")
        ax.text(0.56, y, detail, transform=ax.transAxes, va="center")
        _add_source(panel_source, "A", name.replace(" ", "_").lower(), float(count),
                    numerator=int(count), denominator=int(count))
    ax.text(0.02, 0.05, "Rule and hashes committed before execution; no held-out fitting",
            transform=ax.transAxes, color="#555555")

    ax = axes[0, 1]
    metrics = ("sensitivity", "false_positive_rate", "precision")
    labels = ("Sensitivity", "False-positive rate", "Precision")
    x = np.arange(len(metrics))
    width = 0.34
    for offset, method in zip((-width / 2, width / 2), methods):
        values = [float(overall[method][metric]) for metric in metrics]
        bars = ax.bar(x + offset, values, width, label=METHOD_LABELS[method], color=COLORS[method])
        ax.bar_label(bars, labels=[f"{value:.3f}" for value in values], padding=2, fontsize=7)
        for metric, value in zip(metrics, values):
            _add_source(panel_source, "B", metric, value, method=method)
    ax.set_xticks(x, ("Sensitivity", "False-positive\nrate", "Precision"))
    ax.set_ylim(0, 1.12)
    ax.set_title("B   Overall classification metrics", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.15)

    ax = axes[0, 2]
    outcomes = ("TP", "FN", "FP", "TN")
    x = np.arange(len(outcomes))
    for offset, method in zip((-width / 2, width / 2), methods):
        values = [int(overall[method][name]) for name in outcomes]
        bars = ax.bar(x + offset, values, width, label=METHOD_LABELS[method], color=COLORS[method])
        ax.bar_label(bars, fontsize=7, padding=2)
        for name, value in zip(outcomes, values):
            _add_source(panel_source, "C", name, float(value), method=method,
                        numerator=value, denominator=len(grouped[method]))
    ax.set_xticks(x, outcomes)
    ax.set_ylabel("Family conditions")
    ax.set_title("C   Confusion counts", loc="left", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.15)

    for panel, axis, metric, title in (
        ("D", axes[1, 0], "sensitivity", "D   Sensitivity by fresh genome"),
        ("E", axes[1, 1], "false_positive_rate", "E   False-positive rate by fresh genome"),
    ):
        seed_x = np.arange(len(seeds))
        for method in methods:
            values = []
            for seed in seeds:
                subset = [row for row in grouped[method] if row["seed"] == seed]
                value = float(confusion(subset)[metric])
                values.append(value)
                _add_source(panel_source, panel, metric, value, method=method, seed=seed)
            axis.plot(seed_x, values, marker="o", lw=2, color=COLORS[method],
                      label=METHOD_LABELS[method])
        axis.set_xticks(seed_x, [f"s{seed}" for seed in seeds])
        axis.set_ylim(-0.03, 1.05)
        axis.set_title(title, loc="left", fontweight="bold")
        axis.grid(axis="y", alpha=0.15)
        axis.legend(frameon=False, loc="best")

    ax = axes[1, 2]
    coverages = (1.0, 5.0, 20.0)
    errors = (0.0, 0.001, 0.01)
    gains = np.zeros((len(coverages), len(errors)))
    annotations: list[list[str]] = []
    for row_index, coverage in enumerate(coverages):
        labels_row = []
        for column_index, error in enumerate(errors):
            values = {}
            for method in methods:
                subset = [row for row in grouped[method]
                          if float(row["coverage"]) == coverage
                          and float(row["substitution_rate"]) == error
                          and float(row["assembly_fraction"]) == 0.5]
                values[method] = confusion(subset)
            baseline_value = float(values[BASELINE]["sensitivity"])
            frozen_value = float(values[FROZEN]["sensitivity"])
            gains[row_index, column_index] = frozen_value - baseline_value
            labels_row.append(f"{baseline_value:.2f}->{frozen_value:.2f}")
            _add_source(panel_source, "F", "sensitivity_gain", gains[row_index, column_index],
                        method=FROZEN, coverage=f"{coverage:g}",
                        substitution_rate=f"{error:g}", assembly_fraction="0.5")
        annotations.append(labels_row)
    limit = max(0.01, float(np.max(np.abs(gains))))
    mesh = ax.pcolormesh(
        np.arange(len(errors) + 1), np.arange(len(coverages) + 1), gains,
        cmap="PiYG", vmin=-limit, vmax=limit, shading="flat", rasterized=False,
    )
    ax.set_ylim(len(coverages), 0)
    for i in range(len(coverages)):
        for j in range(len(errors)):
            ax.text(j + 0.5, i + 0.5, annotations[i][j], ha="center", va="center", fontsize=7,
                    color="white" if abs(gains[i, j]) > limit * 0.55 else "black")
    ax.set_xticks(np.arange(len(errors)) + 0.5, ["0", "0.1%", "1%"])
    ax.set_yticks(np.arange(len(coverages)) + 0.5, ["1x", "5x", "20x"])
    ax.set_xlabel("Substitution rate")
    ax.set_ylabel("Coverage")
    ax.set_title("F   Recall gain at 50% retention", loc="left", fontweight="bold")
    colorbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    colorbar.set_label("Sensitivity difference")

    fig.suptitle("Frozen multi-k rule: fresh conditional assembly validation",
                 fontsize=14, fontweight="bold")
    fig.text(0.5, 0.003,
             "Exact-copy known-catalogue simulations; improvement does not establish biological collapse performance.",
             ha="center", fontsize=8, color="#555555")
    write_table(outdir / "panel_source.tsv", panel_source, list(panel_source[0]))
    outputs = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir / f"collapse_validation.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        outputs[path.name] = digest_file(path)
    plt.close(fig)
    provenance = {
        "complete": True,
        "methods": list(methods),
        "paired_conditions": len(grouped[BASELINE]),
        "fresh_seeds": seeds,
        "input_sha256": {
            str(metrics_path.resolve()): digest_file(metrics_path),
            str(validation_path.resolve()): digest_file(validation_path),
        },
        "script_sha256": digest_file(Path(__file__)),
        "panel_source_sha256": digest_file(outdir / "panel_source.tsv"),
        "outputs": outputs,
        "interpretation": (
            "predeclared exact-copy known-catalogue simulation; not biological collapse truth"
        ),
        "validation_method_confusion": validation["method_confusion"],
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.metrics, args.validation, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
