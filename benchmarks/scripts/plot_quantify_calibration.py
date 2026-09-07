"""Render a six-panel audit of quantify depth-calibration development data."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, write_table


matplotlib.rcParams["svg.fonttype"] = "none"

METHODS = (
    "baseline_total_bases",
    "oracle_error_survival",
    "empirical_controls",
    "empirical_controls_plus_oracle_error",
)
LABELS = {
    "baseline_total_bases": "Total bases",
    "oracle_error_survival": "Oracle error",
    "empirical_controls": "Controls",
    "empirical_controls_plus_oracle_error": "Controls + oracle",
    "depth_gated_controls": "Depth-gated controls",
}
COLORS = {
    "baseline_total_bases": "#6b6b6b",
    "oracle_error_survival": "#4c78a8",
    "empirical_controls": "#59a14f",
    "empirical_controls_plus_oracle_error": "#f28e2b",
    "depth_gated_controls": "#b07aa1",
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot(evidence: Path, outdir: Path) -> dict[str, object]:
    evidence = evidence.resolve()
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    required = {
        "metrics": evidence / "run" / "metrics.tsv",
        "executions": evidence / "run" / "executions.tsv",
        "summary": evidence / "run" / "summary.tsv",
        "decision": evidence / "decision.json",
        "validation": evidence / "run" / "validation.json",
        "manifest": evidence / "archive_manifest.json",
    }
    if any(not path.is_file() for path in required.values()):
        raise ValueError("Quantify calibration evidence is incomplete")
    metrics = _read(required["metrics"])
    executions = _read(required["executions"])
    summaries = _read(required["summary"])
    decision = json.loads(required["decision"].read_text())
    validation = json.loads(required["validation"].read_text())
    if not validation.get("complete") or len(metrics) != 5940 or len(executions) != 108:
        raise ValueError("Quantify calibration matrix is incomplete")

    seeds = sorted({row["seed"] for row in metrics}, key=int)
    coverages = sorted({int(float(row["coverage"])) for row in metrics})
    by_key = {
        (row["seed"], row["condition_id"], row["family_id"], row["method"]): row
        for row in metrics
    }
    base_keys = sorted(
        key[:3] for key in by_key if key[3] == "baseline_total_bases"
    )
    threshold = float(decision["posthoc_candidate"]["control_mean_depth_threshold"])
    candidate_rows: list[dict[str, str]] = []
    for key in base_keys:
        control = by_key[(*key, "empirical_controls")]
        selected_method = (
            "empirical_controls"
            if float(control["control_mean_depth"]) >= threshold
            else "baseline_total_bases"
        )
        selected = dict(by_key[(*key, selected_method)])
        selected["method"] = "depth_gated_controls"
        candidate_rows.append(selected)

    by_method: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        by_method[row["method"]].append(row)
    by_method["depth_gated_controls"] = candidate_rows
    plotted_methods = (*METHODS, "depth_gated_controls")
    source_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(3, 2, figsize=(16, 15), constrained_layout=True)
    x = np.arange(len(plotted_methods))
    aggregate_mare = [
        statistics.mean(float(row["absolute_relative_error"]) for row in by_method[method])
        for method in plotted_methods
    ]
    axes[0, 0].bar(x, aggregate_mare, color=[COLORS[m] for m in plotted_methods], alpha=0.88)
    for index, method in enumerate(plotted_methods):
        for seed in seeds:
            values = [
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if row["seed"] == seed
            ]
            value = statistics.mean(values)
            axes[0, 0].scatter(index, value, color="black", s=22, zorder=3)
            source_rows.append({"panel": "A", "method": method, "stratum": seed, "metric": "mean_absolute_relative_error", "value": value})
    axes[0, 0].set_xticks(x, [LABELS[m] for m in plotted_methods], rotation=25, ha="right")
    axes[0, 0].set_ylabel("mean absolute relative error")
    axes[0, 0].set_title("A  Overall error (points are independent genomes)", loc="left", fontweight="bold")

    biases = [
        statistics.mean(float(row["signed_relative_error"]) for row in by_method[method])
        for method in METHODS
    ]
    axes[0, 1].bar(np.arange(len(METHODS)), biases, color=[COLORS[m] for m in METHODS], alpha=0.88)
    axes[0, 1].axhline(0, color="black", linewidth=0.8)
    axes[0, 1].set_xticks(np.arange(len(METHODS)), [LABELS[m] for m in METHODS], rotation=25, ha="right")
    axes[0, 1].set_ylabel("mean signed relative error")
    axes[0, 1].set_title("B  Residual copy-number bias", loc="left", fontweight="bold")
    for method, value in zip(METHODS, biases):
        source_rows.append({"panel": "B", "method": method, "stratum": "all", "metric": "mean_signed_relative_error", "value": value})

    for method in plotted_methods:
        values = []
        for coverage in coverages:
            rows = [row for row in by_method[method] if int(float(row["coverage"])) == coverage]
            value = statistics.mean(float(row["absolute_relative_error"]) for row in rows)
            values.append(value)
            source_rows.append({"panel": "C", "method": method, "stratum": coverage, "metric": "mean_absolute_relative_error", "value": value})
        axes[1, 0].plot(coverages, values, marker="o", color=COLORS[method], label=LABELS[method])
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xticks(coverages, [f"{value}×" for value in coverages])
    axes[1, 0].set_xlabel("nominal read coverage")
    axes[1, 0].set_ylabel("mean absolute relative error")
    axes[1, 0].set_title("C  Error depends on coverage", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)

    paired = {row["method"]: row for row in decision["paired_against_baseline"]}
    alternatives = tuple(method for method in METHODS if method != "baseline_total_bases")
    bottoms = np.zeros(len(alternatives))
    for label, field, color in (
        ("improved", "improved", "#59a14f"),
        ("equal", "equal", "#bab0ac"),
        ("worse", "worse", "#e15759"),
    ):
        fractions = np.array([
            float(paired[method][field]) / 1485 for method in alternatives
        ])
        axes[1, 1].bar(np.arange(len(alternatives)), fractions, bottom=bottoms, color=color, label=label)
        bottoms += fractions
        for method, value in zip(alternatives, fractions):
            source_rows.append({"panel": "D", "method": method, "stratum": field, "metric": "paired_family_fraction", "value": value})
    axes[1, 1].set_xticks(np.arange(len(alternatives)), [LABELS[m] for m in alternatives], rotation=25, ha="right")
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_ylabel("fraction of 1,485 paired family conditions")
    axes[1, 1].set_title("D  Improvements coexist with regressions", loc="left", fontweight="bold")
    axes[1, 1].legend(frameon=False, fontsize=8)

    execution_by_method: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in executions:
        execution_by_method[row["method"]].append(row)
    for method in METHODS:
        values = []
        for coverage in coverages:
            rows = [row for row in execution_by_method[method] if int(float(row["coverage"])) == coverage]
            value = statistics.median(float(row["runtime_seconds"]) for row in rows)
            values.append(value)
            source_rows.append({"panel": "E", "method": method, "stratum": coverage, "metric": "median_runtime_seconds", "value": value})
        axes[2, 0].plot(coverages, values, marker="o", color=COLORS[method], label=LABELS[method])
    axes[2, 0].set_xscale("log")
    axes[2, 0].set_yscale("log")
    axes[2, 0].set_xticks(coverages, [f"{value}×" for value in coverages])
    axes[2, 0].set_xlabel("nominal read coverage")
    axes[2, 0].set_ylabel("median wall time (s; log scale)")
    axes[2, 0].set_title("E  Original error-survival path adds runtime", loc="left", fontweight="bold")
    axes[2, 0].legend(frameon=False, fontsize=8, ncol=2)

    width = 0.18
    for index, method in enumerate(METHODS):
        values = []
        for coverage in coverages:
            rows = [
                row for row in summaries
                if row["method"] == method and int(float(row["coverage"])) == coverage
            ]
            value = statistics.mean(float(row["diagnostic_spread_truth_coverage"]) for row in rows)
            values.append(value)
            source_rows.append({"panel": "F", "method": method, "stratum": coverage, "metric": "diagnostic_spread_truth_coverage", "value": value})
        positions = np.arange(len(coverages)) + (index - 1.5) * width
        axes[2, 1].bar(positions, values, width=width, color=COLORS[method], label=LABELS[method])
    axes[2, 1].axhline(0.95, color="#d62728", linestyle="--", linewidth=1.2, label="nominal 0.95")
    axes[2, 1].set_xticks(np.arange(len(coverages)), [f"{value}×" for value in coverages])
    axes[2, 1].set_ylim(0, 1)
    axes[2, 1].set_xlabel("nominal read coverage")
    axes[2, 1].set_ylabel("truth within diagnostic 10th–90th spread")
    axes[2, 1].set_title("F  Diagnostic spread is not a calibrated CI", loc="left", fontweight="bold")
    axes[2, 1].legend(frameon=False, fontsize=8, ncol=2)

    fig.suptitle(
        "Public quantify calibration: lower aggregate error, low-depth trade-offs, and an uncalibrated spread",
        fontsize=15,
        fontweight="bold",
    )
    outdir.mkdir(parents=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"quantify_calibration.{extension}", dpi=180)
    plt.close(fig)
    write_table(outdir / "panel_source.tsv", source_rows, ["panel", "method", "stratum", "metric", "value"])
    (outdir / "figure_legend.md").write_text(
        "**Quantify calibration development audit.** A, mean absolute relative "
        "copy-number error across 1,485 family conditions per method; points are "
        "the three independently simulated genomes. The depth-gated control rule "
        "was selected post hoc from this development matrix. B, signed error, "
        "showing residual underestimation. C, error by nominal coverage; empirical "
        "controls regress at 1× and improve 5×/20× results. D, paired outcomes "
        "against total-bases normalization. E, descriptive single-execution wall "
        "time from the original implementation; the oracle rate is unavailable in "
        "blind data. F, truth inclusion within exported 10th–90th diagnostic-k-mer "
        "spread, demonstrating that it is not a calibrated sampling confidence "
        "interval. All conditions use known catalogues and simulation-truth-assisted "
        "control selection. No method is promoted without untouched validation.\n"
    )
    provenance = {
        "benchmark_id": "quantify_calibration_development_v1",
        "inputs": {name: {"path": str(path), "sha256": digest_file(path)} for name, path in required.items()},
        "outputs": {
            path.name: digest_file(path)
            for path in sorted(outdir.iterdir())
            if path.name != "figure_provenance.json"
        },
        "warning": "development_only;posthoc_candidate;oracle_not_blind;diagnostic_spread_not_sampling_CI",
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.evidence, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
