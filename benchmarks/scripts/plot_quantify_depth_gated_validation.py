"""Render a six-panel audit of frozen depth-gated quantify validation."""
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
    "empirical_controls",
    "depth_gated_controls",
)
LABELS = {
    "baseline_total_bases": "Total bases",
    "empirical_controls": "Controls",
    "depth_gated_controls": "Depth-gated",
}
COLORS = {
    "baseline_total_bases": "#6b6b6b",
    "empirical_controls": "#59a14f",
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
        "metrics": evidence / "run/metrics.tsv",
        "paired": evidence / "run/paired.tsv",
        "executions": evidence / "run/executions.tsv",
        "gates": evidence / "run/gate_results.json",
        "validation": evidence / "run/validation.json",
        "decision": evidence / "decision.json",
        "manifest": evidence / "archive_manifest.json",
    }
    if any(not path.is_file() for path in required.values()):
        raise ValueError("Depth-gated validation evidence is incomplete")
    metrics = _read(required["metrics"])
    paired = _read(required["paired"])
    executions = _read(required["executions"])
    gates = json.loads(required["gates"].read_text())
    validation = json.loads(required["validation"].read_text())
    decision = json.loads(required["decision"].read_text())
    if (
        not validation.get("complete")
        or validation.get("promotion_status") != "passed"
        or gates.get("status") != "passed"
        or len(metrics) != 4455
        or len(paired) != 1485
        or len(executions) != 54
    ):
        raise ValueError("Depth-gated validation matrix is not the complete passing run")

    by_method: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        by_method[row["method"]].append(row)
    if set(by_method) != set(METHODS):
        raise ValueError("Unexpected methods in validation metrics")
    seeds = sorted({row["seed"] for row in metrics}, key=int)
    coverages = sorted({int(float(row["coverage"])) for row in metrics})
    source_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(3, 2, figsize=(16, 15), constrained_layout=True)

    aggregates = [
        statistics.mean(float(row["absolute_relative_error"]) for row in by_method[method])
        for method in METHODS
    ]
    x = np.arange(len(METHODS))
    axes[0, 0].bar(x, aggregates, color=[COLORS[method] for method in METHODS], alpha=0.88)
    for index, (method, aggregate) in enumerate(zip(METHODS, aggregates)):
        axes[0, 0].text(index, aggregate + 0.012, f"{aggregate:.3f}", ha="center", fontsize=9)
        source_rows.append(
            {
                "panel": "A",
                "method": method,
                "stratum": "aggregate",
                "metric": "mean_absolute_relative_error",
                "value": aggregate,
            }
        )
        for seed in seeds:
            value = statistics.mean(
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if row["seed"] == seed
            )
            axes[0, 0].scatter(index, value, color="black", s=24, zorder=3)
            source_rows.append(
                {
                    "panel": "A",
                    "method": method,
                    "stratum": seed,
                    "metric": "seed_mean_absolute_relative_error",
                    "value": value,
                }
            )
    axes[0, 0].set_xticks(x, [LABELS[method] for method in METHODS])
    axes[0, 0].set_ylim(0, 0.5)
    axes[0, 0].set_ylabel("mean absolute relative error")
    axes[0, 0].set_title("A  Frozen validation error", loc="left", fontweight="bold")

    seed_colors = dict(zip(seeds, ("#4c78a8", "#f28e2b", "#e15759")))
    for seed in seeds:
        values = []
        for method in ("baseline_total_bases", "depth_gated_controls"):
            value = statistics.mean(
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if row["seed"] == seed
            )
            values.append(value)
            source_rows.append(
                {
                    "panel": "B",
                    "method": method,
                    "stratum": seed,
                    "metric": "seed_mean_absolute_relative_error",
                    "value": value,
                }
            )
        axes[0, 1].plot(
            (0, 1), values, marker="o", linewidth=1.8, color=seed_colors[seed], label=f"seed {seed}"
        )
    axes[0, 1].set_xticks((0, 1), ("Total bases", "Depth-gated"))
    axes[0, 1].set_xlim(-0.25, 1.25)
    axes[0, 1].set_ylim(0.25, 0.47)
    axes[0, 1].set_ylabel("mean absolute relative error")
    axes[0, 1].set_title("B  Every validation genome improves", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False, fontsize=8)

    for method in METHODS:
        values = []
        for coverage in coverages:
            value = statistics.mean(
                float(row["absolute_relative_error"])
                for row in by_method[method]
                if int(float(row["coverage"])) == coverage
            )
            values.append(value)
            source_rows.append(
                {
                    "panel": "C",
                    "method": method,
                    "stratum": coverage,
                    "metric": "coverage_mean_absolute_relative_error",
                    "value": value,
                }
            )
        axes[1, 0].plot(
            coverages, values, marker="o", linewidth=2, color=COLORS[method], label=LABELS[method]
        )
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xticks(coverages, [f"{coverage}×" for coverage in coverages])
    axes[1, 0].set_xlabel("nominal read coverage")
    axes[1, 0].set_ylabel("mean absolute relative error")
    axes[1, 0].set_title("C  Conservative low-depth routing", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False, fontsize=8)

    outcomes = [sum(row["outcome"] == name for row in paired) for name in ("improved", "equal", "worse")]
    bar_position = -0.35
    bottom = 0.0
    for name, count, color in zip(
        ("improved", "equal", "worse"), outcomes, ("#59a14f", "#bab0ac", "#e15759")
    ):
        fraction = count / len(paired)
        axes[1, 1].bar(
            bar_position,
            fraction,
            bottom=bottom,
            width=0.58,
            color=color,
            label=f"{name}: {count}",
        )
        axes[1, 1].text(
            bar_position,
            bottom + fraction / 2,
            f"{fraction:.1%}",
            ha="center",
            va="center",
        )
        source_rows.append(
            {
                "panel": "D",
                "method": "depth_gated_controls",
                "stratum": name,
                "metric": "paired_family_fraction",
                "value": fraction,
            }
        )
        bottom += fraction
    axes[1, 1].set_xticks((bar_position,), ("Depth-gated vs total bases",))
    axes[1, 1].set_xlim(-1.0, 1.0)
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_ylabel("fraction of 1,485 paired family conditions")
    axes[1, 1].set_title("D  Family-level gains and regressions", loc="left", fontweight="bold")
    axes[1, 1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.66, 0.5))

    candidate_rows = by_method["depth_gated_controls"]
    selections = {
        (row["seed"], row["condition_id"]): (
            int(float(row["coverage"])), row["candidate_selected_method"]
        )
        for row in candidate_rows
    }
    if len(selections) != 27:
        raise ValueError("Expected 27 condition-level branch selections")
    bottom_values = np.zeros(len(coverages))
    for method in ("baseline_total_bases", "empirical_controls"):
        values = np.array(
            [
                sum(value == (coverage, method) for value in selections.values())
                for coverage in coverages
            ]
        )
        axes[2, 0].bar(
            np.arange(len(coverages)),
            values,
            bottom=bottom_values,
            color=COLORS[method],
            label=LABELS[method],
        )
        for coverage, value in zip(coverages, values):
            source_rows.append(
                {
                    "panel": "E",
                    "method": method,
                    "stratum": coverage,
                    "metric": "selected_condition_count",
                    "value": int(value),
                }
            )
        bottom_values += values
    axes[2, 0].set_xticks(np.arange(len(coverages)), [f"{coverage}×" for coverage in coverages])
    axes[2, 0].set_ylim(0, 10)
    axes[2, 0].set_xlabel("nominal read coverage")
    axes[2, 0].set_ylabel("condition count across three genomes")
    axes[2, 0].set_title("E  Both frozen branches are exercised", loc="left", fontweight="bold")
    axes[2, 0].legend(frameon=False, fontsize=8)

    marker_by_coverage = {1: "o", 5: "s", 20: "^"}
    for method in ("baseline_total_bases", "empirical_controls"):
        for coverage in coverages:
            rows = [
                row
                for row in executions
                if row["method"] == method and int(float(row["coverage"])) == coverage
            ]
            runtimes = [float(row["runtime_seconds"]) for row in rows]
            memories = [float(row["peak_rss_mib"]) for row in rows]
            axes[2, 1].scatter(
                runtimes,
                memories,
                color=COLORS[method],
                marker=marker_by_coverage[coverage],
                alpha=0.55,
                s=34,
            )
            axes[2, 1].scatter(
                statistics.median(runtimes),
                statistics.median(memories),
                color=COLORS[method],
                marker=marker_by_coverage[coverage],
                edgecolor="black",
                linewidth=0.8,
                s=90,
                label=f"{LABELS[method]}, {coverage}×",
            )
            for row in rows:
                execution_key = (
                    f"s{row['seed']}:{row['condition_id']}:{coverage}x"
                )
                source_rows.append(
                    {
                        "panel": "F",
                        "method": method,
                        "stratum": execution_key,
                        "metric": "runtime_seconds",
                        "value": float(row["runtime_seconds"]),
                    }
                )
                source_rows.append(
                    {
                        "panel": "F",
                        "method": method,
                        "stratum": execution_key,
                        "metric": "peak_rss_mib",
                        "value": float(row["peak_rss_mib"]),
                    }
                )
    axes[2, 1].set_xscale("log")
    axes[2, 1].set_xlabel("wall time per public command (s; log scale)")
    axes[2, 1].set_ylabel("direct-child peak RSS (MiB)")
    axes[2, 1].set_title("F  Descriptive time and memory", loc="left", fontweight="bold")
    axes[2, 1].legend(frameon=False, fontsize=7, ncol=2)

    fig.suptitle(
        "Frozen depth-gated quantification validation: passed gates with retained limitations",
        fontsize=15,
        fontweight="bold",
    )
    outdir.mkdir(parents=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"quantify_depth_gated_validation.{extension}", dpi=180)
    plt.close(fig)
    write_table(
        outdir / "panel_source.tsv",
        source_rows,
        ["panel", "method", "stratum", "metric", "value"],
    )
    (outdir / "figure_legend.md").write_text(
        "**Frozen depth-gated quantification validation.** A, aggregate mean "
        "absolute relative copy-number error across 1,485 known-catalogue family "
        "conditions per method; points are three independently simulated genomes. "
        "B, paired genome means for total-bases and the frozen depth-gated rule. "
        "C, coverage strata, including the conservative choice to retain the "
        "total-bases estimate at 1×. Ungated empirical controls had slightly lower "
        "aggregate error in this validation split. D, paired family outcomes; all "
        "329 regressions are retained. E, condition-level branch use across the "
        "three genomes. F, wall time and direct-child peak RSS for all 54 public "
        "commands; large outlined markers are group medians. These single "
        "executions are descriptive engineering measurements rather than repeated "
        "publication timing. The threshold and nine gates were committed and "
        "passed hosted CI before seeds 6401--6403 were generated. Controls were "
        "selected with simulation truth, catalogues were supplied, conditions "
        "within a genome are dependent, and the result is not biological truth.\n"
    )
    provenance = {
        "benchmark_id": "quantify_depth_gated_validation_v1",
        "inputs": {
            name: {"path": str(path), "sha256": digest_file(path)}
            for name, path in required.items()
        },
        "outputs": {
            path.name: digest_file(path)
            for path in sorted(outdir.iterdir())
            if path.name != "figure_provenance.json"
        },
        "warning": decision["warning"],
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
