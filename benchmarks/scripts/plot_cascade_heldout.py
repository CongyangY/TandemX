"""Render a six-panel audit of the frozen cascade held-out benchmark."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

matplotlib.rcParams["svg.fonttype"] = "none"

from benchmarks.challenge.schema import digest_file, write_table


TOOLS = ("tandemx", "trf", "tidehunter")
TOOL_LABELS = {"tandemx": "TandemX", "trf": "TRF", "tidehunter": "TideHunter"}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _heatmap(
    axis: plt.Axes,
    values: np.ndarray,
    rows: list[str],
    columns: list[str],
    title: str,
    *,
    minimum: float = 0.9,
) -> None:
    masked = np.ma.masked_invalid(values)
    axis.pcolormesh(
        np.arange(values.shape[1] + 1),
        np.arange(values.shape[0] + 1),
        masked,
        vmin=minimum,
        vmax=1.0,
        cmap="viridis",
        shading="flat",
    )
    axis.set_xlim(0, values.shape[1])
    axis.set_ylim(values.shape[0], 0)
    axis.set_title(title, loc="left", fontweight="bold")
    axis.set_yticks(np.arange(len(rows)) + 0.5, rows)
    axis.set_xticks(np.arange(len(columns)) + 0.5, columns, rotation=55, ha="right")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            text = "NA" if not math.isfinite(value) else f"{value:.2f}"
            color = "black" if math.isfinite(value) and value > 0.96 else "white"
            axis.text(j + 0.5, i + 0.5, text, ha="center", va="center", fontsize=6.5, color=color)


def plot(evidence: Path, outdir: Path) -> dict[str, object]:
    evidence = evidence.resolve()
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    required = {
        "summary": evidence / "run" / "summary.tsv",
        "raw": evidence / "run" / "raw_runs.tsv",
        "paired": evidence / "evaluation" / "paired_tidehunter.tsv",
        "gates": evidence / "evaluation" / "gate_results.json",
        "config": evidence / "frozen_config.yaml",
        "archive": evidence / "archive_summary.json",
    }
    if any(not path.is_file() for path in required.values()):
        raise ValueError("Held-out evidence is incomplete")
    summary = _read(required["summary"])
    raw = _read(required["raw"])
    paired = _read(required["paired"])
    gates = json.loads(required["gates"].read_text())
    config = yaml.safe_load(required["config"].read_text())
    if gates["matrix"] != {"raw_rows": len(raw), "summary_rows": len(summary), "datasets": len(paired)}:
        raise ValueError("Figure inputs disagree with the gate matrix")
    scenarios = [row["name"] for row in config["scenarios"]]
    positive = [
        row["name"]
        for row in config["scenarios"]
        if float(row.get("positive_fraction", 1.0)) > 0
    ]
    negative = [name for name in scenarios if name not in positive]
    by_group: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in summary:
        by_group.setdefault((row["scenario"], row["tool"]), []).append(row)

    source_rows: list[dict[str, object]] = []
    matrices: dict[str, np.ndarray] = {}
    for panel, field in (("A", "array_recall"), ("B", "array_precision")):
        matrix = np.full((len(TOOLS), len(positive)), np.nan)
        for i, tool in enumerate(TOOLS):
            for j, scenario in enumerate(positive):
                values = [
                    value
                    for row in by_group[(scenario, tool)]
                    for value in [_number(row.get(field))]
                    if value is not None
                ]
                if values:
                    matrix[i, j] = _mean(values)
                source_rows.append(
                    {
                        "panel": panel,
                        "scenario": scenario,
                        "tool": tool,
                        "metric": field,
                        "value": matrix[i, j] if values else None,
                        "seed": "mean_of_3101_3102_3103",
                        "status": "observed" if values else "missing",
                    }
                )
        matrices[panel] = matrix

    negative_matrix = np.full((len(TOOLS), len(negative)), np.nan)
    for i, tool in enumerate(TOOLS):
        for j, scenario in enumerate(negative):
            values = [
                value
                for row in by_group[(scenario, tool)]
                for value in [_number(row.get("negative_read_call_rate"))]
                if value is not None
            ]
            if values:
                negative_matrix[i, j] = _mean(values)
            source_rows.append(
                {
                    "panel": "C",
                    "scenario": scenario,
                    "tool": tool,
                    "metric": "negative_read_call_rate",
                    "value": negative_matrix[i, j] if values else None,
                    "seed": "mean_of_3101_3102_3103",
                    "status": "observed" if values else "all_repetitions_failed",
                }
            )

    fig, axes = plt.subplots(3, 2, figsize=(16, 15), constrained_layout=True)
    _heatmap(
        axes[0, 0],
        matrices["A"],
        [TOOL_LABELS[tool] for tool in TOOLS],
        positive,
        "A  Mean array recall across held-out seeds",
    )
    _heatmap(
        axes[0, 1],
        matrices["B"],
        [TOOL_LABELS[tool] for tool in TOOLS],
        positive,
        "B  Mean raw-call array precision",
        minimum=0.5,
    )
    _heatmap(
        axes[1, 0],
        negative_matrix,
        [TOOL_LABELS[tool] for tool in TOOLS],
        negative,
        "C  Negative-read call rate (NA = process failure)",
        minimum=0.0,
    )

    runtime = axes[1, 1]
    rss = axes[2, 0]
    for panel, axis, field, threshold, title, ylabel in (
        (
            "D",
            runtime,
            "median_runtime_seconds_ratio",
            float(config["acceptance_gates"]["tidehunter_runtime_geometric_mean_ratio_max"]),
            "D  Paired TandemX/TideHunter wall-time ratio",
            "runtime ratio",
        ),
        (
            "E",
            rss,
            "median_peak_rss_mib_ratio",
            float(config["acceptance_gates"]["tidehunter_peak_rss_geometric_mean_ratio_max"]),
            "E  Paired TandemX/TideHunter peak-RSS ratio",
            "peak RSS ratio",
        ),
    ):
        values: list[float] = []
        for index, row in enumerate(paired):
            value = _number(row.get(field))
            if value is None:
                continue
            values.append(value)
            is_positive = row["positive"].lower() == "true"
            axis.scatter(index, value, s=22, color="#2c7fb8" if is_positive else "#f28e2b", alpha=0.85)
            source_rows.append(
                {
                    "panel": panel,
                    "scenario": row["scenario"],
                    "tool": "tandemx_over_tidehunter",
                    "metric": field,
                    "value": value,
                    "seed": row["seed"],
                    "status": "positive" if is_positive else "negative_control",
                }
            )
        geomean = math.exp(sum(math.log(value) for value in values) / len(values))
        axis.scatter([], [], s=22, color="#2c7fb8", label="positive")
        axis.scatter([], [], s=22, color="#f28e2b", label="negative control")
        axis.axhline(threshold, color="#d62728", linestyle="--", linewidth=1.2, label=f"gate {threshold:g}")
        axis.axhline(geomean, color="black", linestyle=":", linewidth=1.2, label=f"geomean {geomean:.3f}")
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel("scenario-seed pair (48 total)")
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(frameon=False, fontsize=8)

    gate_axis = axes[2, 1]
    gate_rows = gates["gates"]
    positions = np.arange(len(gate_rows))
    colors = ["#2ca02c" if row["passed"] else "#d62728" for row in gate_rows]
    gate_axis.barh(positions, [1] * len(gate_rows), color=colors, alpha=0.85)
    gate_axis.set_yticks(positions, [row["name"].replace("_", " ") for row in gate_rows], fontsize=8)
    gate_axis.set_xlim(0, 1.02)
    gate_axis.set_xticks([])
    gate_axis.invert_yaxis()
    gate_axis.set_title("F  Predeclared promotion gates", loc="left", fontweight="bold")
    for index, row in enumerate(gate_rows):
        observed = row.get("observed")
        observed_text = "NA" if observed is None else f"{float(observed):.3g}"
        label = f"{'PASS' if row['passed'] else 'FAIL'}  {observed_text} {row['operator']} {float(row['threshold']):g}"
        gate_axis.text(0.02, index, label, va="center", color="white", fontsize=7.5, fontweight="bold")
        source_rows.append(
            {
                "panel": "F",
                "scenario": row["name"],
                "tool": "gate",
                "metric": "passed",
                "value": 1 if row["passed"] else 0,
                "seed": "heldout_matrix",
                "status": label,
            }
        )
    fig.suptitle(
        "Cascade discovery held-out audit: accuracy and memory pass; promotion fails on runtime and comparator completion",
        fontsize=15,
        fontweight="bold",
    )
    outdir.mkdir(parents=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"cascade_heldout.{extension}", dpi=180)
    plt.close(fig)
    write_table(
        outdir / "panel_source.tsv",
        source_rows,
        ["panel", "scenario", "tool", "metric", "value", "seed", "status"],
    )
    (outdir / "figure_legend.md").write_text(
        "**Cascade held-out audit.** A-B, mean one-to-one array recall and raw-call "
        "precision across independent seeds 3101-3103 for 13 positive scenarios. "
        "C, negative-read call rates for three controls; TRF is NA on the low-"
        "complexity control because all nine processes reached the frozen 180-s "
        "timeout. D-E, paired TandemX/TideHunter median wall-time and direct-child "
        "peak-RSS ratios for all 48 scenario-seed groups; dotted lines are geometric "
        "means and dashed lines are predeclared limits. F, all 12 frozen gates. "
        "TandemX passed accuracy, determinism, false-call, family and RSS gates but "
        "failed promotion because the runtime ratio was 2.458 (>2.0) and comparator "
        "completion was required. Technical repetitions measure execution variation, "
        "not biological replication.\n"
    )
    provenance = {
        "benchmark_id": gates["benchmark_id"],
        "gate_status": gates["status"],
        "failed_gate_names": gates["failed_gate_names"],
        "inputs": {name: {"path": str(path), "sha256": digest_file(path)} for name, path in required.items()},
        "outputs": {
            path.name: digest_file(path)
            for path in sorted(outdir.iterdir())
            if path.name != "figure_provenance.json"
        },
        "warning": "heldout_failure_retained;mean_cells_average_three_independent_seeds;timing_repetitions_not_biological_replicates",
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
