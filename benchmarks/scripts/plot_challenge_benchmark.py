#!/usr/bin/env python3
"""Render a four-panel, source-backed challenge diagnostic figure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from benchmarks.challenge.schema import digest_file, read_table, write_table

TOOLS = ["tandemx", "trf", "tidehunter"]
LABELS = {"tandemx": "TX legacy", "trf": "TRF", "tidehunter": "TideHunter", "elastic": "TX elastic"}
SCENARIO_LABELS = {"indel_01pct": "Indels (0.1%)", "indel_1pct": "Indels (1%)", "indel_4pct": "Indels (4%)",
                   "substitution_1pct": "Substitutions (1%)", "substitution_5pct": "Substitutions (5%)"}


def make_figure(run: Path, output: Path, seed: int = 1101, elastic_run: Path | None = None) -> None:
    validation = json.loads((run / "validation.json").read_text())
    if not validation.get("complete"):
        raise ValueError("The benchmark must complete before rendering")
    rows = [r for r in read_table(run / "summary.tsv") if int(r["seed"]) == seed and "control" not in r["scenario"]]
    if not rows:
        raise ValueError("No positive-scenario rows for the requested seed")
    scenario_order = list(dict.fromkeys(r["scenario"] for r in read_table(run / "raw_runs.tsv")
                                       if int(r["seed"]) == seed and "control" not in r["scenario"]))
    tools = TOOLS + (["elastic"] if elastic_run else [])
    table_paths = [run / "summary.tsv", run / "raw_runs.tsv"]
    if elastic_run:
        receipt = json.loads((elastic_run / "validation.json").read_text())
        if not receipt.get("complete") or receipt.get("failed_runs"):
            raise ValueError("Elastic run must complete successfully")
        for scenario in scenario_order:
            relative = Path("datasets") / f"{scenario}_s{seed}" / "reads.fa"
            if digest_file(run / relative) != digest_file(elastic_run / relative):
                raise ValueError("Cannot compare different input reads in one development matrix")
        rows.extend(dict(r, tool="elastic") for r in read_table(elastic_run / "summary.tsv")
                    if int(r["seed"]) == seed and r["tool"] == "tandemx")
        table_paths.extend([elastic_run / "summary.tsv", elastic_run / "raw_runs.tsv"])
    lookup = {(r["scenario"], r["tool"]): r for r in rows}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 9.0), layout="constrained")
    color = LinearSegmentedColormap.from_list("tandemx_blue", ["#F5F7FA", "#A8C1E2", "#3B6FB6"])
    color.set_bad("#D9D9D9")
    source = []
    for ax, panel, field, title in zip(axes.flat, "abc", ["array_recall", "array_precision", "sequence_family_recall"],
                                     ["Array recall", "Array precision", "Sequence-supported family recovery"]):
        matrix = np.array([[float(lookup[(s, t)][field]) if lookup[(s, t)][field] != "NA" else np.nan
                            for t in tools] for s in scenario_order])
        # pcolormesh preserves vector-editable cells; imshow embeds a bitmap in SVG.
        ax.pcolormesh(np.arange(len(tools) + 1), np.arange(len(scenario_order) + 1), matrix,
                      vmin=0, vmax=1, cmap=color, shading="flat", rasterized=False)
        ax.set_xlim(0, len(tools))
        ax.set_ylim(len(scenario_order), 0)
        ax.set_xticks(np.arange(len(tools)) + 0.5, [LABELS[t] for t in tools])
        ax.set_yticks(np.arange(len(scenario_order)) + 0.5,
                      [SCENARIO_LABELS.get(s, s.replace("_", " ")) for s in scenario_order])
        ax.tick_params(length=0, pad=5)
        for y, scenario in enumerate(scenario_order):
            for x, tool in enumerate(tools):
                value = matrix[y, x]
                ax.text(x + 0.5, y + 0.5, "NA" if not np.isfinite(value) else f"{value:.2f}", ha="center", va="center",
                        color="white" if np.isfinite(value) and value > 0.8 else "#222222", fontsize=8)
                source.append({"panel": panel, "scenario": scenario, "tool": tool, "seed": seed,
                               "metric": field, "value": value, "source_dataset": lookup[(scenario, tool)]["dataset_id"]})
        ax.set_title(f"{panel}   {title}", loc="left", fontweight="bold", fontsize=10, pad=10)

    ax = axes[1, 1]
    dataset = f"indel_1pct_s{seed}"
    truth_path = run / "datasets" / dataset / "truth_arrays.tsv"
    truth = read_table(truth_path)
    matches = read_table(run / "runs" / dataset / "tandemx" / "rep1" / "matches.tsv")
    mismatch_ids = {r["read_id"] for r in matches if r["status"] == "unmatched"}
    read_id = next((r["read_id"] for r in truth if r["read_id"] in mismatch_ids), truth[0]["read_id"])
    colors = ["#444444", "#3B6FB6", "#BE8540", "#647A4D", "#7C5A97"]
    hatches = ["", "//", "..", "xx", "++"]
    interval_sources = [truth_path] + [run / "runs" / dataset / t / "rep1" / "predictions.tsv" for t in TOOLS]
    if elastic_run:
        interval_sources.append(elastic_run / "runs" / dataset / "tandemx" / "rep1" / "predictions.tsv")
    interval_rows = []
    for y, (label, path) in enumerate(zip(["Planted truth"] + [LABELS[t] for t in tools], interval_sources)):
        spans = [r for r in read_table(path) if r["read_id"] == read_id]
        ax.broken_barh([(int(r["start"]), int(r["end"]) - int(r["start"])) for r in spans], (len(tools) - y - 0.22, 0.44),
                       facecolors=colors[y], edgecolors="#333333", linewidths=0.5, hatch=hatches[y])
        interval_rows.extend({"track": label, "read_id": read_id, "start": int(r["start"]), "end": int(r["end"]),
                              "period": int(r["period"]), "source": str(path.resolve())} for r in spans)
    lengths = {r["read_id"]: int(r["length_bp"]) for r in read_table(run / "datasets" / dataset / "truth_reads.tsv")}
    ax.set(xlim=(0, lengths[read_id]), ylim=(-0.7, len(tools) + 0.7), xlabel="Position in read (bp; 0-based)")
    ax.set_yticks(list(range(len(tools), -1, -1)), ["Planted truth"] + [LABELS[t] for t in tools])
    ax.set_title(f"d   Array intervals: 1% indels\n{read_id}; seed {seed}", loc="left", fontweight="bold", fontsize=10)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    ax.set_axisbelow(True)
    repetitions = {tool: sorted({r["attempted_runs"] for r in rows if r["tool"] == tool}) for tool in tools}
    subtitle = (f"Matched inputs; development seed {seed}; baseline runs {','.join(repetitions['tandemx'])}, elastic {','.join(repetitions['elastic'])}"
                if elastic_run else "100 reads per scenario; 70 planted-positive reads; 3 timing repetitions")
    figure.suptitle("TandemX development challenge benchmark\n" + subtitle,
                     fontsize=11)
    output.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        figure.savefig(output / f"challenge_diagnostic.{suffix}", dpi=220)
    plt.close(figure)
    write_table(output / "heatmap_source.tsv", source, list(source[0]))
    write_table(output / "interval_source.tsv", interval_rows, list(interval_rows[0]))
    receipt = {"run_directory": str(run.resolve()), "seed": seed, "figure_type": "development_diagnostic_not_final_method_benchmark",
               "source_sha256": {str(p.resolve()): digest_file(p) for p in [*table_paths, *interval_sources]},
               "panel_count": 4, "svg_editable_text": True,
               "warning": "one_development_seed;repetitions_are_not_independent_datasets;negative_controls_not_shown_in_this_figure"}
    (output / "figure_provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1101)
    parser.add_argument("--elastic-run", type=Path)
    args = parser.parse_args()
    make_figure(args.run, args.outdir, args.seed, args.elastic_run)


if __name__ == "__main__":
    main()
