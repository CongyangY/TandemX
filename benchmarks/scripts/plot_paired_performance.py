"""Plot all same-tool performance pairs, retaining output parity and RSS trade-offs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmarks.challenge.schema import digest_file, read_table, write_table


def plot(source: Path, outdir: Path) -> None:
    validation = json.loads((source / "validation.json").read_text())
    rows = read_table(source / "summary.tsv")
    if not validation["all_pairs_pass"] or not rows or any(
        r["all_successful"] != "True" or r["all_six_outputs_identical"] != "True" for r in rows
    ):
        raise ValueError("Output equality and successful execution are required")
    outdir.mkdir(parents=True, exist_ok=False)
    labels = [r["dataset"].removesuffix("_s2101").replace("_", " ") for r in rows]
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 8,
                                "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")
    colours = ("#756BB1", "#158E87")
    sources = []
    for panel, metric, xlabel in (("A", "runtime_seconds", "Median wall time (s)"),
                                  ("D", "cpu_user_seconds", "Median user CPU time (s)")):
        ax = axes[0, 0] if panel == "A" else axes[1, 1]
        for i, row in enumerate(rows):
            values = [float(row[f"{v}_median_{metric}"]) for v in ("baseline", "native_seed")]
            ax.plot(values, [i, i], color="#BBBBBB", linewidth=1)
            for j, value in enumerate(values):
                ax.scatter(value, i, color=colours[j], s=20,
                           label=("Python seeds" if j == 0 else "Native seeds") if i == 0 else None)
                sources.append(dict(panel=panel, dataset=row["dataset"], metric=metric,
                                    variant=("baseline", "native_seed")[j], value=value))
        ax.set_title(f"{panel}   Identical outputs, lower computation time", loc="left", weight="bold")
        ax.set_xlabel(xlabel)
        ax.legend(frameon=False, fontsize=7)
    for panel, ax, field, xlabel, ref in (
        ("B", axes[0, 1], "speedup_baseline_over_native", "Wall-time speedup (baseline / native)", 1),
        ("C", axes[1, 0], "rss_ratio_native_over_baseline", "Peak RSS change (%)", 0),
    ):
        for i, row in enumerate(rows):
            value = float(row[field])
            if panel == "C":
                value = (value - 1) * 100
            ax.barh(i, value, color="#CB614D" if panel == "C" and value > 0 else colours[1], height=.65)
            sources.append(dict(panel=panel, dataset=row["dataset"], metric=field,
                                variant="ratio" if panel == "B" else "percent_change", value=value))
        ax.axvline(ref, color="#666666", linestyle="--", linewidth=.8)
        ax.set_xlabel(xlabel)
        ax.set_title("B   Every scenario, including negative controls" if panel == "B" else
                     "C   Memory does not improve in every scenario", loc="left", weight="bold")
    for ax in axes.flat:
        ax.set(yticks=range(len(rows)), yticklabels=labels)
        ax.invert_yaxis()
    fig.suptitle("Native seed processing: paired development experiment", fontsize=14, weight="bold")
    fig.supxlabel("16 datasets × 2 implementations × 3 repeats; 100 reads per dataset, seed 2101; six outputs byte-identical\n"
                  "Same-tool optimization only. Technical repeats are not biological replicates.", fontsize=8)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"paired_performance.{ext}", dpi=180)
    plt.close(fig)
    write_table(outdir / "panel_source.tsv", sources, ["panel", "dataset", "metric", "variant", "value"])
    speed = [float(r["speedup_baseline_over_native"]) for r in rows]
    rss = [float(r["rss_ratio_native_over_baseline"]) for r in rows]
    receipt = {"input_sha256": {name: digest_file(source / name) for name in
                               ("summary.tsv", "raw_runs.tsv", "validation.json", "environment.json")},
               "script_sha256": digest_file(Path(__file__)),
               "speedup_range": [min(speed), max(speed)], "median_dataset_speedup": statistics.median(speed),
               "rss_change_percent_range": [(min(rss) - 1) * 100, (max(rss) - 1) * 100],
               "rss_increased_datasets": sum(r > 1 for r in rss),
               "outputs": {p.name: digest_file(p) for p in outdir.iterdir()}}
    (outdir / "figure_provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.outdir)
