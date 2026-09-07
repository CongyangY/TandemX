#!/usr/bin/env python3
"""Render a six-panel TideCluster real-reference scaling audit."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, write_table


matplotlib.rcParams["svg.fonttype"] = "none"

STAGES = ("tidehunter", "clustering")
COLORS = {"tidehunter": "#4c78a8", "clustering": "#f28e2b"}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot(evidence: Path, outdir: Path) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    evidence = evidence.resolve()
    scaling_path = evidence / "scaling_summary.tsv"
    summary_path = evidence / "summary.json"
    manifest_path = evidence / "archive_manifest.json"
    if not all(path.is_file() for path in (scaling_path, summary_path, manifest_path)):
        raise ValueError("TideCluster scaling evidence is incomplete")
    summary = json.loads(summary_path.read_text())
    scaling = _rows(scaling_path)
    if (
        summary.get("complete") is not True
        or summary.get("accuracy")
        != "not_assessed_without_independent_real_array_and_family_truth"
        or len(scaling) < 4
    ):
        raise ValueError("TideCluster scaling evidence lacks its declared boundary")
    sample_ids = sorted(
        {row["sample_id"] for row in scaling},
        key=lambda sample_id: int(
            next(row["input_bases"] for row in scaling if row["sample_id"] == sample_id)
        ),
    )
    by_key = {(row["sample_id"], row["stage"]): row for row in scaling}
    if set(by_key) != {(sample_id, stage) for sample_id in sample_ids for stage in STAGES}:
        raise ValueError("Each TideCluster size must have exactly two stages")
    labels = [f"{float(by_key[(sample_id, STAGES[0])]['input_mb']):g} Mb" for sample_id in sample_ids]
    unique = {sample_id: by_key[(sample_id, STAGES[0])] for sample_id in sample_ids}
    normalized = {}
    for sample_id in sample_ids:
        path = evidence / "runs" / sample_id / "normalized/normalized_arrays.tsv"
        if not path.is_file():
            raise ValueError(f"Normalized arrays are missing for {sample_id}")
        normalized[sample_id] = _rows(path)
        if len(normalized[sample_id]) != int(unique[sample_id]["predicted_array_count"]):
            raise ValueError(f"Normalized array count differs for {sample_id}")

    outdir.mkdir(parents=True)
    source_rows: list[dict[str, object]] = []
    fig, axes = plt.subplots(3, 2, figsize=(15, 14), constrained_layout=True)
    x = np.arange(len(sample_ids))
    width = 0.34

    for offset, stage in zip((-width / 2, width / 2), STAGES):
        values = [float(by_key[(sample_id, stage)]["wall_seconds"]) for sample_id in sample_ids]
        axes[0, 0].bar(x + offset, values, width, color=COLORS[stage], label=stage)
        for sample_id, value in zip(sample_ids, values):
            source_rows.append(_source("A", sample_id, stage, "wall_seconds", value))
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("wall time (s)")
    axes[0, 0].set_title("A  End-to-end stage time", loc="left", fontweight="bold")
    axes[0, 0].legend(frameon=False)

    for offset, stage in zip((-width / 2, width / 2), STAGES):
        values = [
            float(by_key[(sample_id, stage)]["maximum_rss_kb"]) / 1024**2
            for sample_id in sample_ids
        ]
        axes[0, 1].bar(x + offset, values, width, color=COLORS[stage], label=stage)
        for sample_id, value in zip(sample_ids, values):
            source_rows.append(_source("B", sample_id, stage, "maximum_rss_gib", value))
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_ylabel("maximum RSS (GiB)")
    axes[0, 1].set_title("B  Container stage memory", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False)

    calls = ("predicted_array_count", "predicted_family_count")
    call_colors = ("#59a14f", "#b07aa1")
    for offset, metric, color in zip((-width / 2, width / 2), calls, call_colors):
        values = [int(unique[sample_id][metric]) for sample_id in sample_ids]
        axes[1, 0].bar(x + offset, values, width, color=color, label=metric.replace("predicted_", "").replace("_count", ""))
        for sample_id, value in zip(sample_ids, values):
            source_rows.append(_source("C", sample_id, "all", metric, value))
    axes[1, 0].set_xticks(x, labels)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylabel("calls (log scale)")
    axes[1, 0].set_title("C  Descriptive calls", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False)

    fractions = [float(unique[sample_id]["predicted_union_base_fraction"]) for sample_id in sample_ids]
    axes[1, 1].bar(x, np.array(fractions) * 100, color="#76b7b2", width=0.58)
    axes[1, 1].set_xticks(x, labels)
    axes[1, 1].set_ylabel("union coverage (%)")
    axes[1, 1].set_title("D  Sampled repeat fraction", loc="left", fontweight="bold")
    for sample_id, value in zip(sample_ids, fractions):
        source_rows.append(_source("D", sample_id, "all", "union_base_fraction", value))

    provenance = (
        ("exact_intermediate_interval_count", "exact", "#4c78a8"),
        ("clipped_interval_count", "clipped", "#f2cf5b"),
        ("merged_interval_count", "merged", "#e15759"),
    )
    bottoms = np.zeros(len(sample_ids))
    totals = np.array([int(unique[sample_id]["predicted_array_count"]) for sample_id in sample_ids])
    for metric, label, color in provenance:
        counts = np.array([int(unique[sample_id][metric]) for sample_id in sample_ids])
        values = counts / totals
        axes[2, 0].bar(x, values, bottom=bottoms, color=color, label=label)
        bottoms += values
        for sample_id, count, value in zip(sample_ids, counts, values):
            source_rows.append(_source("E", sample_id, label, "interval_fraction", value, int(count)))
    axes[2, 0].set_xticks(x, labels)
    axes[2, 0].set_ylim(0, 1)
    axes[2, 0].set_ylabel("fraction of final intervals")
    axes[2, 0].set_title("E  Coordinate provenance", loc="left", fontweight="bold")
    axes[2, 0].legend(frameon=False)

    period_values = [
        [int(row["period"]) for row in normalized[sample_id]] for sample_id in sample_ids
    ]
    axes[2, 1].boxplot(period_values, tick_labels=labels, showfliers=False)
    axes[2, 1].set_yscale("log")
    axes[2, 1].set_ylabel("selected representative period (bp, log scale)")
    axes[2, 1].set_title("F  Period distribution", loc="left", fontweight="bold")
    for sample_id, values in zip(sample_ids, period_values):
        for index, value in enumerate(values):
            source_rows.append(_source("F", sample_id, str(index + 1), "period_bp", value))

    fig.suptitle(
        "TideCluster 1.21.2 on nested MorexV3 reference windows\n"
        "Descriptive output; independent biological accuracy truth is unavailable",
        fontsize=14,
        fontweight="bold",
    )
    svg = outdir / "tidecluster_morex_reference_scaling.svg"
    pdf = outdir / "tidecluster_morex_reference_scaling.pdf"
    png = outdir / "tidecluster_morex_reference_scaling.png"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    write_table(outdir / "panel_source.tsv", source_rows, list(source_rows[0]))
    svg_text = svg.read_text(encoding="utf-8")
    receipt = {
        "complete": True,
        "panel_count": 6,
        "input_sha256": {
            str(path.relative_to(evidence)): digest_file(path)
            for path in (
                scaling_path,
                summary_path,
                manifest_path,
                *(
                    evidence / "runs" / sample_id / "normalized/normalized_arrays.tsv"
                    for sample_id in sample_ids
                ),
            )
        },
        "output_sha256": {
            path.name: digest_file(path) for path in (svg, pdf, png, outdir / "panel_source.tsv")
        },
        "svg_text_nodes": svg_text.count("<text"),
        "svg_image_nodes": svg_text.count("<image"),
        "warning": "descriptive_nested_reference_windows_without_independent_accuracy_truth",
    }
    (outdir / "figure_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def _source(
    panel: str,
    sample_id: str,
    category: str,
    metric: str,
    value: float | int,
    count: int | str = "NA",
) -> dict[str, object]:
    return {
        "panel": panel,
        "sample_id": sample_id,
        "category": category,
        "metric": metric,
        "value": value,
        "count": count,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.evidence, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
