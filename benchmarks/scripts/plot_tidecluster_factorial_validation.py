#!/usr/bin/env python3
"""Render the frozen six-panel TideCluster factorial-validation figure."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


BLUE = "#2A6FBB"
ORANGE = "#D97706"
CHARCOAL = "#30343B"
GREY = "#9299A1"
GRID = "#D9DDE2"
SETTING_STYLE = {
    "default_primary": {"color": BLUE, "label": "Default"},
    "matched_period_sensitivity": {
        "color": ORANGE,
        "label": "Matched P30–1000",
    },
}
SOURCE_FIELDS = (
    "panel",
    "seed",
    "setting",
    "metric",
    "stage",
    "value",
    "unit",
    "status",
    "warning",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def optional_float(value: str | None) -> float | None:
    if value is None or value.strip().lower() in {"", "none", "na", "nan"}:
        return None
    return float(value)


def style_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(CHARCOAL)
    axis.spines["bottom"].set_color(CHARCOAL)
    axis.tick_params(colors=CHARCOAL, labelsize=8)
    axis.grid(True, color=GRID, linewidth=0.55, alpha=0.7, zorder=0)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.13,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=12,
        fontweight="bold",
        color=CHARCOAL,
        va="top",
    )


def plot_metric_pair(
    axis: plt.Axes,
    rows: list[dict[str, str]],
    panel: str,
    metrics: tuple[tuple[str, str, str, str], ...],
    source_rows: list[dict[str, Any]],
) -> None:
    seeds = sorted({int(row["seed"]) for row in rows})
    for setting, style in SETTING_STYLE.items():
        selected = {int(row["seed"]): row for row in rows if row["setting"] == setting}
        for field, label, marker, line_style in metrics:
            values = [optional_float(selected[seed].get(field)) for seed in seeds]
            valid = [(seed, value) for seed, value in zip(seeds, values) if value is not None]
            if valid:
                axis.plot(
                    [seed for seed, _value in valid],
                    [value for _seed, value in valid],
                    color=style["color"],
                    marker=marker,
                    markerfacecolor=("white" if line_style != "-" else style["color"]),
                    markeredgecolor=CHARCOAL,
                    linewidth=1.25,
                    linestyle=line_style,
                    markersize=5,
                    label=f"{style['label']}: {label}",
                    zorder=3,
                )
            for seed, value in zip(seeds, values):
                source_rows.append(
                    {
                        "panel": panel,
                        "seed": seed,
                        "setting": setting,
                        "metric": field,
                        "value": "" if value is None else value,
                        "status": "unavailable" if value is None else "ok",
                        "warning": "technical_seed_not_biological_replicate",
                    }
                )
    axis.set_xticks(seeds, [str(seed) for seed in seeds])
    axis.set_xlabel("Frozen simulation seed")


def plot(result_dir: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite figure directory: {outdir}")
    paths = {
        "summary": result_dir / "summary.tsv",
        "run_receipt": result_dir / "run_receipt.json",
        "independent_verification": result_dir / "independent_verification.json",
        "environment": result_dir / "environment.json",
    }
    for name in ("cell_fates", "stage_fates"):
        path = result_dir / f"{name}.tsv"
        if path.is_file():
            paths[name] = path
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"TideCluster figure inputs are incomplete: {missing}")
    receipt = json.loads(paths["run_receipt"].read_text(encoding="utf-8"))
    verification = json.loads(
        paths["independent_verification"].read_text(encoding="utf-8")
    )
    if receipt.get("complete") is not True:
        raise ValueError("TideCluster factorial run is incomplete")
    if verification.get("verification_passed") is not True or verification.get("failures") != []:
        raise ValueError("TideCluster factorial accuracy failed independent verification")
    rows = read_tsv(paths["summary"])
    expected = {
        (seed, setting)
        for seed in (6401, 6402, 6403)
        for setting in SETTING_STYLE
    }
    if {(int(row["seed"]), row["setting"]) for row in rows} != expected:
        raise ValueError("TideCluster figure requires all six frozen runs")
    allowed_statuses = {
        "ok",
        "external_resource_failure",
        "external_resource_failure_parent_v1",
        "normalization_or_evaluation_failure",
    }
    unexpected = sorted({row["status"] for row in rows} - allowed_statuses)
    if unexpected:
        raise ValueError(f"unexpected TideCluster run statuses: {unexpected}")
    success_count = sum(row["status"] == "ok" for row in rows)
    cell_fates = read_tsv(paths["cell_fates"]) if "cell_fates" in paths else []
    fate_by_key = {
        (int(row["seed"]), row["setting"]): row for row in cell_fates
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "text.color": CHARCOAL,
            "axes.labelcolor": CHARCOAL,
            "axes.titlecolor": CHARCOAL,
        }
    )
    figure, axes = plt.subplots(3, 2, figsize=(8.3, 10.4), constrained_layout=False)
    figure.subplots_adjust(
        left=0.10, right=0.97, top=0.855, bottom=0.165, hspace=0.58, wspace=0.31
    )
    figure.suptitle(
        "TideCluster factorial assembly validation",
        x=0.10,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    figure.text(
        0.10,
        0.945,
        "Pinned v1.21.2; three frozen 10-Mb planted-truth assemblies; "
        f"accuracy available for {success_count}/6 cells",
        fontsize=8.5,
        ha="left",
    )
    source_rows: list[dict[str, Any]] = []

    def add_accuracy_availability(axis: plt.Axes) -> None:
        axis.text(
            0.98,
            0.04,
            f"Evaluated: {success_count}/6\nFailed cells: {6 - success_count} (NA)",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.7,
            color=CHARCOAL,
        )

    axis = axes[0, 0]
    plot_metric_pair(
        axis,
        rows,
        "A",
        (
            ("array_recall", "Recall", "o", "-"),
            ("array_precision", "Precision", "s", "--"),
        ),
        source_rows,
    )
    axis.set_ylim(-0.02, 1.04)
    axis.set_ylabel("Fraction")
    axis.set_title("Array calls recover planted intervals")
    style_axis(axis)
    add_accuracy_availability(axis)
    panel_label(axis, "A")

    axis = axes[0, 1]
    plot_metric_pair(
        axis,
        rows,
        "B",
        (
            ("base_union_recall", "Recall", "o", "-"),
            ("base_union_precision", "Precision", "s", "--"),
        ),
        source_rows,
    )
    axis.set_ylim(-0.02, 1.04)
    axis.set_ylabel("Fraction")
    axis.set_title("Base-union overlap measures coverage")
    style_axis(axis)
    add_accuracy_availability(axis)
    panel_label(axis, "B")

    axis = axes[1, 0]
    plot_metric_pair(
        axis,
        rows,
        "C",
        (
            ("matched_boundary_mae_bp", "Boundary MAE", "o", "-"),
            ("matched_period_mae_bp", "Period MAE", "s", "--"),
        ),
        source_rows,
    )
    axis.set_yscale("symlog", linthresh=0.5)
    axis.set_ylim(0, 20)
    axis.set_ylabel("Mean absolute error (bp)")
    axis.set_title("Matched-call errors remain conditional")
    style_axis(axis)
    add_accuracy_availability(axis)
    panel_label(axis, "C")

    axis = axes[1, 1]
    plot_metric_pair(
        axis,
        rows,
        "D",
        (
            ("cyclic_monomer_recall", "Monomer recall", "o", "-"),
            (
                "homologous_consensus_fraction",
                "Homologous consensus fraction",
                "s",
                "--",
            ),
        ),
        source_rows,
    )
    axis.set_ylim(-0.02, 1.04)
    axis.set_ylabel("Fraction")
    axis.set_title("Sequence-supported family recovery")
    style_axis(axis)
    add_accuracy_availability(axis)
    panel_label(axis, "D")

    for panel, axis, suffix, title, ylabel, unit in (
        ("E", axes[2, 0], "wall_seconds", "Stage elapsed time", "Wall time (s)", "seconds"),
        (
            "F",
            axes[2, 1],
            "maximum_rss_kb",
            "Stage memory",
            "Maximum RSS (MiB)",
            "MiB",
        ),
    ):
        seeds = sorted({int(row["seed"]) for row in rows})
        for setting, style in SETTING_STYLE.items():
            selected = {
                int(row["seed"]): row for row in rows if row["setting"] == setting
            }
            for stage, marker, line_style in (
                ("tidehunter", "o", "-"),
                ("clustering", "^", ":"),
            ):
                field = f"{stage}_{suffix}"
                values = [optional_float(selected[seed].get(field)) for seed in seeds]
                if suffix == "maximum_rss_kb":
                    values = [None if value is None else value / 1024 for value in values]
                stage_states = [
                    fate_by_key.get((seed, setting), {}).get(
                        f"{stage}_status",
                        "ok" if selected[seed]["status"] == "ok" else "unavailable",
                    )
                    for seed in seeds
                ]
                valid = [
                    (seed, value)
                    for seed, value, state in zip(seeds, values, stage_states)
                    if value is not None and state == "ok"
                ]
                failed = [
                    (seed, value)
                    for seed, value, state in zip(seeds, values, stage_states)
                    if value is not None and "fail" in state
                ]
                if valid:
                    axis.plot(
                        [seed for seed, _value in valid],
                        [value for _seed, value in valid],
                        color=style["color"],
                        marker=marker,
                        markerfacecolor=("white" if stage == "clustering" else style["color"]),
                        markeredgecolor=CHARCOAL,
                        linewidth=1.25,
                        linestyle=line_style,
                        markersize=5,
                        label=f"{style['label']}: {stage}",
                        zorder=3,
                    )
                if failed:
                    axis.scatter(
                        [seed for seed, _value in failed],
                        [value for _seed, value in failed],
                        color=style["color"],
                        marker="x",
                        linewidths=1.5,
                        s=34,
                        label=f"{style['label']}: failed {stage} resource",
                        zorder=4,
                    )
                for seed, value, state in zip(seeds, values, stage_states):
                    source_rows.append(
                        {
                            "panel": panel,
                            "seed": seed,
                            "setting": setting,
                            "metric": field,
                            "stage": stage,
                            "value": "" if value is None else value,
                            "unit": unit,
                            "status": (
                                "unavailable"
                                if value is None
                                else "failed_resource"
                                if "fail" in state
                                else "ok"
                            ),
                            "warning": (
                                "failed_stage_accuracy_is_NA_not_zero"
                                if "fail" in state
                                else "single_same_host_execution_per_seed_setting"
                            ),
                        }
                    )
        axis.set_xticks(seeds, [str(seed) for seed in seeds])
        axis.set_xlabel("Frozen simulation seed")
        axis.set_yscale("log")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        style_axis(axis)
        panel_label(axis, panel)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        frameon=False,
        fontsize=7.1,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.905),
    )
    handles, labels = axes[2, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        frameon=False,
        fontsize=6.7,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.055),
    )
    figure.text(
        0.10,
        0.012,
        "Same-process simulations and technical seeds; not biological replicates. "
        "Failed stages remain missing rather than zero.",
        fontsize=7.1,
        ha="left",
    )

    outdir.mkdir(parents=True)
    svg = outdir / "tidecluster_factorial_validation.svg"
    pdf = outdir / "tidecluster_factorial_validation.pdf"
    png = outdir / "tidecluster_factorial_validation.png"
    figure.savefig(svg, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=240, bbox_inches="tight")
    plt.close(figure)
    source_path = outdir / "panel_source.tsv"
    with source_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SOURCE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(source_rows)
    svg_text = svg.read_text(encoding="utf-8")
    provenance = {
        "schema_version": 1,
        "complete": True,
        "panel_count": 6,
        "frozen_run_count": len(rows),
        "successful_accuracy_run_count": success_count,
        "unavailable_accuracy_run_count": len(rows) - success_count,
        "panel_source_rows": len(source_rows),
        "svg_text_element_count": svg_text.count("<text"),
        "svg_raster_image_element_count": svg_text.count("<image"),
        "source_files": {
            name: {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}
            for name, path in paths.items()
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": digest(path)}
            for path in (svg, pdf, png, source_path)
        },
        "warning": (
            "same_process_simulated_genomes;technical_seeds_not_biological_replicates;"
            "conditional_errors_on_matched_calls"
        ),
    }
    (outdir / "figure_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.result_dir, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
