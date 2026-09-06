"""Render source-backed one-thread real-read comparator diagnostics."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from benchmarks.challenge.schema import digest_file, iter_table, write_table


TOOLS = ("tandemx", "trf", "tidehunter")
COLORS = {"tandemx": "#4477AA", "trf": "#CC6677", "tidehunter": "#228833"}
DISPLAY = {"tandemx": "TandemX", "trf": "TRF", "tidehunter": "TideHunter"}
MARKERS = ("o", "s", "D", "^", "v", "P", "X", "<", ">")


@dataclass(frozen=True)
class DiagnosticRun:
    run_id: str
    material: str
    path: Path
    read_count: int
    total_bases: int
    rows: tuple[dict, ...]


def parse_run(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("run must be MATERIAL=PATH")
    material, raw_path = value.split("=", 1)
    if not material.strip() or not raw_path:
        raise argparse.ArgumentTypeError("run must be MATERIAL=PATH")
    return material.strip(), Path(raw_path)


def _validate_archive(path: Path) -> None:
    manifest_path = path / "archive_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"Missing archive manifest: {path}")
    manifest = json.loads(manifest_path.read_text())
    indexed = {row.get("file"): row for row in manifest}
    for name in ("environment.json", "summary.tsv"):
        source = path / name
        row = indexed.get(name)
        if (
            not source.is_file()
            or not row
            or row.get("sha256") != digest_file(source)
            or row.get("bytes") != source.stat().st_size
        ):
            raise ValueError(f"Archive hash/size differs for {source}")


def load_run(material: str, path: Path) -> DiagnosticRun:
    path = path.resolve()
    _validate_archive(path)
    environment = json.loads((path / "environment.json").read_text())
    if (
        environment.get("accuracy") != "not_assessed_without_curated_independent_truth"
        or int(environment.get("threads", 1)) != 1
        or int(environment.get("repetitions", 1)) != 1
        or "not publication ranking" not in environment.get("resources", "")
    ):
        raise ValueError(f"Unexpected evidence boundary or execution design: {path}")
    sample = environment.get("input", {})
    read_count = int(sample.get("read_count", 0))
    total_bases = int(sample.get("total_bases", 0))
    if read_count <= 0 or total_bases <= 0:
        raise ValueError(f"Invalid input denominator: {path}")

    rows = tuple(iter_table(path / "summary.tsv", {
        "tool", "exit_code", "runtime_seconds", "peak_rss_mib", "timed_out",
        "observed_in_scope_calls", "observed_positive_reads", "observed_union_bp",
        "observed_union_base_fraction", "normalization",
    }))
    if len(rows) != len(TOOLS) or {row["tool"] for row in rows} != set(TOOLS):
        raise ValueError(f"Require exactly the three declared tools: {path}")
    parsed: list[dict] = []
    for row in rows:
        if row["exit_code"] != "0" or row["timed_out"].lower() != "false" or row["normalization"] != "ok":
            raise ValueError(f"Incomplete tool result: {path}")
        numeric = {
            "tool": row["tool"],
            "runtime_seconds": float(row["runtime_seconds"]),
            "peak_rss_mib": float(row["peak_rss_mib"]),
            "observed_in_scope_calls": int(row["observed_in_scope_calls"]),
            "observed_positive_reads": int(row["observed_positive_reads"]),
            "observed_union_bp": int(row["observed_union_bp"]),
            "observed_union_base_fraction": float(row["observed_union_base_fraction"]),
        }
        if (
            not math.isfinite(numeric["runtime_seconds"])
            or not math.isfinite(numeric["peak_rss_mib"])
            or numeric["runtime_seconds"] <= 0
            or numeric["peak_rss_mib"] <= 0
            or not 0 <= numeric["observed_positive_reads"] <= read_count
            or not 0 <= numeric["observed_union_bp"] <= total_bases
            or abs(numeric["observed_union_bp"] / total_bases - numeric["observed_union_base_fraction"]) > 1e-12
        ):
            raise ValueError(f"Invalid summary metric: {path}")
        parsed.append(numeric)
    return DiagnosticRun(path.name, material, path, read_count, total_bases, tuple(parsed))


def source_rows(runs: list[DiagnosticRun]) -> list[dict]:
    output: list[dict] = []
    for run in runs:
        by_tool = {row["tool"]: row for row in run.rows}
        trf = by_tool["trf"]
        for tool in TOOLS:
            row = by_tool[tool]
            output.append({
                "run_id": run.run_id,
                "material": run.material,
                "tool": tool,
                "read_count": run.read_count,
                "total_bases": run.total_bases,
                "input_Gb": run.total_bases / 1e9,
                "runtime_seconds": row["runtime_seconds"],
                "peak_rss_mib": row["peak_rss_mib"],
                "throughput_Mbp_per_second": run.total_bases / 1e6 / row["runtime_seconds"],
                "runtime_ratio_to_trf": row["runtime_seconds"] / trf["runtime_seconds"],
                "rss_ratio_to_trf": row["peak_rss_mib"] / trf["peak_rss_mib"],
                "observed_in_scope_calls": row["observed_in_scope_calls"],
                "observed_call_density_per_Mbp": row["observed_in_scope_calls"] / (run.total_bases / 1e6),
                "observed_positive_read_fraction": row["observed_positive_reads"] / run.read_count,
                "observed_union_base_fraction": row["observed_union_base_fraction"],
            })
    return output


def _scatter_panel(axis, rows: list[dict], y_field: str, ylabel: str, panel: str, *, log_y: bool = False) -> None:
    materials = sorted({row["material"] for row in rows})
    markers = {material: MARKERS[index % len(MARKERS)] for index, material in enumerate(materials)}
    for tool in TOOLS:
        for material in materials:
            selected = sorted(
                (row for row in rows if row["tool"] == tool and row["material"] == material),
                key=lambda row: row["input_Gb"],
            )
            if not selected:
                continue
            x = [row["input_Gb"] for row in selected]
            y = [row[y_field] for row in selected]
            if len(selected) > 1:
                axis.plot(x, y, color=COLORS[tool], alpha=0.22, lw=0.8)
            axis.scatter(
                x, y, color=COLORS[tool], marker=markers[material], s=38,
                edgecolor="black", linewidth=0.25, alpha=0.9,
            )
    axis.set_xscale("log")
    if log_y:
        axis.set_yscale("log")
    axis.set_xlabel("Input sequence (Gb)")
    axis.set_ylabel(ylabel)
    axis.set_title(panel, loc="left", fontweight="bold")
    axis.grid(alpha=0.16, which="both")


def plot(run_specs: list[tuple[str, Path]], outdir: Path) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    if len(run_specs) < 2:
        raise ValueError("Require at least two completed run archives")
    runs = sorted(
        (load_run(material, path) for material, path in run_specs),
        key=lambda run: (run.material, run.total_bases, run.run_id),
    )
    if len({run.run_id for run in runs}) != len(runs):
        raise ValueError("Run directory names must be unique")
    rows = source_rows(runs)

    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "DejaVu Sans",
        "font.size": 8.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig = plt.figure(figsize=(14.8, 9.4), layout="constrained")
    grid = fig.add_gridspec(3, 3, height_ratios=(0.10, 1, 1))
    legend_axis = fig.add_subplot(grid[0, :])
    legend_axis.axis("off")
    axes = [
        [fig.add_subplot(grid[1, column]) for column in range(3)],
        [fig.add_subplot(grid[2, column]) for column in range(3)],
    ]
    _scatter_panel(axes[0][0], rows, "runtime_seconds", "Wall time (s)", "A   Wall time", log_y=True)
    _scatter_panel(
        axes[0][1], rows, "throughput_Mbp_per_second", "Throughput (Mbp s⁻¹)",
        "B   Observed throughput", log_y=True,
    )
    _scatter_panel(axes[0][2], rows, "peak_rss_mib", "Peak RSS (MiB)", "C   Peak memory")

    ratio_rows = [row for row in rows if row["tool"] != "trf"]
    _scatter_panel(
        axes[1][0], ratio_rows, "runtime_ratio_to_trf", "Wall-time ratio to TRF",
        "D   Paired wall-time ratio", log_y=True,
    )
    axes[1][0].axhline(1, color="#555555", ls="--", lw=0.9)
    _scatter_panel(
        axes[1][1], ratio_rows, "rss_ratio_to_trf", "Peak-RSS ratio to TRF",
        "E   Paired memory ratio", log_y=True,
    )
    axes[1][1].axhline(1, color="#555555", ls="--", lw=0.9)
    _scatter_panel(
        axes[1][2], rows, "observed_union_base_fraction", "Called union / input bases",
        "F   Descriptive called-base fraction",
    )

    tool_handles = [
        Line2D([], [], marker="o", linestyle="", color=COLORS[tool], label=DISPLAY[tool], markersize=6)
        for tool in TOOLS
    ]
    materials = sorted({run.material for run in runs})
    material_handles = [
        Line2D([], [], marker=MARKERS[index % len(MARKERS)], linestyle="", color="#555555",
               label=material, markersize=6)
        for index, material in enumerate(materials)
    ]
    tool_legend = legend_axis.legend(
        handles=tool_handles, loc="center left", bbox_to_anchor=(0.03, 0.32),
        ncols=3, frameon=False, title="Tool",
    )
    legend_axis.add_artist(tool_legend)
    legend_axis.legend(
        handles=material_handles, loc="center right", bbox_to_anchor=(0.97, 0.32),
        ncols=min(5, len(material_handles)), frameon=False, title="Material marker",
    )
    fig.suptitle(
        f"One-thread real-read diagnostics: {len(runs)} nested inputs across {len(materials)} materials",
        fontsize=13,
        fontweight="bold",
    )
    fig.supxlabel(
        "Single executions overlapped acquisition work and are not final resource rankings. "
        "Panel F is output extent, not accuracy; nested inputs are not independent biological replicates.",
        fontsize=8,
    )

    outdir.mkdir(parents=True)
    fields = list(rows[0])
    source_path = outdir / "panel_source.tsv"
    write_table(source_path, rows, fields)
    outputs = []
    for extension in ("pdf", "svg", "png"):
        path = outdir / f"real_comparator_diagnostics.{extension}"
        fig.savefig(path, dpi=220)
        outputs.append(path)
    plt.close(fig)
    provenance = {
        "complete": True,
        "run_count": len(runs),
        "material_count": len(materials),
        "tools": list(TOOLS),
        "materials": {run.run_id: run.material for run in runs},
        "design": "one thread, one execution, nested whole-library random samples; concurrent development diagnostics",
        "accuracy": "not_assessed_without_curated_independent_truth",
        "run_archives": {
            run.run_id: {
                "material": run.material,
                "archive_manifest_sha256": digest_file(run.path / "archive_manifest.json"),
            }
            for run in runs
        },
        "source_sha256": digest_file(source_path),
        "script_sha256": digest_file(Path(__file__)),
        "outputs": {path.name: digest_file(path) for path in outputs},
        "interpretation": "resource and output-extent diagnostics; not an accuracy or final performance ranking",
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, type=parse_run,
                        help="Completed compact archive as MATERIAL=PATH; repeat for each run")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.run, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
