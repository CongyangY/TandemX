"""Render source-backed cross-library QC as an editable six-panel figure."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.challenge.schema import digest_file, read_table, write_table


@dataclass(frozen=True)
class LibraryQC:
    run_accession: str
    species: str
    material: str
    label: str
    read_count: int
    total_bases: int
    median_length: int
    read_n50: int
    gc_fraction: float
    length_histogram: Counter[int]
    gc_histogram: Counter[int]
    quality_histogram: Counter[int]
    inputs: tuple[Path, ...]


def _integer_histogram(path: Path, value_field: str, expected_reads: int) -> Counter[int]:
    result: Counter[int] = Counter()
    for row in read_table(path):
        value = int(row[value_field])
        count = int(row["read_count"])
        if value < 0 or count <= 0:
            raise ValueError(f"Invalid histogram value/count in {path}")
        result[value] += count
    if sum(result.values()) != expected_reads:
        raise ValueError(f"Histogram total differs from cohort receipt: {path}")
    return result


def _joint_histograms(path: Path, expected_reads: int) -> tuple[Counter[int], Counter[int]]:
    gc: Counter[int] = Counter()
    quality: Counter[int] = Counter()
    total = 0
    for row in read_table(path):
        count = int(row["read_count"])
        gc_bin = int(row["gc_bin_percent"])
        quality_bin = int(row["mean_quality_bin_phred"])
        if count <= 0 or gc_bin < 0 or quality_bin < 0:
            raise ValueError(f"Invalid joint-distribution value/count in {path}")
        gc[gc_bin] += count
        quality[quality_bin] += count
        total += count
    if total != expected_reads:
        raise ValueError(f"Joint-distribution total differs from cohort receipt: {path}")
    return gc, quality


def _abbreviate_species(species: str) -> str:
    words = species.split()
    return f"{words[0][0]}. {' '.join(words[1:])}" if len(words) > 1 else species


def load_cohort(cohort_table: Path, repo_root: Path) -> list[LibraryQC]:
    rows = read_table(cohort_table)
    if not rows:
        raise ValueError("Cohort table is empty")
    libraries: list[LibraryQC] = []
    seen_runs: set[str] = set()
    for row in rows:
        required = {
            "reported_species", "reported_material", "run_accession", "read_count",
            "total_bases", "median_read_length", "read_n50", "gc_fraction",
            "raw_sha256", "qc_receipt", "qc_receipt_sha256",
        }
        if required - row.keys():
            raise ValueError("Cohort table lacks required QC columns")
        run = row["run_accession"]
        if not run or run in seen_runs:
            raise ValueError(f"Run accessions must be nonempty and unique: {run}")
        seen_runs.add(run)
        receipt_path = repo_root / row["qc_receipt"]
        if not receipt_path.is_file() or digest_file(receipt_path) != row["qc_receipt_sha256"]:
            raise ValueError(f"QC receipt hash differs for {run}")
        receipt = json.loads(receipt_path.read_text())
        expected = {
            "read_count": int(row["read_count"]),
            "total_bases": int(row["total_bases"]),
            "median_length": int(row["median_read_length"]),
            "read_n50": int(row["read_n50"]),
            "input_sha256": row["raw_sha256"],
        }
        if receipt.get("complete") is not True or any(receipt.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Cohort fields differ from complete QC receipt for {run}")
        if abs(float(receipt["gc_fraction"]) - float(row["gc_fraction"])) > 1e-12:
            raise ValueError(f"GC fraction differs from QC receipt for {run}")
        length_path = receipt_path.parent / "length_histogram.tsv"
        joint_path = receipt_path.parent / "joint_distribution.tsv"
        length = _integer_histogram(length_path, "length_bp", expected["read_count"])
        gc, quality = _joint_histograms(joint_path, expected["read_count"])
        species, material = row["reported_species"], row["reported_material"]
        libraries.append(
            LibraryQC(
                run_accession=run,
                species=species,
                material=material,
                label=f"{material} ({_abbreviate_species(species)})",
                read_count=expected["read_count"],
                total_bases=expected["total_bases"],
                median_length=expected["median_length"],
                read_n50=expected["read_n50"],
                gc_fraction=float(receipt["gc_fraction"]),
                length_histogram=length,
                gc_histogram=gc,
                quality_histogram=quality,
                inputs=(receipt_path, length_path, joint_path),
            )
        )
    return sorted(libraries, key=lambda item: (item.species, item.material, item.run_accession))


def aggregate_histogram(
    histogram: Counter[int], *, width: int, lower: int = 0, upper: int
) -> tuple[list[int], list[int]]:
    if width <= 0 or lower < 0 or upper <= lower or (upper - lower) % width:
        raise ValueError("Require aligned positive histogram bounds")
    bins = list(range(lower, upper + 1, width))
    counts = Counter()
    for value, count in histogram.items():
        clipped = lower if value < lower else upper if value >= upper else value // width * width
        counts[clipped] += count
    return bins, [counts[value] for value in bins]


def _heatmap(
    axis,
    libraries: list[LibraryQC],
    attribute: str,
    *,
    width: int,
    lower: int,
    upper: int,
    title: str,
    unit_scale: float,
    show_y: bool,
):
    bins: list[int] = []
    counts_by_library: list[list[int]] = []
    for library in libraries:
        current_bins, counts = aggregate_histogram(
            getattr(library, attribute), width=width, lower=lower, upper=upper
        )
        if bins and current_bins != bins:
            raise RuntimeError("Histogram bins changed between libraries")
        bins = current_bins
        counts_by_library.append(counts)
    fractions = np.array(
        [[count / library.read_count for count in counts] for library, counts in zip(libraries, counts_by_library)]
    )
    mesh = axis.pcolormesh(
        np.arange(len(bins) + 1), np.arange(len(libraries) + 1), fractions,
        cmap="viridis", shading="flat", rasterized=False,
    )
    axis.set_ylim(len(libraries), 0)
    axis.set_yticks(np.arange(len(libraries)) + 0.5)
    axis.set_yticklabels([library.label for library in libraries] if show_y else [])
    tick_indices = sorted(set([0, len(bins) - 1, *range(0, len(bins), max(1, len(bins) // 5))]))
    tick_labels = []
    for index in tick_indices:
        value = bins[index] * unit_scale
        tick_labels.append(("≥" if index == len(bins) - 1 else "") + f"{value:g}")
    axis.set_xticks(np.array(tick_indices) + 0.5, tick_labels)
    axis.set_title(title, loc="left", fontweight="bold")
    colorbar = axis.figure.colorbar(mesh, ax=axis, fraction=0.046, pad=0.03)
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    colorbar.set_label("Read fraction")
    return bins, counts_by_library, fractions


def plot(cohort_table: Path, outdir: Path, repo_root: Path) -> None:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    libraries = load_cohort(cohort_table, repo_root)
    if len(libraries) < 2:
        raise ValueError("Cross-library figure requires at least two libraries")

    plt.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.family": "DejaVu Sans",
            "font.size": 8.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(15.8, 9.8), layout="constrained")
    labels = [library.label for library in libraries]
    y = np.arange(len(libraries))
    colors = plt.get_cmap("tab10")(np.arange(len(libraries)) % 10)
    source: list[dict] = []

    total_gb = np.array([library.total_bases / 1e9 for library in libraries])
    axes[0, 0].barh(y, total_gb, color=colors, alpha=0.85)
    axes[0, 0].set_yticks(y, labels)
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlabel("Sequence bases (Gb)")
    axes[0, 0].set_title("A   Validated data volume", loc="left", fontweight="bold")
    axes[0, 0].grid(axis="x", alpha=0.15)
    for index, (library, value) in enumerate(zip(libraries, total_gb)):
        axes[0, 0].text(value, index, f" {value:.1f}", va="center", fontsize=7.5)
        source.append(_source_row("A", library, "total_bases_Gb", "NA", "NA", value, library.read_count))

    medians = np.array([library.median_length / 1000 for library in libraries])
    n50s = np.array([library.read_n50 / 1000 for library in libraries])
    for index, library in enumerate(libraries):
        axes[0, 1].plot([medians[index], n50s[index]], [index, index], color="#999999", lw=1.2)
    axes[0, 1].scatter(medians, y, color="#4477AA", marker="o", label="Median", zorder=3)
    axes[0, 1].scatter(n50s, y, color="#CC6677", marker="D", label="Read N50", zorder=3)
    axes[0, 1].set_yticks(y, [])
    axes[0, 1].set_ylim(len(libraries) - 0.5, -0.5)
    axes[0, 1].set_xlabel("Read length (kb)")
    axes[0, 1].set_title("B   Read-length summaries", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False, ncols=2, loc="lower right")
    axes[0, 1].grid(axis="x", alpha=0.15)
    for library, median, n50 in zip(libraries, medians, n50s):
        source.extend(
            [
                _source_row("B", library, "median_read_length_kb", "NA", "NA", median, library.read_count),
                _source_row("B", library, "read_n50_kb", "NA", "NA", n50, library.read_count),
            ]
        )

    gc_percent = np.array([library.gc_fraction * 100 for library in libraries])
    axes[0, 2].hlines(y, gc_percent.min() - 1, gc_percent, color="#BBBBBB", lw=1)
    axes[0, 2].scatter(gc_percent, y, color=colors, edgecolor="black", linewidth=0.3, s=42)
    axes[0, 2].set_yticks(y, [])
    axes[0, 2].set_ylim(len(libraries) - 0.5, -0.5)
    axes[0, 2].set_xlabel("GC bases (%)")
    axes[0, 2].set_title("C   Whole-file GC fraction", loc="left", fontweight="bold")
    axes[0, 2].grid(axis="x", alpha=0.15)
    for library, value in zip(libraries, gc_percent):
        source.append(_source_row("C", library, "gc_percent", "NA", "NA", value, library.read_count))

    heatmaps = [
        ("D", axes[1, 0], "length_histogram", 2000, 0, 80000, "D   Read-length distribution", 0.001, True),
        ("E", axes[1, 1], "gc_histogram", 2, 20, 70, "E   Read-GC distribution", 1.0, False),
        ("F", axes[1, 2], "quality_histogram", 5, 0, 60, "F   Reported read-quality distribution", 1.0, False),
    ]
    for panel, axis, attribute, width, lower, upper, title, scale, show_y in heatmaps:
        bins, counts_by_library, fractions = _heatmap(
            axis, libraries, attribute, width=width, lower=lower, upper=upper,
            title=title, unit_scale=scale, show_y=show_y,
        )
        axis.set_xlabel(
            "Read length bin lower bound (kb)" if panel == "D"
            else "GC bin lower bound (%)" if panel == "E"
            else "Mean reported Phred bin lower bound"
        )
        for row_index, library in enumerate(libraries):
            for column, bin_value in enumerate(bins):
                label = f"at_least_{bin_value}" if column == len(bins) - 1 else str(bin_value)
                source.append(
                    _source_row(
                        panel, library, attribute, bin_value, label,
                        float(fractions[row_index, column]), counts_by_library[row_index][column],
                    )
                )

    fig.suptitle(
        f"Cross-cohort QC of {len(libraries)} complete plant HiFi libraries",
        fontsize=13,
        fontweight="bold",
    )
    fig.supxlabel(
        "All distributions use complete checksum-validated files. Reported quality is not empirical accuracy; "
        "libraries and nested samples are not biological replicates.",
        fontsize=8,
    )
    outdir.mkdir(parents=True)
    source_path = outdir / "panel_source.tsv"
    fields = [
        "panel", "run_accession", "reported_species", "reported_material", "metric",
        "bin_lower_bound", "bin_label", "value", "read_count", "library_read_count",
    ]
    write_table(source_path, source, fields)
    outputs = []
    for extension in ("pdf", "svg", "png"):
        path = outdir / f"multispecies_input_qc.{extension}"
        fig.savefig(path, dpi=220)
        outputs.append(path)
    plt.close(fig)
    inputs = [cohort_table, *(path for library in libraries for path in library.inputs)]
    provenance = {
        "complete": True,
        "library_count": len(libraries),
        "reported_species_count": len({library.species for library in libraries}),
        "read_count": sum(library.read_count for library in libraries),
        "total_bases": sum(library.total_bases for library in libraries),
        "library_order": [library.run_accession for library in libraries],
        "inputs": {str(path.resolve()): digest_file(path) for path in inputs},
        "source_sha256": digest_file(source_path),
        "script_sha256": digest_file(Path(__file__)),
        "outputs": {path.name: digest_file(path) for path in outputs},
        "interpretation": (
            "complete-file input distributions; not biological accuracy, repeat truth, "
            "independent replication or empirical sequencing error"
        ),
    }
    (outdir / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def _source_row(
    panel: str,
    library: LibraryQC,
    metric: str,
    bin_lower_bound: int | str,
    bin_label: str,
    value: float,
    read_count: int,
) -> dict:
    return {
        "panel": panel,
        "run_accession": library.run_accession,
        "reported_species": library.species,
        "reported_material": library.material,
        "metric": metric,
        "bin_lower_bound": bin_lower_bound,
        "bin_label": bin_label,
        "value": value,
        "read_count": read_count,
        "library_read_count": library.read_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-table", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    plot(args.cohort_table, args.outdir, args.repo_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
