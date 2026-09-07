"""Publication-oriented static overview for TandemX cohort outputs."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_source(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render_cohort_overview(outdir: Path, top_families: int = 30) -> list[Path]:
    if top_families < 1:
        raise ValueError("--top-families must be positive")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matplotlib.rcParams["svg.fonttype"] = "none"
    inputs = {
        "families": outdir / "pan_families.tsv",
        "abundance": outdir / "sample_family_abundance.tsv",
        "representation": outdir / "sample_family_representation.tsv",
    }
    if any(not path.is_file() for path in inputs.values()):
        raise ValueError("Cohort overview requires completed cohort tables")
    family_rows = _read(inputs["families"])
    abundance_rows = _read(inputs["abundance"])
    representation_rows = _read(inputs["representation"])
    if not family_rows or not abundance_rows or not representation_rows:
        raise ValueError("Cohort overview inputs must contain records")
    samples = list(dict.fromkeys(row["sample_id"] for row in abundance_rows))
    families = {row["pan_family_id"]: row for row in family_rows}
    abundance = {(row["sample_id"], row["pan_family_id"]): row for row in abundance_rows}
    representation = {
        (row["sample_id"], row["pan_family_id"]): row for row in representation_rows
    }
    expected = {(sample, family) for sample in samples for family in families}
    if set(abundance) != expected or set(representation) != expected:
        raise ValueError("Cohort long tables do not form complete sample-by-family grids")

    def peak_bp(family: str) -> float:
        values = [abundance[(sample, family)]["estimated_bp"] for sample in samples]
        finite = [float(value) for value in values if value != "NA"]
        return max(finite, default=-1.0)

    selected = sorted(families, key=lambda family: (-peak_bp(family), family))[
        :top_families
    ]
    shape = (len(selected), len(samples))
    abundance_matrix = np.full(shape, np.nan)
    representation_matrix = np.full(shape, np.nan)
    interval_width_matrix = np.full(shape, np.nan)
    source_rows: list[dict[str, object]] = []
    for i, family in enumerate(selected):
        for j, sample in enumerate(samples):
            abundance_row = abundance[(sample, family)]
            representation_row = representation[(sample, family)]
            estimated = abundance_row["estimated_bp"]
            low = abundance_row["estimated_bp_interval_low"]
            high = abundance_row["estimated_bp_interval_high"]
            ratio = representation_row["assembly_read_ratio"]
            if estimated != "NA":
                estimated_value = float(estimated)
                abundance_matrix[i, j] = np.log10(estimated_value + 1)
                if low != "NA" and high != "NA" and estimated_value > 0:
                    interval_width_matrix[i, j] = (
                        float(high) - float(low)
                    ) / estimated_value
            if ratio != "NA":
                representation_matrix[i, j] = float(ratio)
            source_rows.append(
                {
                    "sample_id": sample,
                    "pan_family_id": family,
                    "selection_rank": i + 1,
                    "estimated_bp": estimated,
                    "estimated_bp_interval_low": low,
                    "estimated_bp_interval_high": high,
                    "relative_interval_width": (
                        "NA"
                        if np.isnan(interval_width_matrix[i, j])
                        else f"{interval_width_matrix[i, j]:.6f}"
                    ),
                    "assembly_read_ratio": ratio,
                    "abundance_status": abundance_row["status"],
                    "representation_status": representation_row["status"],
                    "sample_count": families[family]["sample_count"],
                    "warning": "top_families_ranked_by_max_sample_estimated_bp",
                }
            )

    height = max(9.5, 0.32 * len(selected) + 5)
    width = max(12.5, 1.0 * len(samples) + 8)
    fig, axes = plt.subplots(2, 2, figsize=(width, height), constrained_layout=True)
    abundance_cmap = plt.colormaps["viridis"].copy()
    abundance_cmap.set_bad("#eeeeee")
    ratio_cmap = plt.colormaps["coolwarm"].copy()
    ratio_cmap.set_bad("#eeeeee")
    width_cmap = plt.colormaps["magma"].copy()
    width_cmap.set_bad("#eeeeee")

    masked_abundance = np.ma.masked_invalid(abundance_matrix)
    image = axes[0, 0].pcolormesh(
        masked_abundance, cmap=abundance_cmap, shading="flat", rasterized=False
    )
    _heatmap_axes(axes[0, 0], samples, selected)
    axes[0, 0].set_title("A  Read-estimated abundance", loc="left", fontweight="bold")
    colorbar = fig.colorbar(
        image, ax=axes[0, 0], label="log10(estimated bp + 1)", shrink=0.8
    )
    colorbar.solids.set_rasterized(False)

    masked_ratio = np.ma.masked_invalid(representation_matrix)
    image = axes[0, 1].pcolormesh(
        masked_ratio, cmap=ratio_cmap, vmin=0, vmax=2, shading="flat", rasterized=False
    )
    _heatmap_axes(axes[0, 1], samples, selected)
    axes[0, 1].set_title("B  Assembly/read representation", loc="left", fontweight="bold")
    colorbar = fig.colorbar(
        image,
        ax=axes[0, 1],
        label="assembly/read ratio (clipped at 2)",
        shrink=0.8,
    )
    colorbar.solids.set_rasterized(False)

    prevalence = [int(families[family]["sample_count"]) for family in selected]
    y = np.arange(len(selected))
    axes[1, 0].barh(y, prevalence, color="#59a14f")
    axes[1, 0].set_yticks(y, selected)
    axes[1, 0].invert_yaxis()
    axes[1, 0].set_xlim(0, max(1, len(samples)))
    axes[1, 0].set_xlabel("samples with a catalogue member")
    axes[1, 0].set_title("C  Pan-family prevalence", loc="left", fontweight="bold")

    masked_width = np.ma.masked_invalid(interval_width_matrix)
    image = axes[1, 1].pcolormesh(
        masked_width, cmap=width_cmap, shading="flat", rasterized=False
    )
    _heatmap_axes(axes[1, 1], samples, selected)
    axes[1, 1].set_title("D  Propagated interval width", loc="left", fontweight="bold")
    colorbar = fig.colorbar(
        image, ax=axes[1, 1], label="(high - low) / estimate", shrink=0.8
    )
    colorbar.solids.set_rasterized(False)

    fig.suptitle(
        f"TandemX cohort overview: top {len(selected)} operational pan families\n"
        "NA is unavailable evidence, not zero abundance",
        fontweight="bold",
    )
    svg = outdir / "cohort_overview.svg"
    pdf = outdir / "cohort_overview.pdf"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    source = outdir / "cohort_plot_source.tsv"
    _write_source(source, source_rows)
    svg_text = svg.read_text(encoding="utf-8")
    receipt_path = outdir / "cohort_figure_receipt.json"
    receipt = {
        "schema_version": 1,
        "complete": True,
        "panel_count": 4,
        "top_families_requested": top_families,
        "top_families_rendered": len(selected),
        "selected_pan_family_ids": selected,
        "input_sha256": {name: _digest(path) for name, path in inputs.items()},
        "output_sha256": {path.name: _digest(path) for path in (svg, pdf, source)},
        "svg_text_nodes": svg_text.count("<text"),
        "svg_image_nodes": svg_text.count("<image"),
        "warning": "operational_pan_families;propagated_endpoints_not_joint_confidence_intervals",
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return [svg, pdf, source, receipt_path]


def _heatmap_axes(ax, samples: list[str], families: list[str]) -> None:
    ax.set_xticks(
        [index + 0.5 for index in range(len(samples))],
        samples,
        rotation=35,
        ha="right",
    )
    ax.set_yticks([index + 0.5 for index in range(len(families))], families)
    ax.set_xlim(0, len(samples))
    ax.set_ylim(len(families), 0)
    ax.set_xlabel("sample")
