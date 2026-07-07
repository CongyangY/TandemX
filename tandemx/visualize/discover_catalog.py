"""Static report plots for ``tandemx discover`` outputs."""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

from tandemx.io.sequences import read_sequence_records


@dataclass(frozen=True)
class DiscoverCatalogConfig:
    families_tsv: Path
    family_similarity_tsv: Path
    monomers_fa: Path
    outdir: Path
    top_n: int = 20


@dataclass(frozen=True)
class FamilyRecord:
    family_id: str
    monomer_id: str
    monomer_length_bp: int
    sequence_length_bp: int
    consensus_md5: str
    gc_fraction: float
    support_read_count: int
    support_span_bp: int
    mean_identity: float
    low_complexity_flag: bool
    confidence: str
    warning: str
    warning_labels: tuple[str, ...]

    @property
    def has_relationship_warning(self) -> bool:
        return any(label.startswith("possible_higher_order_or_partial:") for label in self.warning_labels)

    @property
    def status_label(self) -> str:
        if self.low_complexity_flag and self.has_relationship_warning:
            return "low_complexity + related"
        if self.has_relationship_warning:
            return "related"
        if self.low_complexity_flag:
            return "low_complexity"
        return "clean"

    @property
    def short_label(self) -> str:
        return abbreviate_family_id(self.family_id)


@dataclass(frozen=True)
class SimilarityRecord:
    family_a: str
    family_b: str
    length_a_bp: int
    length_b_bp: int
    kmer_jaccard: float
    shared_kmer_fraction: float
    local_identity: float
    local_overlap_bp: int
    local_overlap_fraction_shorter: float
    length_ratio: float
    orientation: str
    relationship: str
    redundant_candidate: bool
    notes: str

    @property
    def pair_label(self) -> str:
        return f"{self.family_a}-{self.family_b}"

    @property
    def short_pair_label(self) -> str:
        return f"{abbreviate_family_id(self.family_a)}-{abbreviate_family_id(self.family_b)}"


@dataclass(frozen=True)
class DiscoverCatalogSummary:
    family_count: int
    low_complexity_count: int
    relationship_warning_count: int
    warning_count: int
    non_distinct_pair_count: int
    median_length_bp: float
    total_support_reads: int
    total_support_bp: int
    confidence_counts: tuple[tuple[str, int], ...]
    status_counts: tuple[tuple[str, int], ...]
    top_family_id: str
    top_family_span_bp: int


def render_discover_catalog_report(config: DiscoverCatalogConfig) -> list[Path]:
    cache_root = config.outdir / ".plot_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config.outdir.mkdir(parents=True, exist_ok=True)
    families = load_families(config.families_tsv, config.monomers_fa)
    similarities = load_similarities(config.family_similarity_tsv)
    summary = summarize_catalog(families, similarities)

    outputs: list[Path] = []
    outputs.extend(write_tables(config.outdir, families, similarities, config.top_n))
    outputs.extend(plot_family_abundance(plt, config.outdir, families, config.top_n))
    outputs.extend(plot_length_distribution(plt, config.outdir, families))
    outputs.extend(plot_abundance_vs_length(plt, config.outdir, families, config.top_n))
    outputs.extend(plot_quality_overview(plt, config.outdir, families, config.top_n))
    outputs.extend(plot_similarity_space(plt, config.outdir, similarities))
    outputs.extend(plot_flagged_pairs_heatmap(plt, config.outdir, similarities))
    outputs.append(write_markdown_report(config.outdir, summary, families, similarities, config.top_n))
    return outputs


def load_families(families_tsv: Path, monomers_fa: Path) -> list[FamilyRecord]:
    lengths = monomer_lengths_by_family(monomers_fa)
    rows: list[FamilyRecord] = []
    with families_tsv.open("rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels = tuple(item for item in row["warning"].split(";") if item)
            rows.append(
                FamilyRecord(
                    family_id=row["family_id"],
                    monomer_id=row["monomer_id"],
                    monomer_length_bp=int(row["monomer_length_bp"]),
                    sequence_length_bp=lengths.get(row["family_id"], int(row["monomer_length_bp"])),
                    consensus_md5=row["consensus_md5"],
                    gc_fraction=float(row["gc_fraction"]),
                    support_read_count=int(row["support_read_count"]),
                    support_span_bp=int(row["support_span_bp"]),
                    mean_identity=float(row["mean_identity"]),
                    low_complexity_flag=row["low_complexity_flag"].lower() == "true",
                    confidence=row["confidence"],
                    warning=row["warning"],
                    warning_labels=labels,
                )
            )
    return rows


def load_similarities(path: Path) -> list[SimilarityRecord]:
    rows: list[SimilarityRecord] = []
    with path.open("rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                SimilarityRecord(
                    family_a=row["family_a"],
                    family_b=row["family_b"],
                    length_a_bp=int(row["length_a_bp"]),
                    length_b_bp=int(row["length_b_bp"]),
                    kmer_jaccard=float(row["kmer_jaccard"]),
                    shared_kmer_fraction=float(row["shared_kmer_fraction"]),
                    local_identity=float(row["local_identity"]),
                    local_overlap_bp=int(row["local_overlap_bp"]),
                    local_overlap_fraction_shorter=float(row["local_overlap_fraction_shorter"]),
                    length_ratio=float(row["length_ratio"]),
                    orientation=row["orientation"],
                    relationship=row["relationship"],
                    redundant_candidate=row["redundant_candidate"].lower() == "true",
                    notes=row["notes"],
                )
            )
    return rows


def monomer_lengths_by_family(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for record in read_sequence_records(path):
        family_id = parse_header_field(record.description, "family_id")
        lengths[family_id] = len(record.sequence)
    return lengths


def parse_header_field(description: str, field: str) -> str:
    for part in description.split(";"):
        key, _, value = part.partition("=")
        if key == field:
            return value
    raise ValueError(f"Missing {field} in FASTA header: {description}")


def summarize_catalog(
    families: Sequence[FamilyRecord],
    similarities: Sequence[SimilarityRecord],
) -> DiscoverCatalogSummary:
    confidence_counts = count_labels(family.confidence for family in families)
    status_counts = count_labels(family.status_label for family in families)
    top_family = max(families, key=lambda item: item.support_span_bp)
    return DiscoverCatalogSummary(
        family_count=len(families),
        low_complexity_count=sum(family.low_complexity_flag for family in families),
        relationship_warning_count=sum(family.has_relationship_warning for family in families),
        warning_count=sum(bool(family.warning) for family in families),
        non_distinct_pair_count=sum(item.relationship != "distinct" for item in similarities),
        median_length_bp=median(family.monomer_length_bp for family in families),
        total_support_reads=sum(family.support_read_count for family in families),
        total_support_bp=sum(family.support_span_bp for family in families),
        confidence_counts=tuple(sorted(confidence_counts.items())),
        status_counts=tuple(sorted(status_counts.items())),
        top_family_id=top_family.family_id,
        top_family_span_bp=top_family.support_span_bp,
    )


def count_labels(labels: Iterable[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return counts


def write_tables(outdir: Path, families: Sequence[FamilyRecord], similarities: Sequence[SimilarityRecord], top_n: int) -> list[Path]:
    top_path = outdir / "top_families.tsv"
    flagged_path = outdir / "flagged_pairs.tsv"
    label_map_path = outdir / "family_label_map.tsv"
    top_families = sorted(families, key=lambda item: item.support_span_bp, reverse=True)[:top_n]
    flagged_pairs = [item for item in similarities if item.relationship != "distinct"]

    with top_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "family_id",
                "monomer_length_bp",
                "sequence_length_bp",
                "support_read_count",
                "support_span_bp",
                "mean_identity",
                "gc_fraction",
                "low_complexity_flag",
                "confidence",
                "status_label",
                "warning",
            ]
        )
        for family in top_families:
            writer.writerow(
                [
                    family.family_id,
                    family.monomer_length_bp,
                    family.sequence_length_bp,
                    family.support_read_count,
                    family.support_span_bp,
                    f"{family.mean_identity:.4f}",
                    f"{family.gc_fraction:.4f}",
                    str(family.low_complexity_flag).lower(),
                    family.confidence,
                    family.status_label,
                    family.warning,
                ]
            )

    with flagged_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "family_a",
                "family_b",
                "relationship",
                "kmer_jaccard",
                "shared_kmer_fraction",
                "local_identity",
                "local_overlap_bp",
                "local_overlap_fraction_shorter",
                "length_ratio",
                "orientation",
                "notes",
            ]
        )
        for pair in sorted(flagged_pairs, key=lambda item: (item.relationship, -item.shared_kmer_fraction, -item.local_identity)):
            writer.writerow(
                [
                    pair.family_a,
                    pair.family_b,
                    pair.relationship,
                    f"{pair.kmer_jaccard:.4f}",
                    f"{pair.shared_kmer_fraction:.4f}",
                    f"{pair.local_identity:.4f}",
                    pair.local_overlap_bp,
                    f"{pair.local_overlap_fraction_shorter:.4f}",
                    f"{pair.length_ratio:.4f}",
                    pair.orientation,
                    pair.notes,
                ]
            )
    with label_map_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["family_id", "short_label"])
        for family in sorted(families, key=lambda item: item.family_id):
            writer.writerow([family.family_id, family.short_label])
    return [top_path, flagged_path, label_map_path]


def plot_family_abundance(plt, outdir: Path, families: Sequence[FamilyRecord], top_n: int) -> list[Path]:
    top_families = sorted(families, key=lambda item: item.support_span_bp, reverse=True)[:top_n]
    labels = [family.short_label for family in top_families]
    spans_mb = [family.support_span_bp / 1_000_000 for family in top_families]
    read_counts_k = [family.support_read_count / 1_000 for family in top_families]
    colors = [status_color(family.status_label) for family in top_families]

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(13.5, max(7.8, top_n * 0.48)), sharey=True)
    y = list(range(len(top_families)))
    ax1.barh(y, spans_mb, color=colors)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=15)
    ax1.invert_yaxis()
    ax1.set_xlabel("Supporting span (Mb)", fontsize=17)
    ax1.set_title(f"Top {top_n} families by support span", fontsize=22)
    for idx, family in enumerate(top_families):
        ax1.text(spans_mb[idx], idx, f"  {family.monomer_length_bp} bp", va="center", fontsize=12)

    ax2.barh(y, read_counts_k, color=colors)
    ax2.set_xlabel("Supporting reads (thousands)", fontsize=17)
    ax2.set_title("Read support", fontsize=22)
    ax1.tick_params(axis="x", labelsize=14)
    ax2.tick_params(axis="x", labelsize=14)
    add_status_legend(ax2)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_family_abundance")


def plot_length_distribution(plt, outdir: Path, families: Sequence[FamilyRecord]) -> list[Path]:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(13.5, 5.8))
    lengths = [family.monomer_length_bp for family in families]
    weighted_lengths = [family.monomer_length_bp for family in families for _ in range(min(25, max(1, int(math.log10(max(10, family.support_span_bp))))))]
    bins = [1, 5, 10, 20, 50, 100, 200, 400, 800, 1200, 1600, 2000, 2400]
    ax1.hist(lengths, bins=bins, color="#4C78A8", edgecolor="white")
    ax1.set_xscale("log")
    ax1.set_xlabel("Monomer length (bp)", fontsize=17)
    ax1.set_ylabel("Family count", fontsize=17)
    ax1.set_title("Monomer length distribution", fontsize=22)

    ax2.hist(weighted_lengths, bins=bins, color="#F58518", edgecolor="white")
    ax2.set_xscale("log")
    ax2.set_xlabel("Monomer length (bp)", fontsize=17)
    ax2.set_ylabel("Abundance-weighted count", fontsize=17)
    ax2.set_title("Length distribution weighted by support", fontsize=22)
    ax1.tick_params(axis="both", labelsize=14)
    ax2.tick_params(axis="both", labelsize=14)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_length_distribution")


def plot_abundance_vs_length(plt, outdir: Path, families: Sequence[FamilyRecord], top_n: int) -> list[Path]:
    fig, ax = plt.subplots(figsize=(9.8, 7.4))
    families_sorted = sorted(families, key=lambda item: item.support_span_bp, reverse=True)
    for family in families_sorted:
        ax.scatter(
            family.monomer_length_bp,
            family.support_span_bp,
            s=30 + 8 * math.log10(max(10, family.support_read_count)),
            c=status_color(family.status_label),
            edgecolors="#222222" if family.has_relationship_warning else "none",
            linewidths=0.6 if family.has_relationship_warning else 0.0,
            alpha=0.85,
        )
    for family in families_sorted[:top_n]:
        ax.annotate(
            family.short_label,
            (family.monomer_length_bp, family.support_span_bp),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=11,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Monomer length (bp)", fontsize=17)
    ax.set_ylabel("Supporting span (bp)", fontsize=17)
    ax.set_title("Family abundance versus monomer length", fontsize=22)
    ax.tick_params(axis="both", labelsize=14)
    add_status_legend(ax)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_abundance_vs_length")


def plot_quality_overview(plt, outdir: Path, families: Sequence[FamilyRecord], top_n: int) -> list[Path]:
    fig, ax = plt.subplots(figsize=(9.8, 7.4))
    ranked = sorted(families, key=lambda item: item.support_span_bp, reverse=True)
    for family in ranked:
        ax.scatter(
            family.gc_fraction,
            family.mean_identity,
            s=35 + 10 * math.log10(max(10, family.support_span_bp)),
            c=status_color(family.status_label),
            alpha=0.85,
        )
    for family in ranked[:top_n]:
        ax.annotate(
            family.short_label,
            (family.gc_fraction, family.mean_identity),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=11,
        )
    ax.set_xlabel("GC fraction", fontsize=17)
    ax.set_ylabel("Mean identity", fontsize=17)
    ax.set_title("Family sequence composition and support quality", fontsize=22)
    ax.tick_params(axis="both", labelsize=14)
    add_status_legend(ax)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_quality_overview")


def plot_similarity_space(plt, outdir: Path, similarities: Sequence[SimilarityRecord]) -> list[Path]:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(14.5, 7.0))
    distinct = [pair for pair in similarities if pair.relationship == "distinct"]
    flagged = [pair for pair in similarities if pair.relationship != "distinct"]
    annotated = select_similarity_annotations(flagged, max_labels=12, per_family_cap=3)

    if distinct:
        ax1.scatter(
            [pair.shared_kmer_fraction for pair in distinct],
            [pair.local_identity for pair in distinct],
            c="#9D9D9D",
            alpha=0.22,
            s=10,
            rasterized=True,
        )
        ax2.scatter(
            [pair.length_ratio for pair in distinct],
            [pair.local_overlap_fraction_shorter for pair in distinct],
            c="#9D9D9D",
            alpha=0.22,
            s=10,
            rasterized=True,
        )
    if flagged:
        sizes = [35 + pair.local_overlap_bp / 12 for pair in flagged]
        ax1.scatter(
            [pair.shared_kmer_fraction for pair in flagged],
            [pair.local_identity for pair in flagged],
            c="#E45756",
            alpha=0.9,
            s=sizes,
        )
        ax2.scatter(
            [pair.length_ratio for pair in flagged],
            [pair.local_overlap_fraction_shorter for pair in flagged],
            c="#E45756",
            alpha=0.9,
            s=sizes,
        )
        for index, pair in enumerate(annotated):
            left_offset = (6, 6 + (index % 3) * 10)
            right_offset = (6, -6 - (index % 4) * 10 if index % 2 else 6 + (index % 4) * 10)
            ax1.annotate(
                pair.short_pair_label,
                (pair.shared_kmer_fraction, pair.local_identity),
                xytext=left_offset,
                textcoords="offset points",
                fontsize=10,
                bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.75},
            )
            ax2.annotate(
                pair.short_pair_label,
                (pair.length_ratio, pair.local_overlap_fraction_shorter),
                xytext=right_offset,
                textcoords="offset points",
                fontsize=10,
                bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.75},
            )

    ax1.set_xlabel("Shared k-mer fraction", fontsize=17)
    ax1.set_ylabel("Local identity", fontsize=17)
    ax1.set_title("Pairwise similarity space", fontsize=22)

    ax2.set_xscale("log")
    ax2.set_xlabel("Length ratio", fontsize=17)
    ax2.set_ylabel("Local overlap fraction of shorter monomer", fontsize=17)
    ax2.set_title("Length ratio versus overlap", fontsize=22)
    ax1.tick_params(axis="both", labelsize=14)
    ax2.tick_params(axis="both", labelsize=14)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_similarity_space")


def plot_flagged_pairs_heatmap(plt, outdir: Path, similarities: Sequence[SimilarityRecord]) -> list[Path]:
    flagged = [item for item in similarities if item.relationship != "distinct"]
    if not flagged:
        fig, ax = plt.subplots(figsize=(6, 2.5))
        ax.text(0.5, 0.5, "No non-distinct family pairs detected", ha="center", va="center", fontsize=12)
        ax.axis("off")
        fig.tight_layout()
        return save_figure(fig, outdir, "discover_flagged_pairs_heatmap")

    flagged.sort(key=lambda item: (-item.shared_kmer_fraction, -item.local_identity, item.pair_label))
    matrix = [
        [
            pair.local_identity,
            pair.local_overlap_fraction_shorter,
            pair.shared_kmer_fraction,
            pair.kmer_jaccard,
            normalize_length_ratio(pair.length_ratio),
        ]
        for pair in flagged
    ]
    fig_height = max(4.8, 0.44 * len(flagged) + 2.0)
    fig, ax = plt.subplots(figsize=(12.5, fig_height))
    image = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(5))
    ax.set_xticklabels(
        [
            "local_identity",
            "overlap_shorter",
            "shared_kmer",
            "kmer_jaccard",
            "length_ratio_norm",
        ],
        rotation=20,
        ha="right",
        fontsize=13,
    )
    ax.set_yticks(range(len(flagged)))
    ax.set_yticklabels([pair.short_pair_label for pair in flagged], fontsize=11)
    ax.set_title("Flagged family-pair metrics", fontsize=22)
    for row_index, pair in enumerate(flagged):
        values = [
            f"{pair.local_identity:.2f}",
            f"{pair.local_overlap_fraction_shorter:.2f}",
            f"{pair.shared_kmer_fraction:.2f}",
            f"{pair.kmer_jaccard:.2f}",
            f"{pair.length_ratio:.2f}",
        ]
        for col_index, value in enumerate(values):
            ax.text(col_index, row_index, value, ha="center", va="center", fontsize=9, color="#222222")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.02, pad=0.02, label="Normalized metric")
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label("Normalized metric", fontsize=14)
    fig.tight_layout()
    return save_figure(fig, outdir, "discover_flagged_pairs_heatmap")


def write_markdown_report(
    outdir: Path,
    summary: DiscoverCatalogSummary,
    families: Sequence[FamilyRecord],
    similarities: Sequence[SimilarityRecord],
    top_n: int,
) -> Path:
    top_families = sorted(families, key=lambda item: item.support_span_bp, reverse=True)[:top_n]
    flagged_pairs = [item for item in similarities if item.relationship != "distinct"]
    report_path = outdir / "discover_summary.md"
    lines = [
        "# Discover Catalog Visualization Summary",
        "",
        "## Key Statistics",
        "",
        f"- Families: {summary.family_count}",
        f"- Families with any warning: {summary.warning_count}",
        f"- Low-complexity families: {summary.low_complexity_count}",
        f"- Families with `possible_higher_order_or_partial`: {summary.relationship_warning_count}",
        f"- Non-distinct family pairs: {summary.non_distinct_pair_count}",
        f"- Median monomer length: {summary.median_length_bp:.1f} bp",
        f"- Total supporting reads: {summary.total_support_reads:,}",
        f"- Total supporting span: {summary.total_support_bp:,} bp",
        f"- Largest family by support span: {summary.top_family_id} ({summary.top_family_span_bp:,} bp)",
        "",
        "## Confidence Counts",
        "",
    ]
    lines.extend(f"- {label}: {count}" for label, count in summary.confidence_counts)
    lines.extend(["", "## Status Counts", ""])
    lines.extend(f"- {label}: {count}" for label, count in summary.status_counts)
    lines.extend(["", f"## Top {top_n} Families by Support Span", ""])
    lines.append("| family_id | length_bp | support_reads | support_span_bp | mean_identity | status | warning |")
    lines.append("|---|---:|---:|---:|---:|---|---|")
    for family in top_families:
        lines.append(
            "| "
            + " | ".join(
                [
                    family.family_id,
                    str(family.monomer_length_bp),
                    f"{family.support_read_count:,}",
                    f"{family.support_span_bp:,}",
                    f"{family.mean_identity:.4f}",
                    family.status_label,
                    family.warning or "",
                ]
            )
            + " |"
        )
    lines.extend(["", "## Flagged Family Pairs", ""])
    if not flagged_pairs:
        lines.append("No non-distinct family pairs were detected.")
    else:
        lines.append("| pair | relationship | local_identity | overlap_shorter | shared_kmer | length_ratio | notes |")
        lines.append("|---|---|---:|---:|---:|---:|---|")
        for pair in sorted(flagged_pairs, key=lambda item: (-item.shared_kmer_fraction, -item.local_identity, item.pair_label)):
            lines.append(
                "| "
                + " | ".join(
                    [
                        pair.pair_label,
                        pair.relationship,
                        f"{pair.local_identity:.4f}",
                        f"{pair.local_overlap_fraction_shorter:.4f}",
                        f"{pair.shared_kmer_fraction:.4f}",
                        f"{pair.length_ratio:.4f}",
                        pair.notes or "",
                    ]
                )
                + " |"
            )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `discover_family_abundance.pdf`",
            "- `discover_length_distribution.pdf`",
            "- `discover_abundance_vs_length.pdf`",
            "- `discover_quality_overview.pdf`",
            "- `discover_similarity_space.pdf`",
            "- `discover_flagged_pairs_heatmap.pdf`",
            "- `top_families.tsv`",
            "- `flagged_pairs.tsv`",
            "- `family_label_map.tsv`",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def status_color(status: str) -> str:
    palette = {
        "clean": "#4C78A8",
        "low_complexity": "#B279A2",
        "related": "#F58518",
        "low_complexity + related": "#E45756",
    }
    return palette[status]


def add_status_legend(ax) -> None:
    from matplotlib.lines import Line2D

    handles = [
        Line2D([0], [0], marker="s", linestyle="", color=status_color("clean"), label="clean"),
        Line2D([0], [0], marker="s", linestyle="", color=status_color("low_complexity"), label="low_complexity"),
        Line2D([0], [0], marker="s", linestyle="", color=status_color("related"), label="related"),
        Line2D([0], [0], marker="s", linestyle="", color=status_color("low_complexity + related"), label="low_complexity + related"),
    ]
    ax.legend(handles=handles, loc="best", frameon=False, fontsize=11)


def normalize_length_ratio(value: float) -> float:
    return min(1.0, math.log10(max(1.0, value)) / math.log10(10.0))


def abbreviate_family_id(family_id: str) -> str:
    if family_id.startswith("TXF") and len(family_id) >= 9:
        return f"F{family_id[-3:]}"
    return family_id


def select_similarity_annotations(
    flagged: Sequence[SimilarityRecord],
    *,
    max_labels: int,
    per_family_cap: int,
) -> list[SimilarityRecord]:
    if len(flagged) <= max_labels:
        return list(flagged)

    def score(pair: SimilarityRecord) -> tuple[float, float, int]:
        shorter = min(pair.length_a_bp, pair.length_b_bp)
        shorter_weight = 0.35 if shorter < 20 else 1.0
        relationship_score = (
            pair.shared_kmer_fraction
            * pair.local_identity
            * pair.local_overlap_fraction_shorter
            * shorter_weight
        )
        return (relationship_score, pair.shared_kmer_fraction, shorter)

    selected: list[SimilarityRecord] = []
    per_family_counts: dict[str, int] = {}
    for pair in sorted(flagged, key=score, reverse=True):
        if len(selected) >= max_labels:
            break
        family_a_count = per_family_counts.get(pair.family_a, 0)
        family_b_count = per_family_counts.get(pair.family_b, 0)
        if family_a_count >= per_family_cap or family_b_count >= per_family_cap:
            continue
        selected.append(pair)
        per_family_counts[pair.family_a] = family_a_count + 1
        per_family_counts[pair.family_b] = family_b_count + 1
    return selected


def save_figure(fig, outdir: Path, stem: str) -> list[Path]:
    output = outdir / f"{stem}.pdf"
    fig.savefig(output, bbox_inches="tight")
    fig.clf()
    try:
        import matplotlib.pyplot as plt

        plt.close(fig)
    except Exception:
        pass
    return [output]
