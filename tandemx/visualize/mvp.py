"""Basic static visualization MVP for TandemX toy outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class VisualizeConfig:
    copy_number: Path
    comparison: Path | None
    probes: Path | None
    fish: Path | None
    outdir: Path


def render_static_plots(config: VisualizeConfig) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config.outdir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    copy_rows = read_tsv(config.copy_number)
    comparison_rows = read_tsv(config.comparison) if config.comparison else []
    probe_rows = read_tsv(config.probes) if config.probes else []
    fish_rows = read_tsv(config.fish) if config.fish else []

    outputs.extend(plot_catalogue_summary(plt, config.outdir, copy_rows, comparison_rows, probe_rows))
    if comparison_rows:
        outputs.extend(plot_assembly_vs_read(plt, config.outdir, comparison_rows))
    if fish_rows:
        outputs.extend(plot_in_silico_fish(plt, config.outdir, fish_rows))
    return outputs


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def plot_catalogue_summary(plt, outdir: Path, copy_rows: Sequence[dict[str, str]], comparison_rows: Sequence[dict[str, str]], probe_rows: Sequence[dict[str, str]]) -> list[Path]:
    family_ids = [row["family_id"] for row in copy_rows]
    read_bp = [float(row.get("estimated_bp", 0.0)) for row in copy_rows]
    assembly_by_family = {row["family_id"]: float(row.get("assembly_estimated_bp", 0.0)) for row in comparison_rows}
    probe_by_family: dict[str, float] = {}
    for row in probe_rows:
        family_id = row["family_id"]
        probe_by_family[family_id] = max(probe_by_family.get(family_id, 0.0), float(row.get("probe_score", 0.0)))
    assembly_bp = [assembly_by_family.get(family_id, 0.0) for family_id in family_ids]
    probe_scores = [probe_by_family.get(family_id, 0.0) for family_id in family_ids]

    fig, ax1 = plt.subplots(figsize=(7, 4))
    x = list(range(len(family_ids)))
    width = 0.35
    ax1.bar([value - width / 2 for value in x], read_bp, width=width, label="read bp", color="#4C78A8")
    ax1.bar([value + width / 2 for value in x], assembly_bp, width=width, label="assembly bp", color="#F58518")
    ax1.set_ylabel("Estimated bp")
    ax1.set_xticks(x)
    ax1.set_xticklabels(family_ids, rotation=30, ha="right")
    ax1.legend(loc="upper left")
    ax2 = ax1.twinx()
    ax2.plot(x, probe_scores, marker="o", color="#54A24B", label="probe score")
    ax2.set_ylabel("Probe score")
    ax2.set_ylim(0, max([1.0, *probe_scores]))
    ax2.legend(loc="upper right")
    fig.tight_layout()
    return save_figure(fig, outdir, "catalogue_summary")


def plot_assembly_vs_read(plt, outdir: Path, rows: Sequence[dict[str, str]]) -> list[Path]:
    fig, ax = plt.subplots(figsize=(5, 5))
    colors = ["#E45756" if row.get("status") == "possible_collapse" else "#4C78A8" for row in rows]
    x = [float(row.get("read_estimated_bp", 0.0)) for row in rows]
    y = [float(row.get("assembly_estimated_bp", 0.0)) for row in rows]
    ax.scatter(x, y, c=colors)
    for row, x_value, y_value in zip(rows, x, y):
        ax.text(x_value, y_value, row["family_id"], fontsize=8)
    max_value = max([1.0, *x, *y])
    ax.plot([0, max_value], [0, max_value], linestyle="--", color="#777777", linewidth=1)
    ax.set_xlabel("Read-estimated bp")
    ax.set_ylabel("Assembly-estimated bp")
    ax.set_title("Assembly vs read repeat abundance")
    fig.tight_layout()
    return save_figure(fig, outdir, "assembly_vs_read")


def plot_in_silico_fish(plt, outdir: Path, rows: Sequence[dict[str, str]]) -> list[Path]:
    chroms = sorted({row["chrom"] for row in rows})
    chrom_index = {chrom: index for index, chrom in enumerate(chroms)}
    fig, ax = plt.subplots(figsize=(7, max(2.5, len(chroms) * 0.8)))
    for chrom in chroms:
        chrom_rows = [row for row in rows if row["chrom"] == chrom]
        max_end = max(int(row["end"]) for row in chrom_rows)
        y = chrom_index[chrom]
        ax.hlines(y, 0, max_end, color="#CCCCCC", linewidth=6)
    for row in rows:
        y = chrom_index[row["chrom"]]
        start = int(row["start"])
        end = int(row["end"])
        signal = float(row.get("predicted_signal", 0.0))
        ax.hlines(y, start, end, color="#E45756", linewidth=max(2, 8 * min(1.0, signal)))
    ax.set_yticks(list(chrom_index.values()))
    ax.set_yticklabels(chroms)
    ax.set_xlabel("Position (bp)")
    ax.set_title("Toy in silico FISH predicted signal")
    fig.tight_layout()
    return save_figure(fig, outdir, "in_silico_fish")


def save_figure(fig, outdir: Path, stem: str) -> list[Path]:
    outputs = [outdir / f"{stem}.svg", outdir / f"{stem}.pdf"]
    for path in outputs:
        fig.savefig(path)
    fig.clf()
    return outputs
