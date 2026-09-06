"""Four-panel diagnostic of scoring definitions and operational monomer resolution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.schema import read_table, write_table, digest_file
from benchmarks.challenge.sequence_metrics import cyclic_edit_similarity
from tandemx.discover.clustering import cluster_monomers
from tandemx.discover.mvp import CandidateRepeat


def exported_candidates(folder: Path) -> list[CandidateRepeat]:
    """Replay exported evidence, including the TSV's rounded alignment scores."""
    sequences = {header.split(";")[0].split("=", 1)[1]: sequence
                 for header, sequence in read_fasta(folder / "candidate_monomers.fa").items()}
    records = read_table(folder / "candidate_reads.tsv")
    if len(sequences) != len(records) or set(sequences) != {r["candidate_id"] for r in records}:
        raise ValueError("Candidate FASTA/TSV membership differs")
    return [CandidateRepeat(r["read_id"], r["candidate_id"], sequences[r["candidate_id"]],
                            int(r["read_start"]), int(r["read_end"]), r["strand"], int(r["period_bp"]),
                            int(r["repeat_span_bp"]), float(r["unit_count"]), float(r["score"]),
                            r["low_complexity_flag"] == "true", r["confidence"], r["warning"])
            for r in records]


def build_figure(rescore: Path, validation: Path, previous: Path, outdir: Path) -> None:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Choose a new empty figure directory")
    outdir.mkdir(parents=True, exist_ok=True)
    rescored = read_table(rescore / "rescored_metrics.tsv")
    summary = read_table(validation / "summary.tsv")
    if not all(r["successful_runs"] == r["attempted_runs"] and r["deterministic"] == "True" for r in summary):
        raise ValueError("Validation requires successful deterministic repetitions")
    dataset = "related_families_s2101"
    relative = Path("runs") / dataset / "tandemx" / "rep1" / "discover"
    folder, old_folder = validation / relative, previous / relative
    new_reads, old_reads = [p / "datasets" / dataset / "reads.fa" for p in (validation, previous)]
    if digest_file(new_reads) != digest_file(old_reads):
        raise ValueError("Previous and current related-monomer data differ")
    candidates = exported_candidates(folder)
    sources = [rescore / "rescored_metrics.tsv", validation / "summary.tsv", new_reads, old_reads,
               folder / "candidate_reads.tsv", folder / "candidate_monomers.fa",
               folder / "monomers.fa", old_folder / "monomers.fa",
               validation / "datasets" / dataset / "truth_monomers.fa"]
    table = []
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9,
                                "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.3), layout="constrained")
    colours = ("#756BB1", "#158E87")

    # A: include every completed ULTRA pilot configuration; ten reads per dataset.
    ax = axes[0, 0]
    pilot = sorted((r for r in rescored if r["tool"] == "ultra"),
                   key=lambda r: ("tuned" in r["source_run"], r["scenario"]))
    labels = []
    for i, row in enumerate(pilot):
        label = ("Tuned" if "tuned" in row["source_run"] else "Default") + ": " + row["scenario"]
        labels.append(label.replace("clean_171", "clean").replace("indel_4pct", "4% indels"))
        for j, metric in enumerate(("sequence_family_recall", "cyclic_monomer_recall")):
            value = float(row[metric])
            ax.barh(i + (j - 0.5) * 0.3, value, height=0.28, color=colours[j],
                    label=("Equal length, ungapped" if j == 0 else "Cyclic edit distance") if i == 0 else None)
            ax.text(max(0.03, value - 0.05), i + (j - 0.5) * 0.3, f"{value:.0%}", va="center",
                    ha="right" if value else "left", color="white" if value else colours[j], fontsize=8)
            table.append(dict(panel="A", item=label, metric=metric, value=value, source=row["source_run"]))
    ax.set(yticks=range(len(labels)), yticklabels=labels, xlim=(0, 1.06), xlabel="Planted-monomer recall")
    ax.invert_yaxis()
    ax.set_title("A   Small indels change the sequence endpoint\nULTRA: 10 reads per dataset, seed 1101", loc="left", weight="bold")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), frameon=False, ncol=1, fontsize=8)

    # B: all positive development scenarios, not selected for a preferred outcome.
    ax = axes[0, 1]
    trf = sorted((r for r in rescored if r["tool"] == "trf" and int(r["truth_array_count"]) > 0),
                 key=lambda r: (float(r["array_precision"]), r["scenario"]))
    for i, row in enumerate(trf):
        values = [float(row[m]) for m in ("array_precision", "base_union_precision")]
        ax.plot(values, [i, i], color="#AAAAAA", linewidth=1)
        for j, (metric, value) in enumerate(zip(("array_precision", "base_union_precision"), values)):
            ax.scatter(value, i, facecolors="none" if j == 0 else colours[j], edgecolors=colours[j],
                       s=60 if j == 0 else 22, zorder=3,
                       label=("Raw array calls" if j == 0 else "Union of bases") if i == 0 else None)
            table.append(dict(panel="B", item=row["scenario"], metric=metric, value=value, source=row["source_run"]))
    ax.set(yticks=range(len(trf)), yticklabels=[r["scenario"].replace("_", " ") for r in trf],
           xlim=(0, 1.05), xlabel="Precision (different counting units)")
    ax.tick_params(axis="y", labelsize=8)
    ax.invert_yaxis()
    ax.set_title("B   Duplicate calls differ from wrong bases\nTRF: all 13 positive scenarios, 100 reads each", loc="left", weight="bold")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, fontsize=8)

    # C: resolution sensitivity with fixed exported candidates, never rediscovery.
    ax = axes[1, 0]
    thresholds, counts = [], []
    for percent in range(90, 101):
        families, _ = cluster_monomers(candidates, 1, percent / 100, "rust")
        if percent == 95 and {f.monomer_sequence for f in families} != set(read_fasta(folder / "monomers.fa").values()):
            raise ValueError("Rounded exported-score replay differs from actual default catalogue")
        thresholds.append(percent)
        counts.append(len(families))
        table.append(dict(panel="C", item=str(percent), metric="operational_cluster_count",
                          value=len(families), source=str(folder)))
    ax.plot(thresholds, counts, "o-", color=colours[1])
    ax.axhline(3, color="#555555", linestyle="--", linewidth=0.8, label="Three planted monomers")
    ax.axvline(95, color="#888888", linestyle=":", linewidth=1)
    ax.set(xlabel="Minimum circular edit similarity (%)", ylabel="Operational monomer clusters",
           xticks=[90, 92, 94, 95, 96, 98, 100], yticks=sorted(set([2, 3, *counts])),
           ylim=(min(counts) - 0.35, max(counts) + 0.45))
    ax.set_title("C   Sequence resolution is explicit\nRelated-monomer case: seed 2101, 70 candidate arrays", loc="left", weight="bold")
    ax.legend(frameon=False, loc="upper left")

    # D: independent optimal similarities; catalogue lengths differ, no forced names.
    ax = axes[1, 1]
    truth = read_fasta(validation / "datasets" / dataset / "truth_monomers.fa")
    old = list(read_fasta(old_folder / "monomers.fa").values())
    new = list(read_fasta(folder / "monomers.fa").values())
    predictions = old + new
    normalization = Normalize(0.8, 1.0)
    labels = [f"Old {i + 1}" for i in range(len(old))] + [f"New {i + 1}" for i in range(len(new))]
    for i, (name, sequence) in enumerate(truth.items()):
        for j, predicted in enumerate(predictions):
            value = cyclic_edit_similarity(sequence, predicted)
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor=plt.cm.YlGnBu(normalization(value)), edgecolor="white"))
            ax.text(j, i, f"{value:.1%}", ha="center", va="center", color="white" if value > .92 else "black")
            table.append(dict(panel="D", item=f"{name}/{labels[j]}", metric="cyclic_edit_similarity",
                              value=value, source=str(old_folder if j < len(old) else folder)))
    ax.axvline(len(old) - .5, color="#222222", linestyle="--", linewidth=1)
    ax.set(xlim=(-.5, len(predictions) - .5), ylim=(len(truth) - .5, -.5),
           xticks=range(len(labels)), xticklabels=labels,
           yticks=range(len(truth)), yticklabels=list(truth), xlabel="Previous versus sequence-cluster catalogue")
    ax.set_title("D   Independent sequence evidence for recovery\nSame seed-2101 reads; similarity includes indels", loc="left", weight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    colourbar = fig.colorbar(plt.cm.ScalarMappable(norm=normalization, cmap="YlGnBu"), ax=ax,
                            orientation="horizontal", shrink=0.7, pad=0.08, aspect=35)
    colourbar.set_label("Cyclic edit similarity", fontsize=8)
    colourbar.set_ticks([0.8, 0.85, 0.9, 0.95, 1.0])
    colourbar.solids.set_rasterized(False)

    fig.suptitle("Evaluation definitions and monomer resolution", fontsize=15, weight="bold")
    fig.supxlabel("Development/diagnostic evidence; no held-out superiority or biological-family taxonomy claim", fontsize=9)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"sequence_audit.{extension}", dpi=180)
    plt.close(fig)
    write_table(outdir / "panel_source.tsv", table, ["panel", "item", "metric", "value", "source"])
    receipt = {"input_sha256": {str(p): digest_file(p) for p in sources},
               "script_sha256": digest_file(Path(__file__)),
               "source_commit": json.loads((validation / "environment.json").read_text())["git_head"],
               "sensitivity_replay": "rounded exported TSV scores; default catalogue equality verified",
               "limits": "first-repetition sequences; three-repeat validation digests agree; reused development/validation seeds",
               "outputs": {p.name: digest_file(p) for p in outdir.iterdir() if p.is_file()}}
    (outdir / "figure_provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("rescore", "validation", "previous", "outdir"):
        parser.add_argument("--" + argument, type=Path, required=True)
    args = parser.parse_args()
    build_figure(args.rescore, args.validation, args.previous, args.outdir)
