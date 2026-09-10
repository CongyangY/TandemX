"""Evaluate the pre-specified TandemX discovery-saturation validation."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.schema import digest_file, read_table, write_table
from benchmarks.challenge.sequence_metrics import score_threshold_recovery
from benchmarks.discovery.saturation import (
    saturation_decision,
    transition_metrics,
    unique_canonical,
)


def tier(copies: int) -> str:
    if copies == 20:
        return "low"
    if copies == 80:
        return "medium"
    if copies >= 200:
        return "high"
    raise ValueError(f"Unexpected planted copy count: {copies}")


def run(run_dir: Path, outdir: Path) -> None:
    run_dir = run_dir.resolve()
    receipt = json.loads((run_dir / "run_receipt.json").read_text())
    environment = json.loads((run_dir / "environment.json").read_text())
    if not receipt.get("complete") or receipt.get("completed_runs") != 21:
        raise ValueError("Require a complete 21-run validation")
    config = json.loads(Path(environment["config"]).read_text())
    if digest_file(Path(environment["config"])) != environment["config_sha256"]:
        raise ValueError("Configuration changed after execution")
    seeds = [int(value) for value in config["seeds"]["validation"]]
    coverages = [float(value) for value in config["coverages"]]
    datasets = {int(row["seed"]): Path(row["path"]) for row in environment["datasets"]}
    run_rows = read_table(run_dir / "run_summary.tsv")
    outdir.mkdir(parents=True, exist_ok=False)

    metrics = []
    family_rows = []
    catalogues: dict[tuple[int, float], list[str]] = {}
    recovered: dict[tuple[int, float], set[str]] = {}
    for row in run_rows:
        seed, coverage = int(row["seed"]), float(row["coverage"])
        if row["status"] != "ok" or seed not in datasets:
            raise ValueError("Incomplete or unknown run in summary")
        folder = Path(row["run_path"]) / "discover"
        summary = json.loads((folder / "discovery_summary.json").read_text())
        for name, expected in summary["output_sha256"].items():
            if digest_file(folder / name) != expected:
                raise ValueError(f"Discovery output changed: {folder / name}")
        candidates = list(read_fasta(folder / "candidate_monomers.fa").values())
        families = list(read_fasta(folder / "monomers.fa").values())
        unique_candidates = unique_canonical(candidates)
        unique_families = unique_canonical(families)
        if int(summary["family_count"]) != len(families) or len(families) != len(unique_families):
            raise ValueError(f"Family count or canonical uniqueness differs: {folder}")
        dataset = datasets[seed]
        founders = read_fasta(dataset / "genome/catalogue.fa")
        truth = {
            truth_row["family_id"]: {
                "copies": int(truth_row["copies"]),
                "period": int(truth_row["period"]),
            }
            for truth_row in read_table(dataset / "genome/truth_copy_number.tsv")
        }
        if set(founders) != set(truth):
            raise ValueError(f"Founder catalogue differs from truth rows: {dataset}")
        recovery_metrics, details = score_threshold_recovery(unique_families, founders, 0.90)
        recovered_ids = {item["truth_id"] for item in details if item["recovered"] == 1}
        recovered[(seed, coverage)] = recovered_ids
        catalogues[(seed, coverage)] = unique_families
        tier_recall = {}
        for label in ("low", "medium", "high"):
            identifiers = {name for name, item in truth.items() if tier(item["copies"]) == label}
            tier_recall[label] = len(identifiers & recovered_ids) / len(identifiers)
        metrics.append(
            {
                "seed": seed,
                "coverage": coverage,
                "actual_bases": int(row["total_bases"]),
                "read_count": int(row["read_count"]),
                "raw_candidate_count": int(summary["candidate_count"]),
                "unique_candidate_count": len(unique_candidates),
                "operational_family_count": len(unique_families),
                "recovered_family_count": recovery_metrics["recovered_family_count"],
                "overall_recall": recovery_metrics["cyclic_monomer_recall"],
                "low_abundance_recall": tier_recall["low"],
                "medium_abundance_recall": tier_recall["medium"],
                "high_abundance_recall": tier_recall["high"],
                "homologous_family_fraction": recovery_metrics["homologous_consensus_fraction"],
                "unmatched_operational_families": recovery_metrics[
                    "unmatched_distinct_consensus_count"
                ],
                "runtime_seconds": float(row["runtime_seconds"]),
                "peak_rss_mib": float(row["peak_rss_mib"]),
            }
        )
        for item in details:
            name = item["truth_id"]
            family_rows.append(
                {
                    "seed": seed,
                    "coverage": coverage,
                    "truth_id": name,
                    "abundance_tier": tier(truth[name]["copies"]),
                    "copies": truth[name]["copies"],
                    "period": truth[name]["period"],
                    "recovered": item["recovered"],
                }
            )
    metrics.sort(key=lambda row: (row["seed"], row["coverage"]))
    family_rows.sort(key=lambda row: (row["seed"], row["coverage"], row["truth_id"]))
    write_table(outdir / "depth_metrics.tsv", metrics, list(metrics[0]))
    write_table(outdir / "family_recovery.tsv", family_rows, list(family_rows[0]))

    transitions = []
    for seed in seeds:
        for low, high in zip(coverages, coverages[1:]):
            observed = transition_metrics(
                catalogues[(seed, low)], catalogues[(seed, high)], low, high, 0.90
            )
            new_truth = recovered[(seed, high)] - recovered[(seed, low)]
            lost_truth = recovered[(seed, low)] - recovered[(seed, high)]
            transitions.append(
                {
                    "seed": seed,
                    "lower_depth": low,
                    "higher_depth": high,
                    **observed,
                    "new_truth_family_count": len(new_truth),
                    "lost_truth_family_count": len(lost_truth),
                    "new_truth_families_per_added_x": len(new_truth) / (high - low),
                    "new_truth_family_ids": ",".join(sorted(new_truth)) or "none",
                    "lost_truth_family_ids": ",".join(sorted(lost_truth)) or "none",
                }
            )
    transitions.sort(key=lambda row: (row["seed"], row["higher_depth"]))
    write_table(outdir / "transition_metrics.tsv", transitions, list(transitions[0]))
    decision = saturation_decision(transitions, seeds, coverages)
    decision.update(
        seeds=seeds,
        coverages=coverages,
        family_match_threshold=0.90,
        main_text_eligible=bool(decision["saturation_reached"]),
        scientific_interpretation="pending_plot_and_cross_metric_review",
    )
    (outdir / "saturation_decision.json").write_text(json.dumps(decision, indent=2) + "\n")

    if len(metrics) != 21 or len(family_rows) != 21 * 55 or len(transitions) != 18:
        raise RuntimeError("Evaluation cardinality differs from protocol")
    (outdir / "evaluation_receipt.json").write_text(
        json.dumps(
            {
                "complete": True,
                "depth_rows": len(metrics),
                "family_rows": len(family_rows),
                "transition_rows": len(transitions),
                "saturation_reached": decision["saturation_reached"],
                "saturation_depth": decision["saturation_depth"],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    run(args.run_dir, args.outdir)
