"""Final allowed phase revision, reusing immutable first-round development inputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import time
from dataclasses import asdict
from pathlib import Path

from benchmarks.challenge.local_phase import LocalPhaseConfig, LocalPhasePrototype
from benchmarks.challenge.phase_gate_metrics import score_predictions


def summarize(rows: list[dict]) -> dict:
    result = {}
    for model in sorted({row["model"] for row in rows}):
        selected = [row for row in rows if row["model"] == model]
        positive = [row for row in selected if row["absolute_relative_error"] is not None]
        tp = sum(row["binary_call"] and row["binary_truth"] for row in selected)
        fp = sum(row["binary_call"] and not row["binary_truth"] for row in selected)
        pos = sum(row["binary_truth"] for row in selected)
        indel = [row["recall"] for row in selected if row["condition"].startswith("indels") and row["recall"] is not None]
        result[model] = {
            "family_conditions": len(selected),
            "MARE": statistics.mean(row["absolute_relative_error"] for row in positive),
            "MARE_by_source_seed": {seed: statistics.mean(row["absolute_relative_error"] for row in positive if row["seed"] == seed) for seed in sorted({row["seed"] for row in positive})},
            "MARE_by_condition": {condition: statistics.mean(row["absolute_relative_error"] for row in positive if row["condition"] == condition) for condition in sorted({row["condition"] for row in positive})},
            "indel_recall": statistics.mean(indel) if indel else None,
            "clean_MARE": statistics.mean(row["absolute_relative_error"] for row in positive if row["condition"] == "clean"),
            "background_excess_bp": sum(row["excess_abundance_bp"] for row in selected if row["condition"] == "background_homology"),
            "sensitivity": tp / pos,
            "FPR": fp / (len(selected)-pos),
            "precision": tp/(tp+fp) if tp+fp else None,
            "mean_case_seconds": statistics.mean(row["seconds"] for row in selected),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round-one", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    config = LocalPhaseConfig(bridge_accepted_short_gaps=True)
    original = args.round_one / "rows.json"
    rows = json.loads(original.read_text())
    protocol = {
        "split": "development", "algorithm_round": 2,
        "only_change": "bridge short gaps on accepted chains; all other thresholds unchanged",
        "config": asdict(config), "round_one_rows_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
        "sources": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (
            Path(__file__), Path("benchmarks/challenge/local_phase.py"),
            Path("benchmarks/challenge/phase_gate_metrics.py"), Path("tandemx/quantify/mvp.py"),
            Path("tandemx/utils/kmers.py"), Path("benchmarks/challenge/context_prototype.py"))},
    }
    (args.outdir / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    for row in rows:
        if row["model"] == "A3":
            row["model"] = "A3_round_one"
    max_rounding_difference = 0
    for root in sorted(args.round_one.glob("dev-*")):
        case = json.loads((root / "truth.json").read_text())
        seed = root.name.split("_", 1)[0]
        condition, coverage = case["condition"], case["coverage"]
        group = [row for row in rows if row["seed"] == seed and row["condition"] == condition and row["coverage"] == coverage]
        with (root / "A0/copy_number.tsv").open() as handle:
            native = {row["family_id"]: float(row["estimated_bp"]) for row in csv.DictReader(handle, delimiter="\t")}
        for row in group:
            if row["model"] == "A0":
                bp = native[row["family"]]
                max_rounding_difference = max(max_rounding_difference, abs(bp-row["predicted_unique_bp"]))
                row["predicted_unique_bp"] = bp
                truth = row["read_truth_bp"]
                row["absolute_relative_error"] = abs(bp-truth)/truth if truth else None
                row["excess_abundance_bp"] = max(0,bp-truth)
        start = time.perf_counter()
        model = LocalPhasePrototype(case["catalogue"], config)
        predictions = {rid: model.intervals(seq) for rid, seq in case["reads"]}
        seconds = time.perf_counter() - start
        score = score_predictions(predictions, case["truth_intervals"], {rid: len(seq) for rid, seq in case["reads"]}, set(case["catalogue"]))
        for family, values in score["per_family"].items():
            row = dict(seed=seed, condition=condition, coverage=coverage, model="A3_round_two", family=family, seconds=seconds, **values)
            group.append(row)
            rows.append(row)
        for row in group:
            family = row["family"]
            assembled = case["source_truth_bp"][family] * (0.3 if family == "F1" else 1)
            estimate = row["predicted_unique_bp"] / (case["read_bases"] / case["genome_size"])
            row.update(binary_truth=family == "F1", binary_call=assembled/estimate < .6 if estimate else False,
                       excess_abundance_bp=max(0,row["predicted_unique_bp"]-row["read_truth_bp"]))
        (args.outdir / f"{root.name}.json").write_text(json.dumps(score, indent=2) + "\n")
    summary = summarize(rows)
    a0,a1,a2,a3=(summary[name] for name in ("A0","A1","A2","A3_round_two"))
    gates = {
        "MARE_15pct_better_than_A0": a3["MARE"] <= .85*a0["MARE"],
        "MARE_15pct_better_than_A1": a3["MARE"] <= .85*a1["MARE"],
        "background_excess_halved_vs_A0": a3["background_excess_bp"] <= .5*a0["background_excess_bp"],
        "indel_recall_better_than_A2": a3["indel_recall"] > a2["indel_recall"],
        "clean_MARE_no_more_than_002_worse_than_A0": a3["clean_MARE"] <= a0["clean_MARE"]+.02,
        "binary_FPR_at_most_002": a3["FPR"] <= .02,
        "preferred_runtime_at_most_twice_A0": a3["mean_case_seconds"] <= 2*a0["mean_case_seconds"],
    }
    result = {"summary": summary, "gates": gates, "A0_native_rounding_max_difference_bp": max_rounding_difference,
              "heldout_status": "not_generated_or_inspected", "decision": "stop_accuracy_branch" if not all(gates.values()) else "eligible_for_heldout_preregistration"}
    (args.outdir / "rows.json").write_text(json.dumps(rows, indent=2) + "\n")
    (args.outdir / "assessment.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
