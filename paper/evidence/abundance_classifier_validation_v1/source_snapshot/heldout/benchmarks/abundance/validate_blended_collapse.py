"""Validate one frozen blend on untouched held-out abundance simulations."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil

from benchmarks.abundance.evaluate_blended_collapse import KEY_FIELDS, blend_row, confusion
from benchmarks.abundance.evaluate_multik_collapse import summarize
from benchmarks.abundance.run import validate_frozen_classifier
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


BASELINE = "single_k21"
MULTIK = "multik_loglinear"


def heldout_acceptance(
    baseline_rows: list[dict],
    selected_rows: list[dict],
    seeds: list[int],
    gates: dict,
) -> tuple[dict, dict[int, dict]]:
    baseline = confusion(baseline_rows)
    selected = confusion(selected_rows)
    by_seed = {}
    for seed in seeds:
        seed_baseline = confusion([
            row for row in baseline_rows if int(row["seed"]) == seed
        ])
        seed_selected = confusion([
            row for row in selected_rows if int(row["seed"]) == seed
        ])
        if any(
            not isinstance(metrics[name], (int, float)) or not math.isfinite(metrics[name])
            for metrics in (seed_baseline, seed_selected)
            for name in ("sensitivity", "false_positive_rate", "precision")
        ):
            raise ValueError("Held-out seed metrics must be finite")
        by_seed[seed] = {
            "baseline": seed_baseline,
            "selected": seed_selected,
            "sensitivity_delta": seed_selected["sensitivity"] - seed_baseline["sensitivity"],
            "false_positive_rate_delta": (
                seed_selected["false_positive_rate"] - seed_baseline["false_positive_rate"]
            ),
            "precision_delta": seed_selected["precision"] - seed_baseline["precision"],
        }
    acceptance = {
        "full_sensitivity_delta": selected["sensitivity"] - baseline["sensitivity"],
        "full_false_positive_rate_delta": (
            selected["false_positive_rate"] - baseline["false_positive_rate"]
        ),
        "full_precision_delta": selected["precision"] - baseline["precision"],
        "minimum_seed_sensitivity_delta": min(
            value["sensitivity_delta"] for value in by_seed.values()
        ),
        "maximum_seed_false_positive_rate_delta": max(
            value["false_positive_rate_delta"] for value in by_seed.values()
        ),
        "minimum_seed_precision_delta": min(
            value["precision_delta"] for value in by_seed.values()
        ),
    }
    expected_gates = {
        "full_sensitivity_min_delta_vs_single_k21",
        "full_false_positive_rate_max_delta_vs_single_k21",
        "full_precision_min_delta_vs_single_k21",
        "minimum_seed_sensitivity_min_delta_vs_single_k21",
        "maximum_seed_false_positive_rate_max_delta_vs_single_k21",
        "minimum_seed_precision_min_delta_vs_single_k21",
    }
    if set(gates) != expected_gates:
        raise ValueError("Invalid held-out classifier gates")
    acceptance["passed"] = (
        acceptance["full_sensitivity_delta"]
            >= float(gates["full_sensitivity_min_delta_vs_single_k21"])
        and acceptance["full_false_positive_rate_delta"]
            <= float(gates["full_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["full_precision_delta"]
            >= float(gates["full_precision_min_delta_vs_single_k21"])
        and acceptance["minimum_seed_sensitivity_delta"]
            >= float(gates["minimum_seed_sensitivity_min_delta_vs_single_k21"])
        and acceptance["maximum_seed_false_positive_rate_delta"]
            <= float(gates["maximum_seed_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["minimum_seed_precision_delta"]
            >= float(gates["minimum_seed_precision_min_delta_vs_single_k21"])
    )
    return acceptance, by_seed


def _seed_rows(baseline_rows: list[dict], selected_rows: list[dict], seeds: list[int]) -> list[dict]:
    output = []
    for seed in seeds:
        for method, rows in ((BASELINE, baseline_rows),
                             (selected_rows[0]["candidate_id"], selected_rows)):
            metrics = confusion([row for row in rows if int(row["seed"]) == seed])
            output.append({"seed": seed, "method": method, **metrics})
    return output


def comparison_output_rows(
    baseline_rows: list[dict], selected_rows: list[dict]
) -> list[dict]:
    """Give baseline and selected rows one stable TSV schema."""
    normalized_baseline = [
        {**row, "candidate_id": BASELINE, "blend_alpha": 0.0}
        for row in baseline_rows
    ]
    combined = normalized_baseline + selected_rows
    fields = set(combined[0])
    if any(set(row) != fields for row in combined):
        raise ValueError("Baseline and selected held-out rows have different fields")
    return combined


def run(previous: Path, outdir: Path) -> None:
    previous = previous.resolve()
    prior_validation = json.loads((previous / "validation.json").read_text())
    prior_environment = json.loads((previous / "environment.json").read_text())
    config = json.loads((previous / "run_config.json").read_text())
    development_metrics = validate_frozen_classifier(config, "heldout")
    seeds = [int(seed) for seed in config.get("seeds", {}).get("heldout", [])]
    if (
        prior_validation.get("complete") is not True
        or prior_validation.get("split") != "heldout"
        or prior_validation.get("evaluation_mode")
            != "raw_single_multik_for_predeclared_blend_grid"
        or set(prior_validation.get("method_confusion", {})) != {BASELINE, MULTIK}
        or prior_environment.get("split") != "heldout"
        or digest_file(previous / "run_config.json")
            != prior_environment.get("previous_config_sha256")
        or digest_file(previous / "comparison_metrics.tsv")
            != prior_validation.get("comparison_metrics_sha256")
        or development_metrics is None
        or len(seeds) != 3 or len(set(seeds)) != 3
    ):
        raise ValueError("Require complete frozen-model held-out single/multi-k evidence")
    rows = read_table(previous / "comparison_metrics.tsv", {
        *KEY_FIELDS, "period", "truth_assembly_read_ratio", "assembly_predicted_bp",
        "read_estimated_copies", "method", "outcome", "decision_threshold",
    })
    paired: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        if row["method"] not in {BASELINE, MULTIK}:
            raise ValueError("Unexpected method in held-out paired input")
        key = tuple(row[field] for field in KEY_FIELDS)
        if row["method"] in paired[key]:
            raise ValueError("Duplicate held-out method row")
        paired[key][row["method"]] = row
    if (
        not paired
        or any(set(pair) != {BASELINE, MULTIK} for pair in paired.values())
        or {int(pair[BASELINE]["seed"]) for pair in paired.values()} != set(seeds)
        or len(rows) != len(paired) * 2
    ):
        raise ValueError("Incomplete held-out single/multi-k pairs")
    model = config["classifier_model"]
    baseline_rows = [pair[BASELINE] for pair in paired.values()]
    selected_rows = [
        blend_row(pair, float(model["blend_alpha"]), float(model["decision_threshold"]))
        for pair in paired.values()
    ]
    selected_id = selected_rows[0]["candidate_id"]
    expected_id = f"blend_a{float(model['blend_alpha']):g}_t{float(model['decision_threshold']):g}"
    if selected_id != expected_id:
        raise ValueError("Held-out classifier identifier differs from frozen model")
    gates = config.get("heldout_rule", {}).get("primary_classifier_gates")
    if not isinstance(gates, dict):
        raise ValueError("Missing predeclared held-out classifier gates")
    acceptance, by_seed = heldout_acceptance(baseline_rows, selected_rows, seeds, gates)

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir / "source_snapshot"
    provenance = source_manifest(root, snapshot)
    relative = Path(__file__).resolve().relative_to(root)
    target = snapshot / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(__file__), target)
    provenance.update(
        benchmark_source_sha256={str(relative): digest_file(target)},
        previous_result=str(previous),
        previous_validation_sha256=digest_file(previous / "validation.json"),
        previous_environment_sha256=digest_file(previous / "environment.json"),
        previous_config_sha256=digest_file(previous / "run_config.json"),
        split="heldout",
        heldout_seeds=seeds,
        classifier_model=selected_id,
        frozen_development_metrics=development_metrics,
        scope="predeclared_heldout_classifier_validation_not_biological_validation",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(previous / "run_config.json", outdir / "run_config.json")
    combined = comparison_output_rows(baseline_rows, selected_rows)
    write_table(outdir / "comparison_metrics.tsv", combined, list(combined[0]))
    summaries = summarize(combined)
    write_table(outdir / "comparison_summary.tsv", summaries, list(summaries[0]))
    seed_rows = _seed_rows(baseline_rows, selected_rows, seeds)
    write_table(outdir / "per_seed_metrics.tsv", seed_rows, list(seed_rows[0]))
    write_table(outdir / "selected_metrics.tsv", selected_rows, list(selected_rows[0]))
    result = {
        "complete": len(selected_rows) == len(paired),
        "heldout": True,
        "independent_genomes": len(seeds),
        "paired_family_conditions": len(paired),
        "selected_candidate": selected_id,
        "fallback_rows": sum(row["fit_status"].startswith("fallback_single_k21")
                             for row in selected_rows),
        "baseline_confusion": confusion(baseline_rows),
        "selected_confusion": confusion(selected_rows),
        "per_seed": by_seed,
        "acceptance": acceptance,
        "comparison_metrics_sha256": digest_file(outdir / "comparison_metrics.tsv"),
        "selected_metrics_sha256": digest_file(outdir / "selected_metrics.tsv"),
        "scientific_acceptance": "heldout_gate_reported_not_assumed_from_execution_success",
    }
    (outdir / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["complete"]:
        raise RuntimeError("Incomplete held-out classifier validation")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run(args.previous, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
