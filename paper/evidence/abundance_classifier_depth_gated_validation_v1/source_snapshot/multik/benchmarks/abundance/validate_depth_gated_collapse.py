"""Validate the frozen depth-gated classifier on untouched held-out simulations."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import shutil

from benchmarks.abundance.evaluate_blended_collapse import KEY_FIELDS, blend_row, confusion
from benchmarks.abundance.evaluate_depth_gated_classifier import MODEL, apply_depth_gate
from benchmarks.abundance.evaluate_multik_collapse import summarize
from benchmarks.abundance.run import validate_frozen_classifier
from benchmarks.abundance.validate_blended_collapse import heldout_acceptance
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


BASELINE = "single_k21"
MULTIK = "multik_loglinear"


def select_depth_gated_rows(
    rows: list[dict], model: dict, seeds: list[int]
) -> tuple[list[dict], list[dict], list[dict]]:
    """Pair raw estimates and apply the frozen rule without outcome-based selection."""
    paired: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        method = row.get("method")
        if method not in {BASELINE, MULTIK}:
            raise ValueError("Unexpected method in held-out paired input")
        key = tuple(row[field] for field in KEY_FIELDS)
        if method in paired[key]:
            raise ValueError("Duplicate held-out method row")
        paired[key][method] = row
    if (
        not paired
        or any(set(pair) != {BASELINE, MULTIK} for pair in paired.values())
        or {int(pair[BASELINE]["seed"]) for pair in paired.values()} != set(seeds)
        or len(rows) != len(paired) * 2
    ):
        raise ValueError("Incomplete held-out single/multi-k pairs")

    baseline_rows = []
    high_depth_rows = []
    selected_rows = []
    for pair in paired.values():
        single = dict(pair[BASELINE])
        single.update(method=BASELINE, candidate_id=BASELINE, blend_alpha=0.0)
        high = blend_row(
            pair,
            float(model["standard_depth_strategy"]["blend_alpha"]),
            float(model["standard_depth_strategy"]["decision_threshold"]),
        )
        selected = apply_depth_gate(single, high, model)
        if selected["candidate_id"] != MODEL or selected["method"] != MODEL:
            raise ValueError("Depth-gated identifier differs from frozen model")
        if set(single) != set(high) or set(single) != set(selected):
            raise ValueError("Depth-gated held-out rows have different fields")
        baseline_rows.append(single)
        high_depth_rows.append(high)
        selected_rows.append(selected)
    return baseline_rows, selected_rows, high_depth_rows


def _seed_rows(
    baseline_rows: list[dict], selected_rows: list[dict], seeds: list[int]
) -> list[dict]:
    output = []
    for seed in seeds:
        for method, rows in ((BASELINE, baseline_rows), (MODEL, selected_rows)):
            metrics = confusion([row for row in rows if int(row["seed"]) == seed])
            output.append({"seed": seed, "method": method, **metrics})
    return output


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
        raise ValueError("Require complete frozen depth-gated held-out single/multi-k evidence")

    rows = read_table(previous / "comparison_metrics.tsv", {
        *KEY_FIELDS, "period", "truth_assembly_read_ratio", "assembly_predicted_bp",
        "read_estimated_copies", "estimated_haploid_depth", "method", "outcome",
        "decision_threshold", "fit_status",
    })
    model = config["classifier_model"]
    baseline_rows, selected_rows, high_depth_rows = select_depth_gated_rows(
        rows, model, seeds
    )
    gates = config.get("heldout_rule", {}).get("primary_classifier_gates")
    if not isinstance(gates, dict):
        raise ValueError("Missing predeclared held-out classifier gates")
    acceptance, by_seed = heldout_acceptance(
        baseline_rows, selected_rows, seeds, gates
    )

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
        classifier_model=MODEL,
        frozen_development_metrics=development_metrics,
        no_heldout_fit_or_selection=True,
        scope="predeclared_depth_gated_heldout_validation_not_biological_validation",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(previous / "run_config.json", outdir / "run_config.json")
    combined = baseline_rows + selected_rows
    write_table(outdir / "comparison_metrics.tsv", combined, list(combined[0]))
    summaries = summarize(combined)
    write_table(outdir / "comparison_summary.tsv", summaries, list(summaries[0]))
    seed_rows = _seed_rows(baseline_rows, selected_rows, seeds)
    write_table(outdir / "per_seed_metrics.tsv", seed_rows, list(seed_rows[0]))
    write_table(outdir / "selected_metrics.tsv", selected_rows, list(selected_rows[0]))
    cutoff = float(model["low_depth_cutoff"])
    result = {
        "complete": len(selected_rows) * 2 == len(rows),
        "heldout": True,
        "independent_genomes": len(seeds),
        "paired_family_conditions": len(selected_rows),
        "selected_candidate": MODEL,
        "low_depth_rows": sum(
            float(row["estimated_haploid_depth"]) < cutoff for row in baseline_rows
        ),
        "standard_depth_rows": sum(
            float(row["estimated_haploid_depth"]) >= cutoff for row in baseline_rows
        ),
        "standard_depth_fallback_rows": sum(
            float(single["estimated_haploid_depth"]) >= cutoff
            and high["fit_status"].startswith("fallback_single_k21")
            for single, high in zip(baseline_rows, high_depth_rows)
        ),
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
        raise RuntimeError("Incomplete depth-gated held-out classifier validation")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run(args.previous, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
