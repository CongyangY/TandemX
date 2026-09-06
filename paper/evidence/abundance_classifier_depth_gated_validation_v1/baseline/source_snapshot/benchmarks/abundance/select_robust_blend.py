"""Select a seed-robust blend from a frozen development candidate grid."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil

from benchmarks.abundance.evaluate_blended_collapse import KEY_FIELDS, confusion
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


SOURCE_FILES = {
    "validation_sha256": "validation.json",
    "candidate_metrics_sha256": "candidate_metrics.tsv",
    "candidate_summary_sha256": "candidate_summary.tsv",
    "selection_sha256": "selection.tsv",
    "environment_sha256": "environment.json",
    "run_config_sha256": "run_config.json",
}


def _candidate_id(alpha: float, threshold: float) -> str:
    return f"blend_a{alpha:g}_t{threshold:g}"


def _validate_config(config: dict) -> tuple[list[int], list[int], dict, list[str]]:
    development = config.get("development_seeds")
    heldout = config.get("reserved_heldout_seeds")
    source = config.get("candidate_source")
    grid = config.get("candidate_grid")
    rule = config.get("selection_rule")
    if (
        not isinstance(development, list) or len(development) != 3
        or len(set(development)) != 3
        or not all(isinstance(seed, int) and not isinstance(seed, bool) for seed in development)
        or not isinstance(heldout, list) or len(heldout) != 3
        or len(set(heldout)) != 3 or set(development) & set(heldout)
        or not all(isinstance(seed, int) and not isinstance(seed, bool) for seed in heldout)
        or not isinstance(source, dict) or set(source) != {"result", *SOURCE_FILES}
        or any(not isinstance(source[name], str) for name in source)
        or any(len(source[name]) != 64 or any(c not in "0123456789abcdef" for c in source[name])
               for name in SOURCE_FILES)
        or not isinstance(grid, dict) or set(grid) != {"alpha", "decision_threshold"}
    ):
        raise ValueError("Invalid robust-selection provenance or seed groups")
    alphas = grid["alpha"]
    thresholds = grid["decision_threshold"]
    if (
        alphas != [0, 0.25, 0.5, 0.75, 1]
        or thresholds != [0.45, 0.5, 0.55, 0.6]
        or not isinstance(rule, dict)
        or rule.get("method") != "seed_robust_minimax"
        or rule.get("baseline_candidate") != "blend_a0_t0.6"
        or rule.get("selection_objective") != (
            "maximize_minimum_seed_sensitivity_delta_then_full_sensitivity_then_"
            "minimize_maximum_seed_fpr_delta_then_maximize_minimum_seed_precision_"
            "delta_then_minimize_alpha_then_threshold"
        )
        or set(rule.get("constraints", {})) != {
            "per_seed_false_positive_rate_max_delta", "per_seed_precision_min_delta"
        }
        or set(rule.get("acceptance_gates", {})) != {
            "minimum_seed_sensitivity_delta", "maximum_seed_false_positive_rate_delta",
            "minimum_seed_precision_delta", "full_sensitivity_delta",
            "full_false_positive_rate_delta", "full_precision_delta",
        }
    ):
        raise ValueError("Invalid robust-selection candidate grid or rule")
    candidates = [_candidate_id(float(alpha), float(threshold))
                  for alpha in alphas for threshold in thresholds]
    return development, heldout, rule, candidates


def _finite(metrics: dict) -> bool:
    return all(isinstance(metrics[name], (int, float)) and math.isfinite(metrics[name])
               for name in ("sensitivity", "false_positive_rate", "precision"))


def select_robust_candidate(
    candidates: dict[str, list[dict]],
    seeds: list[int],
    rule: dict,
) -> tuple[str, list[dict], dict]:
    baseline_id = rule["baseline_candidate"]
    if baseline_id not in candidates:
        raise ValueError("Missing declared baseline candidate")
    baseline = candidates[baseline_id]
    baseline_full = confusion(baseline)
    baseline_seed = {
        seed: confusion([row for row in baseline if int(row["seed"]) == seed])
        for seed in seeds
    }
    if not _finite(baseline_full) or not all(_finite(value) for value in baseline_seed.values()):
        raise ValueError("Baseline metrics require finite sensitivity, FPR and precision")
    constraints = rule["constraints"]
    summaries = []
    selectable = []
    for candidate_id, rows in candidates.items():
        first = rows[0]
        full = confusion(rows)
        per_seed = {
            seed: confusion([row for row in rows if int(row["seed"]) == seed])
            for seed in seeds
        }
        if not _finite(full) or not all(_finite(value) for value in per_seed.values()):
            raise ValueError(f"Candidate has non-finite metrics: {candidate_id}")
        sensitivity_deltas = [
            per_seed[seed]["sensitivity"] - baseline_seed[seed]["sensitivity"]
            for seed in seeds
        ]
        fpr_deltas = [
            per_seed[seed]["false_positive_rate"]
            - baseline_seed[seed]["false_positive_rate"]
            for seed in seeds
        ]
        precision_deltas = [
            per_seed[seed]["precision"] - baseline_seed[seed]["precision"]
            for seed in seeds
        ]
        eligible = (
            full["unavailable"] == 0
            and all(value["unavailable"] == 0 for value in per_seed.values())
            and max(fpr_deltas)
                <= float(constraints["per_seed_false_positive_rate_max_delta"]) + 1e-15
            and min(precision_deltas)
                >= float(constraints["per_seed_precision_min_delta"]) - 1e-15
        )
        summary = {
            "candidate_id": candidate_id,
            "blend_alpha": float(first["blend_alpha"]),
            "decision_threshold": float(first["decision_threshold"]),
            "minimum_seed_sensitivity_delta": min(sensitivity_deltas),
            "maximum_seed_false_positive_rate_delta": max(fpr_deltas),
            "minimum_seed_precision_delta": min(precision_deltas),
            "full_sensitivity": full["sensitivity"],
            "full_false_positive_rate": full["false_positive_rate"],
            "full_precision": full["precision"],
            "full_sensitivity_delta": full["sensitivity"] - baseline_full["sensitivity"],
            "full_false_positive_rate_delta": (
                full["false_positive_rate"] - baseline_full["false_positive_rate"]
            ),
            "full_precision_delta": full["precision"] - baseline_full["precision"],
            "eligible": eligible,
        }
        summaries.append(summary)
        if eligible:
            selectable.append(summary)
    if not selectable:
        raise ValueError("No candidate satisfies every seed-level constraint")
    selected = max(selectable, key=lambda row: (
        row["minimum_seed_sensitivity_delta"],
        row["full_sensitivity"],
        -row["maximum_seed_false_positive_rate_delta"],
        row["minimum_seed_precision_delta"],
        -row["blend_alpha"],
        -row["decision_threshold"],
    ))
    return selected["candidate_id"], summaries, {
        "baseline_full": baseline_full,
        "baseline_by_seed": baseline_seed,
    }


def run(previous: Path, config_path: Path, outdir: Path) -> None:
    previous = previous.resolve()
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text())
    development, heldout, rule, expected_candidates = _validate_config(config)
    source = config["candidate_source"]
    if Path(source["result"]).resolve() != previous:
        raise ValueError("Candidate source differs from frozen selection config")
    for digest_name, filename in SOURCE_FILES.items():
        if digest_file(previous / filename) != source[digest_name]:
            raise ValueError(f"Frozen candidate source hash differs: {filename}")
    prior_validation = json.loads((previous / "validation.json").read_text())
    prior_environment = json.loads((previous / "environment.json").read_text())
    if (
        prior_validation.get("complete") is not True
        or prior_validation.get("development_only") is not True
        or prior_validation.get("candidate_models") != len(expected_candidates)
        or prior_validation.get("acceptance", {}).get("passed") is not False
        or prior_validation.get("acceptance", {}).get(
            "leave_one_seed_out_selection_consistent") is not False
        or prior_environment.get("development_seeds") != development
        or prior_environment.get("reserved_heldout_seeds") != heldout
    ):
        raise ValueError("Require the frozen failed v1 development selection")
    rows = read_table(previous / "candidate_metrics.tsv", {
        *KEY_FIELDS, "candidate_id", "blend_alpha", "decision_threshold", "outcome"
    })
    candidates: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        candidates[row["candidate_id"]].append(row)
    if set(candidates) != set(expected_candidates):
        raise ValueError("Candidate IDs differ from the frozen 5-by-4 grid")
    expected_rows = int(prior_validation["paired_family_conditions"])
    for candidate_id, group in candidates.items():
        if (
            len(group) != expected_rows
            or {int(row["seed"]) for row in group} != set(development)
            or len({tuple(row[field] for field in KEY_FIELDS) for row in group}) != expected_rows
        ):
            raise ValueError(f"Incomplete candidate rows: {candidate_id}")
        counts = {seed: sum(int(row["seed"]) == seed for row in group)
                  for seed in development}
        if len(set(counts.values())) != 1 or sum(counts.values()) != expected_rows:
            raise ValueError(f"Unbalanced candidate seeds: {candidate_id}")
    selected_id, summaries, baseline = select_robust_candidate(candidates, development, rule)
    selected_summary = next(row for row in summaries if row["candidate_id"] == selected_id)
    gates = rule["acceptance_gates"]
    acceptance = {
        name: selected_summary[name] for name in gates
    }
    acceptance["passed"] = (
        acceptance["minimum_seed_sensitivity_delta"]
            >= float(gates["minimum_seed_sensitivity_delta"])
        and acceptance["maximum_seed_false_positive_rate_delta"]
            <= float(gates["maximum_seed_false_positive_rate_delta"])
        and acceptance["minimum_seed_precision_delta"]
            >= float(gates["minimum_seed_precision_delta"])
        and acceptance["full_sensitivity_delta"] >= float(gates["full_sensitivity_delta"])
        and acceptance["full_false_positive_rate_delta"]
            <= float(gates["full_false_positive_rate_delta"])
        and acceptance["full_precision_delta"] >= float(gates["full_precision_delta"])
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
        config_sha256=digest_file(config_path),
        candidate_source=str(previous),
        candidate_source_hashes={name: source[name] for name in SOURCE_FILES},
        development_seeds=development,
        reserved_heldout_seeds=heldout,
        scope="post_v1_seed_robust_development_refinement_not_independent_validation",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(config_path, outdir / "run_config.json")
    write_table(outdir / "candidate_robust_summary.tsv", summaries, list(summaries[0]))
    selection = [{**selected_summary, "selection_method": rule["method"]}]
    write_table(outdir / "selection.tsv", selection, list(selection[0]))
    selected_rows = candidates[selected_id]
    selected_by_seed = {
        seed: confusion([row for row in selected_rows if int(row["seed"]) == seed])
        for seed in development
    }
    write_table(outdir / "selected_metrics.tsv", selected_rows, list(selected_rows[0]))
    result = {
        "complete": True,
        "development_only": True,
        "development_seeds": development,
        "reserved_heldout_seeds": heldout,
        "candidate_models": len(candidates),
        "candidate_family_rows": len(rows),
        "eligible_candidates": sum(row["eligible"] for row in summaries),
        "selected_candidate": selected_id,
        "baseline_confusion": baseline["baseline_full"],
        "baseline_confusion_by_seed": baseline["baseline_by_seed"],
        "selected_confusion": confusion(selected_rows),
        "selected_confusion_by_seed": selected_by_seed,
        "acceptance": acceptance,
        "scientific_acceptance": "development_refinement_reported_not_independent_validation",
    }
    (outdir / "validation.json").write_text(json.dumps(result, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run(args.previous, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
