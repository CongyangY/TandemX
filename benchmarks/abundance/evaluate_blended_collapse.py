"""Select a transparent single/multi-k assembly classifier on development seeds."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import shutil

from benchmarks.abundance.evaluate_multik_collapse import score_estimate
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


BASELINE = "single_k21"
MULTIK = "multik_loglinear"
KEY_FIELDS = (
    "seed", "unit_substitution_rate", "array_fragments", "coverage",
    "substitution_rate", "assembly_fraction", "family_id",
)


def confusion(rows: list[dict]) -> dict[str, float | int]:
    counts = Counter(row["outcome"] for row in rows)
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    available = tp + fn + fp + tn
    return {
        "available": available,
        "unavailable": counts["NA"],
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else math.nan,
        "false_positive_rate": fp / (fp + tn) if fp + tn else math.nan,
        "precision": tp / (tp + fp) if tp + fp else math.nan,
    }


def _estimate(single: dict, multi: dict, alpha: float) -> tuple[float, str]:
    single_estimate = float(single["read_estimated_copies"])
    if alpha == 0:
        return single_estimate, "single_k21_anchor"
    if multi["outcome"] == "NA" or multi["read_estimated_copies"] == "NA":
        return single_estimate, "fallback_single_k21:multik_unavailable"
    multi_estimate = float(multi["read_estimated_copies"])
    if single_estimate <= 0 or multi_estimate <= 0:
        return single_estimate, "fallback_single_k21:nonpositive_blend_input"
    estimate = math.exp((1 - alpha) * math.log(single_estimate) + alpha * math.log(multi_estimate))
    return estimate, f"log_space_blend:alpha={alpha:g}"


def blend_row(pair: dict[str, dict], alpha: float, decision_threshold: float) -> dict:
    single, multi = pair[BASELINE], pair[MULTIK]
    estimate, fit_status = _estimate(single, multi, alpha)
    candidate_id = f"blend_a{alpha:g}_t{decision_threshold:g}"
    scored = score_estimate(
        family_id=single["family_id"],
        period=int(single["period"]),
        estimate=estimate,
        assembly_bp=float(single["assembly_predicted_bp"]),
        truth_ratio=float(single["truth_assembly_read_ratio"]),
        threshold=float(single["decision_threshold"]),
        decision_threshold=decision_threshold,
        method=candidate_id,
        fit_status=fit_status,
    )
    scored.pop("family_id")
    context = {
        key: value for key, value in single.items()
        if key not in {
            "method", "read_estimated_copies", "read_estimated_bp",
            "predicted_assembly_read_ratio", "native_status",
            "predicted_underrepresented", "outcome", "fit_status",
            "decision_threshold",
        }
    }
    return {
        **context,
        "candidate_id": candidate_id,
        "blend_alpha": alpha,
        **scored,
    }


def select_candidate(
    candidates: dict[str, list[dict]],
    baseline_rows: list[dict],
    rule: dict,
) -> tuple[str, dict[str, float | int]]:
    baseline = confusion(baseline_rows)
    constraints = rule["selection_constraints"]
    eligible = []
    for candidate_id, rows in candidates.items():
        metrics = confusion(rows)
        if (
            metrics["unavailable"] == 0
            and metrics["false_positive_rate"] <= baseline["false_positive_rate"]
                + float(constraints["false_positive_rate_max_delta_vs_single_k21"]) + 1e-15
            and metrics["precision"] >= baseline["precision"]
                + float(constraints["precision_min_delta_vs_single_k21"]) - 1e-15
        ):
            first = rows[0]
            eligible.append((
                metrics["sensitivity"],
                -metrics["false_positive_rate"],
                metrics["precision"],
                -float(first["blend_alpha"]),
                -float(first["decision_threshold"]),
                candidate_id,
                metrics,
            ))
    if not eligible:
        raise ValueError("No candidate satisfies the predeclared baseline constraints")
    selected = max(eligible)
    return selected[-2], selected[-1]


def _metric_fields(prefix: str, values: dict[str, float | int]) -> dict:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def _validate_rule(config: dict) -> dict:
    rule = config.get("classifier_development_rule")
    if not isinstance(rule, dict) or rule.get("method") != "log_space_single_multik_blend":
        raise ValueError("Missing predeclared classifier development rule")
    alphas = rule.get("alpha_grid")
    thresholds = rule.get("decision_threshold_grid")
    if (
        not isinstance(alphas, list) or not alphas or len(set(alphas)) != len(alphas)
        or any(not isinstance(value, (int, float)) or isinstance(value, bool)
               or not math.isfinite(value) or not 0 <= value <= 1 for value in alphas)
        or not isinstance(thresholds, list) or not thresholds
        or len(set(thresholds)) != len(thresholds)
        or any(not isinstance(value, (int, float)) or isinstance(value, bool)
               or not math.isfinite(value) or not 0 < value < 1 for value in thresholds)
        or rule.get("multik_k_values") != [15, 21, 27, 31]
        or rule.get("unavailable_or_nonpositive_rule") != "fallback_single_k21"
        or rule.get("baseline_anchor") != {"alpha": 0, "decision_threshold": 0.6}
        or rule.get("selection_objective")
            != "maximize_sensitivity_then_minimize_fpr_then_maximize_precision_then_minimize_alpha_then_threshold"
    ):
        raise ValueError("Invalid predeclared classifier candidate grid")
    constraints = rule.get("selection_constraints")
    gates = rule.get("development_acceptance_gates")
    required_constraints = {
        "false_positive_rate_max_delta_vs_single_k21",
        "precision_min_delta_vs_single_k21",
    }
    required_gates = {
        "full_sensitivity_min_delta_vs_single_k21",
        "full_false_positive_rate_max_delta_vs_single_k21",
        "full_precision_min_delta_vs_single_k21",
        "leave_one_seed_out_selection_consistent",
        "cross_validated_sensitivity_min_delta_vs_single_k21",
        "cross_validated_false_positive_rate_max_delta_vs_single_k21",
        "cross_validated_precision_min_delta_vs_single_k21",
    }
    if (
        not isinstance(constraints, dict) or set(constraints) != required_constraints
        or not isinstance(gates, dict) or set(gates) != required_gates
        or gates["leave_one_seed_out_selection_consistent"] is not True
    ):
        raise ValueError("Invalid classifier selection constraints or acceptance gates")
    return rule


def run(previous: Path, outdir: Path) -> None:
    previous = previous.resolve()
    validation = json.loads((previous / "validation.json").read_text())
    environment = json.loads((previous / "environment.json").read_text())
    config = json.loads((previous / "run_config.json").read_text())
    rule = _validate_rule(config)
    development_seeds = [int(seed) for seed in config.get("seeds", {}).get("development", [])]
    heldout_seeds = [int(seed) for seed in config.get("seeds", {}).get("heldout", [])]
    if (
        validation.get("complete") is not True
        or validation.get("split") != "development"
        or validation.get("evaluation_mode") != "raw_single_multik_for_predeclared_blend_grid"
        or environment.get("split") != "development"
        or digest_file(previous / "run_config.json") != environment.get("previous_config_sha256")
        or digest_file(previous / "comparison_metrics.tsv")
            != validation.get("comparison_metrics_sha256")
        or set(validation.get("method_confusion", {})) != {BASELINE, MULTIK}
        or len(development_seeds) != 3
        or len(set(development_seeds)) != 3
        or not heldout_seeds
        or set(development_seeds) & set(heldout_seeds)
    ):
        raise ValueError("Require complete three-seed development multi-k evidence")
    rows = read_table(previous / "comparison_metrics.tsv", {
        *KEY_FIELDS, "period", "truth_assembly_read_ratio", "assembly_predicted_bp",
        "read_estimated_copies", "method", "outcome", "decision_threshold",
    })
    paired: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        if row["method"] in {BASELINE, MULTIK}:
            paired[tuple(row[field] for field in KEY_FIELDS)][row["method"]] = row
    if not paired or any(set(pair) != {BASELINE, MULTIK} for pair in paired.values()):
        raise ValueError("Incomplete paired single-k and multi-k development rows")
    pairs = list(paired.values())
    if {int(pair[BASELINE]["seed"]) for pair in pairs} != set(development_seeds):
        raise ValueError("Development rows use unexpected seeds")

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
        split="development",
        development_seeds=development_seeds,
        reserved_heldout_seeds=heldout_seeds,
        scope="predeclared_transparent_classifier_development_not_heldout_validation",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(previous / "run_config.json", outdir / "run_config.json")

    candidate_rows = []
    candidates: dict[str, list[dict]] = defaultdict(list)
    for alpha in rule["alpha_grid"]:
        for threshold in rule["decision_threshold_grid"]:
            for pair in pairs:
                row = blend_row(pair, float(alpha), float(threshold))
                candidate_rows.append(row)
                candidates[row["candidate_id"]].append(row)
    baseline_rows = [pair[BASELINE] for pair in pairs]
    anchor = candidates["blend_a0_t0.6"]
    if [row["outcome"] for row in anchor] != [row["outcome"] for row in baseline_rows]:
        raise ValueError("Declared alpha-zero threshold-0.6 anchor differs from single-k baseline")

    selected_id, _ = select_candidate(candidates, baseline_rows, rule)
    folds = [("full_development", set(development_seeds), set(development_seeds))]
    folds.extend(
        (f"leave_s{heldout}_out", set(development_seeds) - {heldout}, {heldout})
        for heldout in development_seeds
    )
    selection_rows = []
    cross_validated_rows = []
    for fold, training, evaluated in folds:
        train_baseline = [row for row in baseline_rows if int(row["seed"]) in training]
        train_candidates = {
            candidate_id: [row for row in group if int(row["seed"]) in training]
            for candidate_id, group in candidates.items()
        }
        fold_id, train_selected = select_candidate(train_candidates, train_baseline, rule)
        eval_baseline_rows = [row for row in baseline_rows if int(row["seed"]) in evaluated]
        eval_candidate_rows = [row for row in candidates[fold_id] if int(row["seed"]) in evaluated]
        eval_baseline = confusion(eval_baseline_rows)
        eval_selected = confusion(eval_candidate_rows)
        first = train_candidates[fold_id][0]
        selection_rows.append({
            "fold": fold,
            "training_seeds": ",".join(map(str, sorted(training))),
            "evaluation_seeds": ",".join(map(str, sorted(evaluated))),
            "selected_candidate": fold_id,
            "blend_alpha": first["blend_alpha"],
            "decision_threshold": first["decision_threshold"],
            **_metric_fields("training_baseline", confusion(train_baseline)),
            **_metric_fields("training_selected", train_selected),
            **_metric_fields("evaluation_baseline", eval_baseline),
            **_metric_fields("evaluation_selected", eval_selected),
        })
        if fold != "full_development":
            cross_validated_rows.extend({"fold": fold, **row} for row in eval_candidate_rows)
    if selection_rows[0]["selected_candidate"] != selected_id:
        raise ValueError("Full-development selection is inconsistent")

    summaries = []
    for candidate_id, group in sorted(candidates.items()):
        first = group[0]
        summaries.append({
            "candidate_id": candidate_id,
            "blend_alpha": first["blend_alpha"],
            "decision_threshold": first["decision_threshold"],
            **confusion(group),
        })
    full_baseline = confusion(baseline_rows)
    full_selected = confusion(candidates[selected_id])
    cv_baseline = confusion(baseline_rows)
    cv_selected = confusion(cross_validated_rows)
    gates = rule["development_acceptance_gates"]
    stable = all(row["selected_candidate"] == selected_id for row in selection_rows[1:])
    acceptance = {
        "full_sensitivity_delta": full_selected["sensitivity"] - full_baseline["sensitivity"],
        "full_false_positive_rate_delta": (
            full_selected["false_positive_rate"] - full_baseline["false_positive_rate"]
        ),
        "full_precision_delta": full_selected["precision"] - full_baseline["precision"],
        "leave_one_seed_out_selection_consistent": stable,
        "cross_validated_sensitivity_delta": cv_selected["sensitivity"] - cv_baseline["sensitivity"],
        "cross_validated_false_positive_rate_delta": (
            cv_selected["false_positive_rate"] - cv_baseline["false_positive_rate"]
        ),
        "cross_validated_precision_delta": cv_selected["precision"] - cv_baseline["precision"],
    }
    acceptance["passed"] = (
        acceptance["full_sensitivity_delta"]
            >= float(gates["full_sensitivity_min_delta_vs_single_k21"])
        and acceptance["full_false_positive_rate_delta"]
            <= float(gates["full_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["full_precision_delta"]
            >= float(gates["full_precision_min_delta_vs_single_k21"])
        and stable is gates["leave_one_seed_out_selection_consistent"]
        and acceptance["cross_validated_sensitivity_delta"]
            >= float(gates["cross_validated_sensitivity_min_delta_vs_single_k21"])
        and acceptance["cross_validated_false_positive_rate_delta"]
            <= float(gates["cross_validated_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["cross_validated_precision_delta"]
            >= float(gates["cross_validated_precision_min_delta_vs_single_k21"])
    )
    write_table(outdir / "candidate_metrics.tsv", candidate_rows, list(candidate_rows[0]))
    write_table(outdir / "candidate_summary.tsv", summaries, list(summaries[0]))
    write_table(outdir / "selection.tsv", selection_rows, list(selection_rows[0]))
    write_table(outdir / "selected_metrics.tsv", candidates[selected_id], list(candidates[selected_id][0]))
    write_table(outdir / "cross_validated_metrics.tsv", cross_validated_rows,
                list(cross_validated_rows[0]))
    result = {
        "complete": True,
        "development_only": True,
        "independent_development_genomes": len(development_seeds),
        "reserved_heldout_genomes": len(heldout_seeds),
        "paired_family_conditions": len(pairs),
        "candidate_models": len(candidates),
        "candidate_family_rows": len(candidate_rows),
        "selected_candidate": selected_id,
        "baseline_confusion": full_baseline,
        "selected_confusion": full_selected,
        "cross_validated_confusion": cv_selected,
        "acceptance": acceptance,
        "scientific_acceptance": "development_gate_reported_not_assumed_from_execution_success",
    }
    (outdir / "validation.json").write_text(json.dumps(result, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run(args.previous, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
