"""Evaluate a transparent depth-gated blend on consumed development evidence."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil

from benchmarks.abundance.evaluate_blended_collapse import KEY_FIELDS, confusion
from benchmarks.abundance.evaluate_multik_collapse import summarize
from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


BASELINE = "single_k21"
HIGH_DEPTH_SOURCE = "blend_a0.5_t0.5"
MODEL = "depth_gated_blend_v3"
HASH_FILES = {
    "metrics_sha256": None,
    "validation_sha256": "validation.json",
    "environment_sha256": "environment.json",
    "run_config_sha256": "run_config.json",
}


def validate_config(config: dict) -> tuple[list[dict], list[int], dict]:
    sources = config.get("development_sources")
    future = config.get("reserved_future_heldout_seeds")
    rule = config.get("classifier_rule")
    if (
        not isinstance(sources, list) or len(sources) != 2
        or {source.get("name") for source in sources}
            != {"original_development", "consumed_heldout_postmortem"}
        or not isinstance(future, list) or len(future) != 3 or len(set(future)) != 3
        or not all(isinstance(seed, int) and not isinstance(seed, bool) for seed in future)
        or not isinstance(rule, dict)
        or rule.get("method") != "depth_gated_log_space_blend"
        or rule.get("estimated_depth_source") != "single_k21_estimated_haploid_depth"
        or rule.get("low_depth_cutoff") != 2.0
        or rule.get("low_depth_strategy") != {
            "method": "single_k21", "blend_alpha": 0.0, "decision_threshold": 0.6
        }
        or rule.get("standard_depth_strategy") != {
            "method": "log_space_single_multik_blend",
            "blend_alpha": 0.5,
            "decision_threshold": 0.5,
        }
        or rule.get("multik_k_values") != [15, 21, 27, 31]
        or rule.get("unavailable_or_nonpositive_rule") != "fallback_single_k21"
    ):
        raise ValueError("Invalid depth-gated development configuration")
    gates = rule.get("development_acceptance_gates")
    expected_gates = {
        "full_sensitivity_min_delta_vs_single_k21",
        "full_false_positive_rate_max_delta_vs_single_k21",
        "full_precision_min_delta_vs_single_k21",
        "minimum_cohort_sensitivity_min_delta_vs_single_k21",
        "maximum_cohort_false_positive_rate_max_delta_vs_single_k21",
        "minimum_cohort_precision_min_delta_vs_single_k21",
        "minimum_seed_sensitivity_min_delta_vs_single_k21",
        "maximum_seed_false_positive_rate_max_delta_vs_single_k21",
        "minimum_seed_precision_min_delta_vs_single_k21",
    }
    if not isinstance(gates, dict) or set(gates) != expected_gates:
        raise ValueError("Invalid depth-gated development acceptance gates")
    consumed = []
    for source in sources:
        seeds = source.get("seeds")
        if (
            not isinstance(seeds, list) or len(seeds) != 3 or len(set(seeds)) != 3
            or not all(isinstance(seed, int) and not isinstance(seed, bool) for seed in seeds)
            or not isinstance(source.get("result"), str)
            or not isinstance(source.get("metrics_file"), str)
            or set(source) != {
                "name", "role", "seeds", "result", "metrics_file", *HASH_FILES
            }
        ):
            raise ValueError("Invalid depth-gated development source")
        consumed.extend(seeds)
        for field in HASH_FILES:
            value = source.get(field)
            if (
                not isinstance(value, str) or len(value) != 64
                or set(value) - set("0123456789abcdef")
            ):
                raise ValueError("Development source hashes must be lowercase SHA-256")
    if len(set(consumed)) != 6 or set(consumed) & set(future):
        raise ValueError("Consumed and future seed groups must be disjoint")
    return sources, future, rule


def apply_depth_gate(single: dict, high_depth: dict, rule: dict) -> dict:
    if tuple(single[field] for field in KEY_FIELDS) != tuple(
        high_depth[field] for field in KEY_FIELDS
    ):
        raise ValueError("Depth-gated source rows are not paired")
    depth = float(single["estimated_haploid_depth"])
    if not math.isfinite(depth) or depth < 0:
        raise ValueError("Estimated haploid depth must be finite and nonnegative")
    low_depth = depth < float(rule["low_depth_cutoff"])
    source = single if low_depth else high_depth
    expected_threshold = (
        rule["low_depth_strategy"]["decision_threshold"]
        if low_depth else rule["standard_depth_strategy"]["decision_threshold"]
    )
    if not math.isclose(float(source["decision_threshold"]), expected_threshold,
                        rel_tol=0, abs_tol=1e-12):
        raise ValueError("Depth-gated source decision threshold differs")
    row = dict(source)
    row.update(
        method=MODEL,
        candidate_id=MODEL,
        fit_status=(
            "depth_lt_2:single_k21_threshold_0.6"
            if low_depth else "depth_ge_2:blend_alpha_0.5_threshold_0.5"
        ),
    )
    return row


def _paired(rows: list[dict], source_name: str) -> list[dict[str, dict]]:
    pairs: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        candidate = row.get("candidate_id", row.get("method"))
        if candidate in {"blend_a0_t0.6", BASELINE}:
            method = BASELINE
        elif candidate == HIGH_DEPTH_SOURCE:
            method = HIGH_DEPTH_SOURCE
        else:
            continue
        key = tuple(row[field] for field in KEY_FIELDS)
        if method in pairs[key]:
            raise ValueError(f"Duplicate classifier row in {source_name}")
        pairs[key][method] = row
    if not pairs or any(set(pair) != {BASELINE, HIGH_DEPTH_SOURCE} for pair in pairs.values()):
        raise ValueError(f"Incomplete classifier source pair in {source_name}")
    return list(pairs.values())


def _delta(selected: dict, baseline: dict, metric: str) -> float:
    return float(selected[metric]) - float(baseline[metric])


def evaluate_acceptance(
    baseline_rows: list[dict],
    selected_rows: list[dict],
    cohort_by_seed: dict[int, str],
    gates: dict,
) -> tuple[dict, dict, dict]:
    full_baseline = confusion(baseline_rows)
    full_selected = confusion(selected_rows)
    by_seed = {}
    for seed in sorted(cohort_by_seed):
        b = confusion([row for row in baseline_rows if int(row["seed"]) == seed])
        s = confusion([row for row in selected_rows if int(row["seed"]) == seed])
        by_seed[seed] = {
            "cohort": cohort_by_seed[seed], "baseline": b, "selected": s,
            "sensitivity_delta": _delta(s, b, "sensitivity"),
            "false_positive_rate_delta": _delta(s, b, "false_positive_rate"),
            "precision_delta": _delta(s, b, "precision"),
        }
    by_cohort = {}
    for cohort in sorted(set(cohort_by_seed.values())):
        seeds = {seed for seed, value in cohort_by_seed.items() if value == cohort}
        b = confusion([row for row in baseline_rows if int(row["seed"]) in seeds])
        s = confusion([row for row in selected_rows if int(row["seed"]) in seeds])
        by_cohort[cohort] = {
            "baseline": b, "selected": s,
            "sensitivity_delta": _delta(s, b, "sensitivity"),
            "false_positive_rate_delta": _delta(s, b, "false_positive_rate"),
            "precision_delta": _delta(s, b, "precision"),
        }
    acceptance = {
        "full_sensitivity_delta": _delta(full_selected, full_baseline, "sensitivity"),
        "full_false_positive_rate_delta": _delta(
            full_selected, full_baseline, "false_positive_rate"
        ),
        "full_precision_delta": _delta(full_selected, full_baseline, "precision"),
        "minimum_cohort_sensitivity_delta": min(
            row["sensitivity_delta"] for row in by_cohort.values()
        ),
        "maximum_cohort_false_positive_rate_delta": max(
            row["false_positive_rate_delta"] for row in by_cohort.values()
        ),
        "minimum_cohort_precision_delta": min(
            row["precision_delta"] for row in by_cohort.values()
        ),
        "minimum_seed_sensitivity_delta": min(
            row["sensitivity_delta"] for row in by_seed.values()
        ),
        "maximum_seed_false_positive_rate_delta": max(
            row["false_positive_rate_delta"] for row in by_seed.values()
        ),
        "minimum_seed_precision_delta": min(
            row["precision_delta"] for row in by_seed.values()
        ),
    }
    acceptance["passed"] = (
        acceptance["full_sensitivity_delta"]
            >= float(gates["full_sensitivity_min_delta_vs_single_k21"])
        and acceptance["full_false_positive_rate_delta"]
            <= float(gates["full_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["full_precision_delta"]
            >= float(gates["full_precision_min_delta_vs_single_k21"])
        and acceptance["minimum_cohort_sensitivity_delta"]
            >= float(gates["minimum_cohort_sensitivity_min_delta_vs_single_k21"])
        and acceptance["maximum_cohort_false_positive_rate_delta"]
            <= float(gates["maximum_cohort_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["minimum_cohort_precision_delta"]
            >= float(gates["minimum_cohort_precision_min_delta_vs_single_k21"])
        and acceptance["minimum_seed_sensitivity_delta"]
            >= float(gates["minimum_seed_sensitivity_min_delta_vs_single_k21"])
        and acceptance["maximum_seed_false_positive_rate_delta"]
            <= float(gates["maximum_seed_false_positive_rate_max_delta_vs_single_k21"])
        and acceptance["minimum_seed_precision_delta"]
            >= float(gates["minimum_seed_precision_min_delta_vs_single_k21"])
    )
    return acceptance, by_cohort, by_seed


def _summary_rows(grouped: dict, group_name: str) -> list[dict]:
    output = []
    for name, values in grouped.items():
        for method in ("baseline", "selected"):
            output.append({group_name: name, "method": BASELINE if method == "baseline" else MODEL,
                           **values[method]})
    return output


def run(config_path: Path, outdir: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    config = json.loads(config_path.read_text())
    sources, future, rule = validate_config(config)
    baseline_rows = []
    selected_rows = []
    cohort_by_seed = {}
    input_receipts = []
    field_set = None
    for source in sources:
        source_root = Path(source["result"])
        if not source_root.is_absolute():
            source_root = root / source_root
        for field, filename in HASH_FILES.items():
            actual = source["metrics_file"] if filename is None else filename
            if digest_file(source_root / actual) != source[field]:
                raise ValueError(f"Frozen postmortem source differs: {source['name']}/{actual}")
        validation = json.loads((source_root / "validation.json").read_text())
        if validation.get("complete") is not True:
            raise ValueError("Depth-gated development source is incomplete")
        rows = read_table(source_root / source["metrics_file"], {
            *KEY_FIELDS, "estimated_haploid_depth", "candidate_id", "method", "outcome",
            "decision_threshold",
        })
        pairs = _paired(rows, source["name"])
        observed_seeds = {int(pair[BASELINE]["seed"]) for pair in pairs}
        if observed_seeds != set(source["seeds"]):
            raise ValueError("Depth-gated development source seeds differ")
        for pair in pairs:
            single = dict(pair[BASELINE])
            single.update(method=BASELINE, candidate_id=BASELINE, blend_alpha=0.0)
            chosen = apply_depth_gate(single, pair[HIGH_DEPTH_SOURCE], rule)
            current_fields = set(single)
            if field_set is None:
                field_set = current_fields
            if current_fields != field_set or set(chosen) != field_set:
                raise ValueError("Depth-gated cohorts have different metric schemas")
            baseline_rows.append(single)
            selected_rows.append(chosen)
        cohort_by_seed.update({seed: source["name"] for seed in source["seeds"]})
        input_receipts.append({
            "name": source["name"], "role": source["role"],
            "seeds": source["seeds"], "result": str(source_root.resolve()),
            **{field: source[field] for field in HASH_FILES},
        })
    if set(cohort_by_seed) & set(future):
        raise ValueError("Future held-out seeds overlap consumed development")
    acceptance, by_cohort, by_seed = evaluate_acceptance(
        baseline_rows, selected_rows, cohort_by_seed,
        rule["development_acceptance_gates"],
    )
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    snapshot = outdir / "source_snapshot"
    provenance = source_manifest(root, snapshot)
    relative = Path(__file__).resolve().relative_to(root)
    target = snapshot / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(__file__), target)
    provenance.update(
        benchmark_source_sha256={str(relative): digest_file(target)},
        config_sha256=digest_file(config_path),
        development_sources=input_receipts,
        consumed_development_seeds=sorted(cohort_by_seed),
        reserved_future_heldout_seeds=future,
        scope="post_failed_heldout_depth_gated_development_not_independent_validation",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(config_path, outdir / "run_config.json")
    combined = baseline_rows + selected_rows
    write_table(outdir / "comparison_metrics.tsv", combined, list(combined[0]))
    write_table(outdir / "comparison_summary.tsv", summarize(combined),
                list(summarize(combined)[0]))
    cohort_rows = _summary_rows(by_cohort, "cohort")
    seed_rows = _summary_rows(by_seed, "seed")
    write_table(outdir / "per_cohort_metrics.tsv", cohort_rows, list(cohort_rows[0]))
    write_table(outdir / "per_seed_metrics.tsv", seed_rows, list(seed_rows[0]))
    selection = [{
        "candidate_id": MODEL,
        "method": rule["method"],
        "low_depth_cutoff": rule["low_depth_cutoff"],
        "low_depth_method": rule["low_depth_strategy"]["method"],
        "low_depth_blend_alpha": rule["low_depth_strategy"]["blend_alpha"],
        "low_depth_decision_threshold": rule["low_depth_strategy"]["decision_threshold"],
        "standard_method": rule["standard_depth_strategy"]["method"],
        "standard_blend_alpha": rule["standard_depth_strategy"]["blend_alpha"],
        "standard_decision_threshold": rule["standard_depth_strategy"]["decision_threshold"],
        **acceptance,
    }]
    write_table(outdir / "selection.tsv", selection, list(selection[0]))
    result = {
        "complete": True,
        "development_only": True,
        "post_failed_heldout_refinement": True,
        "consumed_development_seeds": sorted(cohort_by_seed),
        "reserved_future_heldout_seeds": future,
        "paired_family_conditions": len(baseline_rows),
        "selected_candidate": MODEL,
        "low_depth_rows": sum(
            float(row["estimated_haploid_depth"]) < rule["low_depth_cutoff"]
            for row in baseline_rows
        ),
        "baseline_confusion": confusion(baseline_rows),
        "selected_confusion": confusion(selected_rows),
        "per_cohort": by_cohort,
        "per_seed": by_seed,
        "acceptance": acceptance,
        "comparison_metrics_sha256": digest_file(outdir / "comparison_metrics.tsv"),
        "selection_sha256": digest_file(outdir / "selection.tsv"),
        "scientific_acceptance": "development_gate_only_requires_fresh_5801_5803_validation",
    }
    (outdir / "validation.json").write_text(json.dumps(result, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run(args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
