"""Run independent conditional copy-number, localization and collapse experiments."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import sys

from benchmarks.abundance.evaluate import score_comparison, score_copy_number, score_localization
from benchmarks.abundance.simulate import (
    GenomeSpec,
    build_genome,
    challenge_scenarios,
    sample_reads,
    scenario_directory,
    write_genome,
)
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table


FROZEN_LOCALIZER_METHOD = "iid_base_exact_anchor_monomer_bridge"
FROZEN_LOCALIZER_FILES = {
    "development_validation_sha256": "validation.json",
    "development_metrics_sha256": "localization_metrics.tsv",
    "development_environment_sha256": "environment.json",
    "development_config_sha256": "run_config.json",
}
FROZEN_CLASSIFIER_METHOD = "log_space_single_multik_blend"
FROZEN_CLASSIFIER_FILES = {
    "development_validation_sha256": "validation.json",
    "development_selection_sha256": "selection.tsv",
    "development_candidate_summary_sha256": "candidate_robust_summary.tsv",
    "development_selected_metrics_sha256": "selected_metrics.tsv",
    "development_environment_sha256": "environment.json",
    "development_config_sha256": "run_config.json",
}
FROZEN_DEPTH_GATED_CLASSIFIER_METHOD = "depth_gated_log_space_blend"
FROZEN_DEPTH_GATED_CLASSIFIER_FILES = {
    "development_validation_sha256": "validation.json",
    "development_selection_sha256": "selection.tsv",
    "development_metrics_sha256": "comparison_metrics.tsv",
    "development_environment_sha256": "environment.json",
    "development_config_sha256": "run_config.json",
}


def aggregate(rows: list[dict], group_fields: list[str], kind: str) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in group_fields)].append(row)
    summaries = []
    for keys, group in sorted(grouped.items()):
        summary = dict(zip(group_fields, keys))
        summary["family_observations"] = len(group)
        if kind == "copy":
            summary.update(mean_signed_relative_error=statistics.mean(r["signed_relative_error"] for r in group),
                           median_absolute_relative_error=statistics.median(r["absolute_relative_error"] for r in group),
                           mean_absolute_relative_error=statistics.mean(r["absolute_relative_error"] for r in group),
                           interval_empirical_coverage=statistics.mean(r["interval_contains_truth"] for r in group),
                           mean_interval_relative_width=statistics.mean(r["interval_relative_width"] for r in group),
                           mean_oracle_relative_error=statistics.mean(r["oracle_relative_error"] for r in group),
                           mean_estimator_minus_oracle=statistics.mean(r["estimator_minus_sampling_oracle"] for r in group))
        else:
            counts = Counter(r["outcome"] for r in group)
            summary.update({key: counts[key] for key in ("TP", "FN", "FP", "TN")})
            for label, numerator, denominator in (
                ("sensitivity", counts["TP"], counts["TP"]+counts["FN"]),
                ("false_positive_rate", counts["FP"], counts["FP"]+counts["TN"]),
                ("precision", counts["TP"], counts["TP"]+counts["FP"]),
            ):
                summary[label] = numerator/denominator if denominator else None
        summaries.append(summary)
    return summaries


def validate_frozen_localizer(config: dict, split: str) -> dict[str, float] | None:
    """Verify the committed development evidence before an IID held-out run."""
    identity_model = config.get("locate_identity_model", "exact_kmer_fraction")
    if split != "heldout" or identity_model != "iid_base":
        return None
    model = config.get("localizer_model")
    heldout_seeds = config.get("seeds", {}).get("heldout", [])
    if (
        not isinstance(model, dict)
        or model.get("method") != FROZEN_LOCALIZER_METHOD
        or model.get("identity_model") != identity_model
        or model.get("k") != config.get("k")
        or model.get("min_identity") != config.get("locate_min_identity")
        or model.get("array_merge_gap_rule") != "max(2*k,monomer_length)"
        or not isinstance(model.get("development_seeds"), list)
        or not model["development_seeds"]
        or set(model["development_seeds"]) & set(heldout_seeds)
    ):
        raise ValueError("IID held-out validation requires a disjoint frozen localizer model")
    if any(
        not isinstance(model.get(field), str)
        or len(model[field]) != 64
        or any(character not in "0123456789abcdef" for character in model[field])
        for field in FROZEN_LOCALIZER_FILES
    ):
        raise ValueError("Frozen localizer evidence hashes must be lowercase SHA-256 values")

    development = Path(str(model.get("development_result", ""))).expanduser()
    if not development.is_dir():
        raise ValueError("Frozen localizer development result is unavailable")
    for field, filename in FROZEN_LOCALIZER_FILES.items():
        path = development / filename
        if not path.is_file() or digest_file(path) != model[field]:
            raise ValueError(f"Frozen localizer evidence differs: {filename}")

    validation = json.loads((development / "validation.json").read_text())
    environment = json.loads((development / "environment.json").read_text())
    development_config = json.loads((development / "run_config.json").read_text())
    if (
        validation.get("complete") is not True
        or validation.get("mode") != "localization_only"
        or validation.get("successful") != validation.get("executions")
        or validation.get("localization_family_rows") <= 0
        or environment.get("split") != "development"
        or environment.get("localization_only") is not True
        or environment.get("locate_identity_model") != identity_model
        or environment.get("locate_min_identity") != config.get("locate_min_identity")
        or environment.get("config_sha256") != model["development_config_sha256"]
        or development_config.get("seeds", {}).get("development") != model["development_seeds"]
        or development_config.get("locate_identity_model") != identity_model
        or development_config.get("locate_min_identity") != config.get("locate_min_identity")
    ):
        raise ValueError("Frozen localizer development provenance is inconsistent")

    rows = read_table(development / "localization_metrics.tsv", {
        "assembly_fraction", "base_recall", "base_precision", "predicted_assembly_bp"
    })
    if len(rows) != validation["localization_family_rows"]:
        raise ValueError("Frozen localizer metric row count differs from validation")
    full = [row for row in rows if float(row["assembly_fraction"]) == 1.0]
    positive = [row for row in rows if float(row["assembly_fraction"]) > 0.0]
    absent = [row for row in rows if float(row["assembly_fraction"]) == 0.0]
    if not full or not positive or not absent:
        raise ValueError("Frozen localizer evidence lacks a required assembly stratum")
    observed = {
        "full_assembly_mean_base_recall": statistics.mean(float(row["base_recall"]) for row in full),
        "positive_assembly_mean_base_precision": statistics.mean(float(row["base_precision"]) for row in positive),
        "absent_family_false_positive_rate": (
            sum(float(row["predicted_assembly_bp"]) > 0 for row in absent) / len(absent)
        ),
    }
    declared = model.get("observed_development_metrics")
    gates = model.get("selection_gates")
    if (
        not isinstance(declared, dict)
        or any(
            key not in declared
            or not math.isclose(float(declared[key]), value, rel_tol=0, abs_tol=1e-12)
            for key, value in observed.items()
        )
        or not isinstance(gates, dict)
        or observed["full_assembly_mean_base_recall"]
            < float(gates.get("full_assembly_mean_base_recall_min", math.inf))
        or observed["positive_assembly_mean_base_precision"]
            < float(gates.get("positive_assembly_mean_base_precision_min", math.inf))
        or observed["absent_family_false_positive_rate"]
            > float(gates.get("absent_family_false_positive_rate_max", -math.inf))
    ):
        raise ValueError("Frozen localizer development metrics do not satisfy the declared gates")
    return observed


def _validate_frozen_depth_gated_classifier(config: dict, model: dict) -> dict[str, float]:
    """Verify post-failure depth-gated development before a new held-out run."""
    raw_rule = config.get("classifier_development_rule")
    heldout_seeds = config.get("seeds", {}).get("heldout", [])
    development_seeds = model.get("development_seeds")
    if (
        not isinstance(raw_rule, dict)
        or raw_rule.get("method") != FROZEN_CLASSIFIER_METHOD
        or raw_rule.get("alpha_grid") != [0, 0.25, 0.5, 0.75, 1]
        or raw_rule.get("decision_threshold_grid") != [0.45, 0.5, 0.55, 0.6]
        or raw_rule.get("multik_k_values") != [15, 21, 27, 31]
        or raw_rule.get("unavailable_or_nonpositive_rule") != "fallback_single_k21"
        or model.get("method") != FROZEN_DEPTH_GATED_CLASSIFIER_METHOD
        or model.get("selection_method") != "post_failed_heldout_coverage_diagnostic"
        or model.get("estimated_depth_source") != "single_k21_estimated_haploid_depth"
        or model.get("low_depth_cutoff") != 2.0
        or model.get("low_depth_strategy") != {
            "method": "single_k21", "blend_alpha": 0.0, "decision_threshold": 0.6
        }
        or model.get("standard_depth_strategy") != {
            "method": "log_space_single_multik_blend",
            "blend_alpha": 0.5,
            "decision_threshold": 0.5,
        }
        or model.get("k_values") != [15, 21, 27, 31]
        or model.get("fallback_rule") != "fallback_single_k21"
        or not isinstance(development_seeds, list) or len(development_seeds) != 6
        or len(set(development_seeds)) != 6
        or not isinstance(heldout_seeds, list) or len(heldout_seeds) != 3
        or len(set(heldout_seeds)) != 3
        or set(development_seeds) & set(heldout_seeds)
    ):
        raise ValueError("Held-out comparison requires a disjoint frozen depth-gated model")
    if any(
        not isinstance(model.get(field), str)
        or len(model[field]) != 64
        or any(character not in "0123456789abcdef" for character in model[field])
        for field in FROZEN_DEPTH_GATED_CLASSIFIER_FILES
    ):
        raise ValueError("Frozen depth-gated evidence hashes must be lowercase SHA-256 values")

    development = Path(str(model.get("development_result", ""))).expanduser()
    if not development.is_dir():
        raise ValueError("Frozen depth-gated development result is unavailable")
    for field, filename in FROZEN_DEPTH_GATED_CLASSIFIER_FILES.items():
        path = development / filename
        if not path.is_file() or digest_file(path) != model[field]:
            raise ValueError(f"Frozen depth-gated evidence differs: {filename}")

    validation = json.loads((development / "validation.json").read_text())
    environment = json.loads((development / "environment.json").read_text())
    development_config = json.loads((development / "run_config.json").read_text())
    development_rule = development_config.get("classifier_rule")
    selection = read_table(development / "selection.tsv", {
        "candidate_id", "method", "low_depth_cutoff", "low_depth_method",
        "low_depth_blend_alpha", "low_depth_decision_threshold", "standard_method",
        "standard_blend_alpha", "standard_decision_threshold", "passed",
    })
    configured_development = config.get("seeds", {}).get("development")
    source_seeds = [
        seed for source in development_config.get("development_sources", [])
        for seed in source.get("seeds", [])
    ]
    if (
        validation.get("complete") is not True
        or validation.get("development_only") is not True
        or validation.get("post_failed_heldout_refinement") is not True
        or validation.get("selected_candidate") != "depth_gated_blend_v3"
        or validation.get("consumed_development_seeds") != development_seeds
        or validation.get("reserved_future_heldout_seeds") != heldout_seeds
        or validation.get("acceptance", {}).get("passed") is not True
        or validation.get("comparison_metrics_sha256")
            != model["development_metrics_sha256"]
        or validation.get("selection_sha256")
            != model["development_selection_sha256"]
        or environment.get("consumed_development_seeds") != development_seeds
        or environment.get("reserved_future_heldout_seeds") != heldout_seeds
        or environment.get("config_sha256") != model["development_config_sha256"]
        or development_config.get("reserved_future_heldout_seeds") != heldout_seeds
        or sorted(source_seeds) != sorted(development_seeds)
        or configured_development != development_seeds
        or not isinstance(development_rule, dict)
        or development_rule.get("method") != model["method"]
        or development_rule.get("estimated_depth_source")
            != model["estimated_depth_source"]
        or development_rule.get("low_depth_cutoff") != model["low_depth_cutoff"]
        or development_rule.get("low_depth_strategy") != model["low_depth_strategy"]
        or development_rule.get("standard_depth_strategy")
            != model["standard_depth_strategy"]
        or development_rule.get("multik_k_values") != model["k_values"]
        or development_rule.get("unavailable_or_nonpositive_rule")
            != model["fallback_rule"]
        or len(selection) != 1
        or selection[0]["candidate_id"] != "depth_gated_blend_v3"
        or selection[0]["method"] != model["method"]
        or not math.isclose(float(selection[0]["low_depth_cutoff"]), 2.0,
                            rel_tol=0, abs_tol=1e-12)
        or selection[0]["low_depth_method"] != "single_k21"
        or not math.isclose(float(selection[0]["low_depth_blend_alpha"]), 0.0,
                            rel_tol=0, abs_tol=1e-12)
        or not math.isclose(float(selection[0]["low_depth_decision_threshold"]), 0.6,
                            rel_tol=0, abs_tol=1e-12)
        or selection[0]["standard_method"] != "log_space_single_multik_blend"
        or not math.isclose(float(selection[0]["standard_blend_alpha"]), 0.5,
                            rel_tol=0, abs_tol=1e-12)
        or not math.isclose(float(selection[0]["standard_decision_threshold"]), 0.5,
                            rel_tol=0, abs_tol=1e-12)
        or selection[0]["passed"].lower() not in {"true", "1"}
    ):
        raise ValueError("Frozen depth-gated development provenance is inconsistent")

    metric_names = (
        "full_sensitivity_delta", "full_false_positive_rate_delta",
        "full_precision_delta", "minimum_cohort_sensitivity_delta",
        "maximum_cohort_false_positive_rate_delta", "minimum_cohort_precision_delta",
        "minimum_seed_sensitivity_delta", "maximum_seed_false_positive_rate_delta",
        "minimum_seed_precision_delta",
    )
    observed = {name: float(validation["acceptance"][name]) for name in metric_names}
    declared = model.get("observed_development_metrics")
    gates = model.get("selection_gates")
    if (
        not isinstance(declared, dict) or set(declared) != set(observed)
        or any(not math.isclose(float(declared[key]), value, rel_tol=0, abs_tol=1e-12)
               for key, value in observed.items())
        or not isinstance(gates, dict) or set(gates) != set(observed)
        or observed["full_sensitivity_delta"] < float(gates["full_sensitivity_delta"])
        or observed["full_false_positive_rate_delta"]
            > float(gates["full_false_positive_rate_delta"])
        or observed["full_precision_delta"] < float(gates["full_precision_delta"])
        or observed["minimum_cohort_sensitivity_delta"]
            < float(gates["minimum_cohort_sensitivity_delta"])
        or observed["maximum_cohort_false_positive_rate_delta"]
            > float(gates["maximum_cohort_false_positive_rate_delta"])
        or observed["minimum_cohort_precision_delta"]
            < float(gates["minimum_cohort_precision_delta"])
        or observed["minimum_seed_sensitivity_delta"]
            < float(gates["minimum_seed_sensitivity_delta"])
        or observed["maximum_seed_false_positive_rate_delta"]
            > float(gates["maximum_seed_false_positive_rate_delta"])
        or observed["minimum_seed_precision_delta"]
            < float(gates["minimum_seed_precision_delta"])
    ):
        raise ValueError("Frozen depth-gated metrics do not satisfy declared gates")
    return observed


def validate_frozen_classifier(config: dict, split: str) -> dict[str, float] | None:
    """Verify classifier development evidence before creating held-out output."""
    rule = config.get("classifier_development_rule")
    model = config.get("classifier_model")
    if split != "heldout" or rule is None and model is None:
        return None
    if isinstance(model, dict) and model.get("method") == FROZEN_DEPTH_GATED_CLASSIFIER_METHOD:
        return _validate_frozen_depth_gated_classifier(config, model)
    heldout_seeds = config.get("seeds", {}).get("heldout", [])
    if (
        not isinstance(rule, dict)
        or rule.get("method") != FROZEN_CLASSIFIER_METHOD
        or rule.get("alpha_grid") != [0, 0.25, 0.5, 0.75, 1]
        or rule.get("decision_threshold_grid") != [0.45, 0.5, 0.55, 0.6]
        or rule.get("multik_k_values") != [15, 21, 27, 31]
        or rule.get("unavailable_or_nonpositive_rule") != "fallback_single_k21"
        or not isinstance(model, dict)
        or model.get("method") != FROZEN_CLASSIFIER_METHOD
        or model.get("selection_method") != "seed_robust_minimax"
        or model.get("k_values") != [15, 21, 27, 31]
        or model.get("fallback_rule") != "fallback_single_k21"
        or model.get("blend_alpha") not in rule["alpha_grid"]
        or model.get("decision_threshold") not in rule["decision_threshold_grid"]
        or not isinstance(model.get("development_seeds"), list)
        or not model["development_seeds"]
        or set(model["development_seeds"]) & set(heldout_seeds)
    ):
        raise ValueError("Held-out comparison requires a disjoint frozen classifier model")
    if any(
        not isinstance(model.get(field), str)
        or len(model[field]) != 64
        or any(character not in "0123456789abcdef" for character in model[field])
        for field in FROZEN_CLASSIFIER_FILES
    ):
        raise ValueError("Frozen classifier evidence hashes must be lowercase SHA-256 values")

    development = Path(str(model.get("development_result", ""))).expanduser()
    if not development.is_dir():
        raise ValueError("Frozen classifier development result is unavailable")
    for field, filename in FROZEN_CLASSIFIER_FILES.items():
        path = development / filename
        if not path.is_file() or digest_file(path) != model[field]:
            raise ValueError(f"Frozen classifier evidence differs: {filename}")

    validation = json.loads((development / "validation.json").read_text())
    environment = json.loads((development / "environment.json").read_text())
    development_config = json.loads((development / "run_config.json").read_text())
    selection = read_table(development / "selection.tsv", {
        "candidate_id", "blend_alpha", "decision_threshold", "selection_method"
    })
    candidate_id = f"blend_a{float(model['blend_alpha']):g}_t{float(model['decision_threshold']):g}"
    if (
        validation.get("complete") is not True
        or validation.get("development_only") is not True
        or validation.get("selected_candidate") != candidate_id
        or validation.get("development_seeds") != model["development_seeds"]
        or validation.get("reserved_heldout_seeds") != heldout_seeds
        or validation.get("acceptance", {}).get("passed") is not True
        or environment.get("development_seeds") != model["development_seeds"]
        or environment.get("reserved_heldout_seeds") != heldout_seeds
        or environment.get("config_sha256") != model["development_config_sha256"]
        or development_config.get("development_seeds") != model["development_seeds"]
        or development_config.get("reserved_heldout_seeds") != heldout_seeds
        or len(selection) != 1
        or selection[0]["candidate_id"] != candidate_id
        or not math.isclose(float(selection[0]["blend_alpha"]),
                            float(model["blend_alpha"]), rel_tol=0, abs_tol=1e-12)
        or not math.isclose(float(selection[0]["decision_threshold"]),
                            float(model["decision_threshold"]), rel_tol=0, abs_tol=1e-12)
        or selection[0]["selection_method"] != model["selection_method"]
    ):
        raise ValueError("Frozen classifier development provenance is inconsistent")

    observed = {
        key: float(validation["acceptance"][key]) for key in (
            "minimum_seed_sensitivity_delta", "maximum_seed_false_positive_rate_delta",
            "minimum_seed_precision_delta", "full_sensitivity_delta",
            "full_false_positive_rate_delta", "full_precision_delta",
        )
    }
    declared = model.get("observed_development_metrics")
    gates = model.get("selection_gates")
    if (
        not isinstance(declared, dict)
        or set(declared) != set(observed)
        or any(not math.isclose(float(declared[key]), value, rel_tol=0, abs_tol=1e-12)
               for key, value in observed.items())
        or not isinstance(gates, dict)
        or observed["minimum_seed_sensitivity_delta"]
            < float(gates.get("minimum_seed_sensitivity_delta", math.inf))
        or observed["maximum_seed_false_positive_rate_delta"]
            > float(gates.get("maximum_seed_false_positive_rate_delta", -math.inf))
        or observed["minimum_seed_precision_delta"]
            < float(gates.get("minimum_seed_precision_delta", math.inf))
        or observed["full_sensitivity_delta"]
            < float(gates.get("full_sensitivity_delta", math.inf))
        or observed["full_false_positive_rate_delta"]
            > float(gates.get("full_false_positive_rate_delta", -math.inf))
        or observed["full_precision_delta"]
            < float(gates.get("full_precision_delta", math.inf))
    ):
        raise ValueError("Frozen classifier development metrics do not satisfy declared gates")
    return observed


def run(config_path: Path, outdir: Path, split: str, localization_only: bool = False) -> None:
    config = json.loads(config_path.read_text())
    seeds = [s for group in config["seeds"].values() for s in group]
    if len(set(seeds)) != len(seeds) or split not in config["seeds"]:
        raise ValueError("Seed groups must be disjoint and requested split must exist")
    for field in ("coverages", "substitution_rates", "assembly_fractions"):
        values = config[field]
        if not values or len(set(values)) != len(values) or any(not math.isfinite(x) for x in values):
            raise ValueError(f"{field} must contain unique finite values")
    if (not config['seeds'][split] or any(x<=0 for x in config['coverages'])
        or any(not 0<=x<1 for x in config['substitution_rates'])
        or any(not 0<=x<=2 for x in config['assembly_fractions'])
        or not 0<config['collapse_threshold']<1 or config['timeout_seconds']<=0
        or not 1<=config['k']<=31 or config['read_length']<1):
        raise ValueError("Invalid abundance experiment configuration")
    scenarios = challenge_scenarios(config)
    identity_model = config.get("locate_identity_model", "exact_kmer_fraction")
    locate_min_identity = config.get("locate_min_identity", 0.8)
    if (
        identity_model not in {"exact_kmer_fraction", "iid_base"}
        or not isinstance(locate_min_identity, (int, float))
        or isinstance(locate_min_identity, bool)
        or not math.isfinite(locate_min_identity)
        or not 0 < locate_min_identity <= 1
    ):
        raise ValueError("Invalid localization identity model or threshold")
    localizer_development_metrics = validate_frozen_localizer(config, split)
    classifier_development_metrics = validate_frozen_classifier(config, split)
    fragment_gap_bp = config.get("fragment_gap_bp", 0)
    if not isinstance(fragment_gap_bp, int) or isinstance(fragment_gap_bp, bool) or fragment_gap_bp < 0:
        raise ValueError("fragment_gap_bp must be a nonnegative integer")
    for unit_rate, fragment_count in scenarios:
        GenomeSpec(
            seed=seeds[0], periods=tuple(config['periods']), copies=tuple(config['copies']),
            flank_bp=config['flank_bp'], unit_substitution_rate=unit_rate,
            array_fragments=fragment_count, fragment_gap_bp=fragment_gap_bp,
        ).validate()
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    source = outdir / "source_snapshot"
    provenance = source_manifest(root, source)
    benchmark_hashes = {}
    for path in sorted(Path(__file__).parent.glob("*.py")):
        relative = path.relative_to(root)
        target = source / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        benchmark_hashes[str(relative)] = digest_file(target)
        if digest_file(path) != digest_file(target):
            raise ValueError("Benchmark source changed during snapshot")
    provenance.update(benchmark_source_sha256=benchmark_hashes, platform=platform.platform(), python=sys.version,
                      config_sha256=digest_file(config_path), split=split,
                      challenge_scenarios=[dict(unit_substitution_rate=rate, array_fragments=count,
                                                fragment_gap_bp=fragment_gap_bp)
                                           for rate, count in scenarios],
                      locate_identity_model=identity_model,
                      locate_min_identity=locate_min_identity,
                      frozen_localizer_development_metrics=localizer_development_metrics,
                      frozen_classifier_development_metrics=classifier_development_metrics,
                      localization_only=localization_only,
                      scope=("known-catalogue conditional quantification/localization with explicit "
                             "unit-divergence and array-fragmentation factors; not discovery or empirical HiFi validation"),
                      resource_note="sequential direct-child wait4; development diagnostics, not external superiority timing")
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2)+"\n")
    shutil.copyfile(config_path, outdir / "run_config.json")
    receipts, cn_scores, locate_scores, compare_scores = [], [], [], []

    def execute(label: str, directory: Path, arguments: list[str]) -> bool:
        directory.mkdir(parents=True)
        command = [sys.executable, "-m", "tandemx.cli", *arguments, "--outdir", str(directory / "output")]
        measured = run_process(command, directory/"stdout.log", directory/"stderr.log", config["timeout_seconds"],
                               {**os.environ, "PYTHONPATH": str(source)}, source)
        receipt = dict(label=label, command=command, **measured)
        receipts.append(receipt)
        (directory / "receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
        with (outdir / "run.log").open("a") as log:
            log.write(f"{label}\texit={measured['exit_code']}\t{directory}\n")
        return measured["exit_code"] == 0 and not measured["timed_out"]

    for seed in config["seeds"][split]:
        for scenario_index, (unit_rate, fragment_count) in enumerate(scenarios, 1):
            spec = GenomeSpec(
                seed=seed, periods=tuple(config["periods"]), copies=tuple(config["copies"]),
                flank_bp=config["flank_bp"], unit_substitution_rate=unit_rate,
                array_fragments=fragment_count, fragment_gap_bp=fragment_gap_bp,
            )
            genome_dir = scenario_directory(
                outdir / "genomes" / f"s{seed}", scenario_index, scenarios, config
            )
            run_dir = scenario_directory(
                outdir / "runs" / f"s{seed}", scenario_index, scenarios, config
            )
            reads_root = scenario_directory(
                outdir / "reads" / f"s{seed}", scenario_index, scenarios, config
            )
            manifest = write_genome(spec, genome_dir, tuple(config["assembly_fractions"]))
            genome, _, truth = build_genome(spec)
            catalogue = genome_dir / "catalogue.fa"
            scenario_context = dict(
                unit_substitution_rate=unit_rate, array_fragments=fragment_count,
                fragment_gap_bp=fragment_gap_bp, locate_identity_model=identity_model,
                locate_min_identity=locate_min_identity,
            )
            located = []
            for variant in manifest["variants"]:
                folder = run_dir / variant["name"] / "locate"
                ok = execute("locate", folder, ["locate", "--assembly", str(genome_dir/f"{variant['name']}.fa"),
                             "--catalog", str(catalogue), "--k", str(config["k"]),
                             "--identity-model", identity_model,
                             "--min-identity", str(locate_min_identity)])
                retained = read_table(genome_dir/f"{variant['name']}.truth.tsv")
                if ok:
                    for row in score_localization(folder/"output"/"arrays.bed", retained, variant["genome_bp"]):
                        locate_scores.append(dict(seed=seed, **scenario_context,
                                                  assembly_fraction=variant["fraction"], **row))
                located.append((variant, retained, folder/"output"/"arrays.bed", ok))
            if localization_only:
                continue
            for coverage in config["coverages"]:
                for error in config["substitution_rates"]:
                    name = f"c{coverage}_e{error}"
                    reads_dir = reads_root / name
                    sampling = sample_reads(genome, truth, reads_dir, seed=seed+1000003,
                                            coverage=coverage, read_length=config["read_length"], substitution_rate=error)
                    folder = run_dir / name / "quantify"
                    ok = execute("quantify", folder, ["quantify", "--reads", str(reads_dir/"reads.fa"), "--catalog", str(catalogue),
                                 "--genome-size", str(len(genome)), "--k", str(config["k"]), "--kmer-backend", "rust", "--no-progress"])
                    if not ok:
                        continue
                    context = dict(seed=seed, **scenario_context, coverage=coverage,
                                   substitution_rate=error, read_length=config["read_length"])
                    cn_scores.extend(dict(**context, **row) for row in score_copy_number(folder/"output"/"copy_number.tsv", truth, sampling))
                    for variant, retained, arrays, success in located:
                        if not success:
                            continue
                        comp = run_dir / name / variant["name"] / "compare"
                        if execute("compare", comp, ["compare", "--copy-number", str(folder/"output"/"copy_number.tsv"),
                                   "--arrays", str(arrays), "--collapse-threshold", str(config["collapse_threshold"])]):
                            compare_scores.extend(dict(**context, assembly_fraction=variant["fraction"], **row) for row in
                                                  score_comparison(comp/"output"/"assembly_vs_read_cn.tsv", truth, retained, config["collapse_threshold"]))
    for filename, rows in (("copy_number_metrics.tsv", cn_scores), ("localization_metrics.tsv", locate_scores),
                           ("comparison_metrics.tsv", compare_scores),
                           ("copy_number_summary.tsv", aggregate(cn_scores, ["unit_substitution_rate", "array_fragments", "coverage", "substitution_rate"], "copy")),
                           ("comparison_summary.tsv", aggregate(compare_scores, ["unit_substitution_rate", "array_fragments", "coverage", "substitution_rate", "assembly_fraction"], "comparison"))):
        if rows:
            write_table(outdir/filename, rows, list(rows[0]))
    validation = dict(complete=True, mode="localization_only" if localization_only else "full",
                      executions=len(receipts), successful=sum(r["exit_code"]==0 for r in receipts),
                      challenge_scenarios=len(scenarios),
                      copy_number_family_rows=len(cn_scores), localization_family_rows=len(locate_scores),
                      comparison_family_rows=len(compare_scores), scientific_acceptance="not_assumed_from_execution_success")
    (outdir / "validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    if validation["executions"] != validation["successful"]:
        raise RuntimeError("Failed stage(s); retain receipts and do not turn missing metrics into zeros")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--split", choices=("development", "heldout"), default="development")
    parser.add_argument("--localization-only", action="store_true",
                        help="Run only assembly localization and skip read simulation, quantification and comparison.")
    args = parser.parse_args()
    run(args.config, args.outdir, args.split, args.localization_only)
