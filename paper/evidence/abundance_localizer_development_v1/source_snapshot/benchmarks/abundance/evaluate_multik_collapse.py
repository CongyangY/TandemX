"""Calibrate or validate assembly comparison using single-k and multi-k estimates."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import re
import shutil

from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table, write_table
from benchmarks.abundance.simulate import challenge_scenarios, scenario_directory
from tandemx.compare.mvp import classify_assembly_read_ratio
from tandemx.discover.mvp import read_fasta
from tandemx.quantify.multik import DEFAULT_K_VALUES, estimate_multik
from tandemx.quantify.mvp import read_monomer_fasta


METHODS = ("single_k21", "multik_loglinear", "multik_fallback_depth_rule")
LOW_DEPTH_CUTOFF = 2.0
LOW_DEPTH_THRESHOLD_GRID = (0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6)


def score_estimate(
    *,
    family_id: str,
    period: int,
    estimate: float | None,
    assembly_bp: float,
    truth_ratio: float,
    threshold: float,
    decision_threshold: float | None = None,
    method: str,
    fit_status: str,
) -> dict:
    if period <= 0 or not math.isfinite(assembly_bp) or assembly_bp < 0:
        raise ValueError("Require a positive period and finite nonnegative assembly extent")
    decision = threshold if decision_threshold is None else decision_threshold
    if not 0 < threshold < 1 or not 0 < decision < 1:
        raise ValueError("Truth and decision thresholds must be in (0,1)")
    truth_positive = truth_ratio < threshold
    if estimate is None:
        return dict(
            family_id=family_id,
            method=method,
            read_estimated_copies="NA",
            read_estimated_bp="NA",
            predicted_assembly_read_ratio="NA",
            native_status="unavailable",
            predicted_underrepresented="NA",
            outcome="NA",
            fit_status=fit_status,
            decision_threshold=decision,
        )
    if not math.isfinite(estimate) or estimate < 0:
        raise ValueError("Read copy estimate must be finite and nonnegative")
    read_bp = estimate * period
    status, _confidence, _warning = classify_assembly_read_ratio(
        read_bp, assembly_bp, collapse_threshold=decision
    )
    call = status in {"possible_collapse", "reads_only"}
    return dict(
        family_id=family_id,
        method=method,
        read_estimated_copies=estimate,
        read_estimated_bp=read_bp,
        predicted_assembly_read_ratio=assembly_bp / read_bp if read_bp else "NA",
        native_status=status,
        predicted_underrepresented=call,
        outcome="TP" if truth_positive and call else "FN" if truth_positive else "FP" if call else "TN",
        fit_status=fit_status,
        decision_threshold=decision,
    )


def summarize(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["method"], row.get("unit_substitution_rate", 0.0),
                 row.get("array_fragments", 1),
                 row["coverage"], row["substitution_rate"], row["assembly_fraction"])].append(row)
    result = []
    for keys, group in sorted(grouped.items()):
        counts = Counter(row["outcome"] for row in group)
        available = sum(row["outcome"] != "NA" for row in group)
        tp, fn, fp, tn = (counts[key] for key in ("TP", "FN", "FP", "TN"))
        result.append(dict(
            method=keys[0], unit_substitution_rate=keys[1], array_fragments=keys[2],
            coverage=keys[3], substitution_rate=keys[4], assembly_fraction=keys[5],
            family_observations=len(group), available=available,
            unavailable=counts["NA"], TP=tp, FN=fn, FP=fp, TN=tn,
            sensitivity=tp / (tp + fn) if tp + fn else "NA",
            false_positive_rate=fp / (fp + tn) if fp + tn else "NA",
            precision=tp / (tp + fp) if tp + fp else "NA",
        ))
    return result


def _hybrid_row(pair: dict[str, dict], low_depth_threshold: float) -> dict:
    single = pair["single_k21"]
    multi = pair["multik_loglinear"]
    selected = multi if multi["outcome"] != "NA" else single
    threshold = low_depth_threshold if float(single["estimated_haploid_depth"]) < LOW_DEPTH_CUTOFF else 0.6
    scored = score_estimate(
        family_id=single["family_id"],
        period=int(single["period"]),
        estimate=float(selected["read_estimated_copies"]),
        assembly_bp=float(single["assembly_predicted_bp"]),
        truth_ratio=float(single["truth_assembly_read_ratio"]),
        threshold=float(single["decision_threshold"]),
        decision_threshold=threshold,
        method="multik_fallback_depth_rule",
        fit_status=(f"multik:{multi['fit_status']}" if selected is multi
                    else "fallback_single_k21:multik_unavailable"),
    )
    scored.pop("family_id")
    context = {key: value for key, value in single.items() if key not in {
        "method", "read_estimated_copies", "read_estimated_bp",
        "predicted_assembly_read_ratio", "native_status", "predicted_underrepresented",
        "outcome", "fit_status", "decision_threshold",
    }}
    return {**context, **scored}


def select_low_depth_threshold(pairs: list[dict[str, dict]], training_seeds: set[int]) -> float:
    training = [pair for pair in pairs if int(pair["single_k21"]["seed"]) in training_seeds]
    if not training:
        raise ValueError("Threshold selection requires training genomes")
    baseline_fp = sum(pair["single_k21"]["outcome"] == "FP" for pair in training)
    candidates = []
    for threshold in LOW_DEPTH_THRESHOLD_GRID:
        counts = Counter(_hybrid_row(pair, threshold)["outcome"] for pair in training)
        if counts["FP"] <= baseline_fp:
            candidates.append((counts["TP"], -counts["FP"], threshold))
    if not candidates:
        raise ValueError("No low-depth threshold satisfies the baseline false-positive constraint")
    return max(candidates)[2]


def calibration_rows(pairs: list[dict[str, dict]], seeds: list[int]) -> tuple[float, list[dict]]:
    selected = select_low_depth_threshold(pairs, set(seeds))
    folds: list[tuple[str, set[int], set[int], float]] = [
        ("full_development", set(seeds), set(seeds), selected)
    ]
    for heldout in seeds:
        training = set(seeds) - {heldout}
        folds.append((f"leave_s{heldout}_out", training, {heldout},
                      select_low_depth_threshold(pairs, training)))
    rows = []
    for fold, training, evaluated, threshold in folds:
        evaluation = [pair for pair in pairs if int(pair["single_k21"]["seed"]) in evaluated]
        baseline = Counter(pair["single_k21"]["outcome"] for pair in evaluation)
        calibrated = Counter(_hybrid_row(pair, threshold)["outcome"] for pair in evaluation)
        rows.append(dict(
            fold=fold,
            training_seeds=",".join(map(str, sorted(training))),
            evaluation_seeds=",".join(map(str, sorted(evaluated))),
            low_depth_cutoff=LOW_DEPTH_CUTOFF,
            selected_low_depth_threshold=threshold,
            standard_depth_threshold=0.6,
            baseline_TP=baseline["TP"], baseline_FN=baseline["FN"],
            baseline_FP=baseline["FP"], baseline_TN=baseline["TN"],
            calibrated_TP=calibrated["TP"], calibrated_FN=calibrated["FN"],
            calibrated_FP=calibrated["FP"], calibrated_TN=calibrated["TN"],
        ))
    return selected, rows


def _frozen_model(config: dict) -> float:
    model = config.get("collapse_model")
    if (
        not isinstance(model, dict)
        or model.get("method") != "multik_loglinear_else_single_k21"
        or model.get("k_values") != list(DEFAULT_K_VALUES)
        or model.get("low_depth_cutoff") != LOW_DEPTH_CUTOFF
        or model.get("standard_depth_threshold") != 0.6
        or model.get("unavailable_rule") != "fallback_single_k21"
        or model.get("calibration_seeds") != config.get("seeds", {}).get("development")
        or any(not re.fullmatch(r"[0-9a-f]{64}", str(model.get(name, ""))) for name in (
            "calibration_validation_sha256", "calibration_table_sha256",
            "calibration_environment_sha256",
        ))
    ):
        raise ValueError("Held-out validation requires the predeclared frozen collapse model")
    threshold = model.get("low_depth_threshold")
    if threshold not in LOW_DEPTH_THRESHOLD_GRID:
        raise ValueError("Frozen low-depth threshold is outside the declared calibration grid")
    return float(threshold)


def run(previous: Path, outdir: Path, backend: str = "rust", split: str = "development") -> None:
    previous = previous.resolve()
    prior_validation = json.loads((previous / "validation.json").read_text())
    prior_environment = json.loads((previous / "environment.json").read_text())
    config = json.loads((previous / "run_config.json").read_text())
    if (
        prior_validation.get("complete") is not True
        or prior_validation.get("successful") != prior_validation.get("executions")
        or prior_environment.get("split") != split
        or digest_file(previous / "run_config.json") != prior_environment.get("config_sha256")
        or backend not in {"python", "rust"}
        or split not in {"development", "heldout"}
    ):
        raise ValueError("Require a complete matching-split conditional abundance run")
    seeds = config.get("seeds", {}).get(split, [])
    other_seeds = [seed for name, group in config.get("seeds", {}).items() if name != split for seed in group]
    if not seeds or any(seed in other_seeds for seed in seeds):
        raise ValueError("Requested and other seed groups must be nonempty and disjoint")
    threshold = float(config["collapse_threshold"])
    if not 0 < threshold < 1:
        raise ValueError("Invalid collapse threshold")
    scenarios = challenge_scenarios(config)
    fragment_gap_bp = config.get("fragment_gap_bp", 0)

    baseline_rows = read_table(previous / "copy_number_metrics.tsv", {
        "seed", "coverage", "substitution_rate", "family_id", "estimate"
    })
    baseline = {
        (int(row["seed"]), float(row.get("unit_substitution_rate", 0)),
         int(row.get("array_fragments", 1)), float(row["coverage"]),
         float(row["substitution_rate"]), row["family_id"]):
        float(row["estimate"])
        for row in baseline_rows
    }
    locations = {
        (int(row["seed"]), float(row.get("unit_substitution_rate", 0)),
         int(row.get("array_fragments", 1)), float(row["assembly_fraction"]),
         row["family_id"]):
        float(row["predicted_assembly_bp"])
        for row in read_table(previous / "localization_metrics.tsv", {
            "seed", "assembly_fraction", "family_id", "predicted_assembly_bp"
        })
    }
    truth_rows = read_table(previous / "comparison_metrics.tsv", {
        "seed", "coverage", "substitution_rate", "assembly_fraction", "family_id",
        "truth_assembly_read_ratio", "truth_underrepresented"
    })
    truth = {
        (int(row["seed"]), float(row.get("unit_substitution_rate", 0)),
         int(row.get("array_fragments", 1)), float(row["coverage"]),
         float(row["substitution_rate"]), float(row["assembly_fraction"]),
         row["family_id"]): row
        for row in truth_rows
    }
    if len(baseline) != prior_validation["copy_number_family_rows"] or len(truth) != prior_validation["comparison_family_rows"]:
        raise ValueError("Prior root metrics contain duplicate or missing keys")

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir / "source_snapshot"
    provenance = source_manifest(root, snapshot)
    benchmark_hashes = {}
    for path in sorted(Path(__file__).parent.glob("*.py")):
        relative = path.relative_to(root)
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        benchmark_hashes[str(relative)] = digest_file(target)
    provenance.update(
        benchmark_source_sha256=benchmark_hashes,
        previous_result=str(previous),
        previous_validation_sha256=digest_file(previous / "validation.json"),
        previous_environment_sha256=digest_file(previous / "environment.json"),
        previous_config_sha256=digest_file(previous / "run_config.json"),
        split=split,
        challenge_scenarios=[dict(unit_substitution_rate=rate, array_fragments=count,
                                  fragment_gap_bp=fragment_gap_bp)
                             for rate, count in scenarios],
        methods=list(METHODS),
        k_values=list(DEFAULT_K_VALUES),
        collapse_threshold=threshold,
        scope=("development_known_catalogue_point_estimate_calibration"
               if split == "development" else
               "predeclared_heldout_known_catalogue_point_estimate_validation"),
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    shutil.copyfile(previous / "run_config.json", outdir / "run_config.json")

    rows = []
    multik_conditions = 0
    for seed in seeds:
        for scenario_index, (unit_rate, fragment_count) in enumerate(scenarios, 1):
            genome_dir = scenario_directory(
                previous / "genomes" / f"s{seed}", scenario_index, scenarios, config
            )
            run_dir = scenario_directory(
                previous / "runs" / f"s{seed}", scenario_index, scenarios, config
            )
            reads_root = scenario_directory(
                previous / "reads" / f"s{seed}", scenario_index, scenarios, config
            )
            genome_manifest = json.loads((genome_dir / "manifest.json").read_text())
            spec = genome_manifest.get("spec", {})
            if (float(spec.get("unit_substitution_rate", 0)) != unit_rate
                    or int(spec.get("array_fragments", 1)) != fragment_count
                    or int(spec.get("fragment_gap_bp", 0)) != fragment_gap_bp):
                raise ValueError("Prior genome scenario differs from configuration")
            for name in ("catalogue.fa", "truth_copy_number.tsv"):
                if digest_file(genome_dir / name) != genome_manifest["files"][name]:
                    raise ValueError(f"Prior genome input differs: s{seed}/{name}")
            catalogue = list(read_monomer_fasta(genome_dir / "catalogue.fa"))
            periods = {record.family_id: len(record.sequence) for record in catalogue}
            for coverage in config["coverages"]:
                for error in config["substitution_rates"]:
                    reads_dir = reads_root / f"c{coverage}_e{error}"
                    read_manifest = json.loads((reads_dir / "manifest.json").read_text())
                    if digest_file(reads_dir / "reads.fa") != read_manifest["files"]["reads.fa"]:
                        raise ValueError("Prior read input hash differs")
                    result = estimate_multik(
                        (record.sequence for record in read_fasta(reads_dir / "reads.fa")),
                        catalogue,
                        int(genome_manifest["genome_bp"]),
                        backend=backend,
                    )
                    multik_conditions += 1
                    by_family = {row["family_id"]: row for row in result.estimates}
                    if set(by_family) != set(periods):
                        raise ValueError("Multi-k families differ from the prior catalogue")
                    native_rows = read_table(
                        run_dir / f"c{coverage}_e{error}" / "quantify" / "output" / "copy_number.tsv",
                        {"family_id", "estimated_copy_number", "haploid_depth"},
                    )
                    native = {row["family_id"]: row for row in native_rows}
                    if set(native) != set(periods):
                        raise ValueError("Prior single-k output families differ from the catalogue")
                    if any(
                        float(native[family_id]["estimated_copy_number"])
                        != baseline[(int(seed), unit_rate, fragment_count,
                                     float(coverage), float(error), family_id)]
                        for family_id in periods
                    ):
                        raise ValueError("Prior root and native single-k estimates differ")
                    for fraction in config["assembly_fractions"]:
                        for family_id, period in periods.items():
                            key = (int(seed), unit_rate, fragment_count, float(coverage),
                                   float(error), float(fraction), family_id)
                            truth_row = truth[key]
                            context = dict(
                                seed=seed, unit_substitution_rate=unit_rate,
                                array_fragments=fragment_count, fragment_gap_bp=fragment_gap_bp,
                                coverage=coverage, substitution_rate=error,
                                assembly_fraction=fraction, family_id=family_id,
                                truth_assembly_read_ratio=float(truth_row["truth_assembly_read_ratio"]),
                                truth_underrepresented=truth_row["truth_underrepresented"],
                                assembly_predicted_bp=locations[(int(seed), unit_rate,
                                                                 fragment_count, float(fraction),
                                                                 family_id)],
                                period=period,
                                estimated_haploid_depth=float(native[family_id]["haploid_depth"]),
                            )
                            for method, estimate, fit_status in (
                                ("single_k21", baseline[(int(seed), unit_rate, fragment_count,
                                                         float(coverage), float(error), family_id)],
                                 "single_k_baseline"),
                                ("multik_loglinear", by_family[family_id]["extrapolated_copy_number"],
                                 by_family[family_id]["status"]),
                            ):
                                scored = score_estimate(
                                    family_id=family_id,
                                    period=period,
                                    estimate=estimate,
                                    assembly_bp=context["assembly_predicted_bp"],
                                    truth_ratio=context["truth_assembly_read_ratio"],
                                    threshold=threshold,
                                    method=method,
                                    fit_status=fit_status,
                                )
                                scored.pop("family_id")
                                rows.append({**context, **scored})
    pairs_by_key: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        key = (row["seed"], row["unit_substitution_rate"], row["array_fragments"],
               row["coverage"], row["substitution_rate"], row["assembly_fraction"],
               row["family_id"])
        pairs_by_key[key][row["method"]] = row
    pairs = list(pairs_by_key.values())
    if any(set(pair) != {"single_k21", "multik_loglinear"} for pair in pairs):
        raise ValueError("Incomplete paired single-k/multi-k results")
    if split == "development":
        selected_threshold, calibration = calibration_rows(pairs, [int(seed) for seed in seeds])
    else:
        selected_threshold = _frozen_model(config)
        calibration = []
    rows.extend(_hybrid_row(pair, selected_threshold) for pair in pairs)
    if split == "heldout":
        baseline_counts = Counter(pair["single_k21"]["outcome"] for pair in pairs)
        calibrated_counts = Counter(row["outcome"] for row in rows
                                    if row["method"] == "multik_fallback_depth_rule")
        model = config["collapse_model"]
        calibration.append(dict(
            fold="predeclared_heldout_validation",
            training_seeds=",".join(map(str, model["calibration_seeds"])),
            evaluation_seeds=",".join(map(str, seeds)),
            low_depth_cutoff=LOW_DEPTH_CUTOFF,
            selected_low_depth_threshold=selected_threshold,
            standard_depth_threshold=0.6,
            baseline_TP=baseline_counts["TP"], baseline_FN=baseline_counts["FN"],
            baseline_FP=baseline_counts["FP"], baseline_TN=baseline_counts["TN"],
            calibrated_TP=calibrated_counts["TP"], calibrated_FN=calibrated_counts["FN"],
            calibrated_FP=calibrated_counts["FP"], calibrated_TN=calibrated_counts["TN"],
        ))
    expected = prior_validation["comparison_family_rows"] * len(METHODS)
    summaries = summarize(rows)
    method_confusion = {}
    for method in METHODS:
        counts = Counter(row["outcome"] for row in rows if row["method"] == method)
        method_confusion[method] = {
            "available": sum(counts[name] for name in ("TP", "FN", "FP", "TN")),
            "unavailable": counts["NA"],
            "TP": counts["TP"],
            "FN": counts["FN"],
            "FP": counts["FP"],
            "TN": counts["TN"],
        }
    write_table(outdir / "comparison_metrics.tsv", rows, list(rows[0]))
    write_table(outdir / "comparison_summary.tsv", summaries, list(summaries[0]))
    write_table(outdir / "calibration.tsv", calibration, list(calibration[0]))
    validation = dict(
        complete=len(rows) == expected,
        split=split,
        independent_genomes=len(seeds),
        challenge_scenarios=len(scenarios),
        multik_read_conditions=multik_conditions,
        comparison_family_rows=len(rows),
        expected_comparison_family_rows=expected,
        unavailable=sum(row["outcome"] == "NA" for row in rows),
        method_confusion=method_confusion,
        selected_low_depth_threshold=selected_threshold,
        low_depth_cutoff=LOW_DEPTH_CUTOFF,
        leave_one_genome_out_folds=len(seeds) if split == "development" else 0,
        heldout_used=split == "heldout",
        scientific_acceptance="not_assumed_from_execution_success",
    )
    (outdir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["complete"]:
        raise RuntimeError("Incomplete multi-k collapse evaluation")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--backend", choices=("python", "rust"), default="rust")
    parser.add_argument("--split", choices=("development", "heldout"), default="development")
    args = parser.parse_args()
    run(args.previous, args.outdir, args.backend, args.split)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
