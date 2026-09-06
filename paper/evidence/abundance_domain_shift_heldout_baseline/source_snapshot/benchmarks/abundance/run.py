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


def run(config_path: Path, outdir: Path, split: str) -> None:
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
                fragment_gap_bp=fragment_gap_bp,
            )
            located = []
            for variant in manifest["variants"]:
                folder = run_dir / variant["name"] / "locate"
                ok = execute("locate", folder, ["locate", "--assembly", str(genome_dir/f"{variant['name']}.fa"),
                             "--catalog", str(catalogue), "--k", str(config["k"])])
                retained = read_table(genome_dir/f"{variant['name']}.truth.tsv")
                if ok:
                    for row in score_localization(folder/"output"/"arrays.bed", retained, variant["genome_bp"]):
                        locate_scores.append(dict(seed=seed, **scenario_context,
                                                  assembly_fraction=variant["fraction"], **row))
                located.append((variant, retained, folder/"output"/"arrays.bed", ok))
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
    validation = dict(complete=True, executions=len(receipts), successful=sum(r["exit_code"]==0 for r in receipts),
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
    args = parser.parse_args()
    run(args.config, args.outdir, args.split)
