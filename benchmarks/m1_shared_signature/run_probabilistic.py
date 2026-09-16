"""Replay frozen M1 cases with a read-likelihood and unknown-rejection candidate.

Example: conda run -n tandemx-dev python -m \
  benchmarks.m1_shared_signature.run_probabilistic --outdir /tmp/tandemx-m1-prob
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import resource
import time

import numpy as np

from . import model, probabilistic, run as original

TRUTH = {"f1": 320, "f2": 180, "decoy_zero": 0}
SETTINGS = probabilistic.Settings()


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def score(values: dict[str, float | None]) -> dict:
    positive = [0.0 if values[name] is None else values[name] for name in ("f1", "f2")]
    truth = [320.0, 180.0]
    return dict(positive_mare_refusal_penalized=float(np.mean([
        abs(estimate / actual - 1) for estimate, actual in zip(positive, truth)])),
        positive_mae_reads=float(np.mean([abs(estimate - actual)
                                          for estimate, actual in zip(positive, truth)])),
        positive_bias_reads=float(np.mean([estimate - actual
                                           for estimate, actual in zip(positive, truth)])),
        individual_refused=sum(values[name] is None for name in ("f1", "f2")),
        zero_decoy_assigned_reads=values["decoy_zero"])


def evaluate_case(seed: int, scenario: tuple, prior: dict) -> dict:
    catalogue, reads, labels = original.make_case(seed, *scenario[1:])
    if (prior["reads_sha256"] != digest(reads)
            or prior["catalogue_sha256"] != digest(catalogue)
            or prior["truth_reads"] != TRUTH):
        raise ValueError("Frozen development inputs or truth changed")
    start = time.perf_counter()
    inference = probabilistic.infer(reads, catalogue, SETTINGS)
    fit_seconds = time.perf_counter() - start
    assignments = inference.pop("read_assignment")
    confusion = Counter()
    for truth, assignment in zip(labels, assignments):
        confusion[f"{truth}->{assignment}"] += 1
    negative_false_attribution = sum(assignment not in ("unknown", "gate", "ambiguous")
                                     for truth, assignment in zip(labels, assignments)
                                     if truth == "unknown")
    wrong_family_attribution = sum(assignment not in (truth, "unknown", "gate", "ambiguous", "f1+f2")
                                   for truth, assignment in zip(labels, assignments)
                                   if truth != "unknown")
    rejected_by_source = {source: sum(assignment in ("unknown", "gate", "ambiguous")
                                      for truth, assignment in zip(labels, assignments)
                                      if truth == source)
                          for source in ("f1", "f2", "unknown")}
    rng = random.Random(200_000 + seed)
    distributions = {name: [] for name, value in inference["family"].items() if value is not None}
    distributions.update({name: [] for name in inference["groups"]})
    start = time.perf_counter()
    for _ in range(original.PARAMETERS["bootstrap_replicates"]):
        sample = rng.choices(reads, k=len(reads))
        replicate = probabilistic.infer(sample, catalogue, SETTINGS)
        for name in distributions:
            distributions[name].append(replicate["family"].get(name, replicate["groups"].get(name)))
    bootstrap_seconds = time.perf_counter() - start
    intervals = {name: [float(np.percentile(values, 2.5)),
                        float(np.percentile(values, 97.5))]
                 for name, values in distributions.items()}
    coverage = {name: lo <= TRUTH[name] <= hi
                for name, (lo, hi) in intervals.items() if name in ("f1", "f2")}
    return dict(seed=seed, scenario=scenario[0], catalogue_sha256=digest(catalogue),
                reads_sha256=digest(reads), n_reads=len(reads), truth_reads=TRUTH,
                ordinary_mapping=prior["ordinary_mapping"],
                shared_signature=prior["shared_signature"],
                probabilistic=inference,
                ordinary_score=score(prior["ordinary_mapping"]),
                shared_score=score(prior["shared_signature"]),
                probabilistic_score=score(inference["family"]),
                probabilistic_confusion=dict(sorted(confusion.items())),
                negative_false_attribution=negative_false_attribution,
                wrong_family_attribution=wrong_family_attribution,
                rejected_by_source=rejected_by_source,
                bootstrap_95pct=intervals, bootstrap_positive_coverage=coverage,
                seconds=dict(ordinary_mapping=prior["seconds"]["mapping"],
                             shared_signature=prior["seconds"]["regression"],
                             probabilistic_fit=fit_seconds,
                             probabilistic_bootstrap=bootstrap_seconds))


def run(outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    prior_path = Path(__file__).with_name("evidence_20260916") / "results.json"
    prior = json.loads(prior_path.read_text())
    for source in (Path(original.__file__), Path(model.__file__)):
        if hashlib.sha256(source.read_bytes()).hexdigest() != prior["source_sha256"][source.name]:
            raise ValueError("Frozen ordinary/shared source changed; independent replay required")
    prior_cases = {(item["scenario"], item["seed"]): item for item in prior["cases"]}
    cases = [evaluate_case(seed, scenario, prior_cases[(scenario[0], seed)])
             for scenario in original.SCENARIOS for seed in original.PARAMETERS["seeds"]]
    sources = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in
               (Path(__file__), Path(probabilistic.__file__), Path(original.__file__), Path(model.__file__))}
    summary = dict(status="development_only_no_heldout", prior_results_sha256=hashlib.sha256(
        prior_path.read_bytes()).hexdigest(), parameters=original.PARAMETERS,
        settings=SETTINGS.__dict__, scenarios=original.SCENARIOS,
        source_sha256=sources, cases=cases,
        process_max_rss_raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        resource_note="macOS ru_maxrss bytes; one process; timing excludes toy generation; no T7 I/O")
    (outdir / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    run(parser.parse_args().outdir)
