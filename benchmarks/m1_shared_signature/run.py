"""Run a small fixed-catalogue M1 development tournament on internal-disk data.

Example: conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.run \
  --outdir /tmp/tandemx-m1-dev
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

from . import model
from .model import align_and_gate, shared_signature_fit

BASES = "ACGT"
SCENARIOS = (
    ("related_low_error", 3, 0.03, 0.0),
    ("related_high_error", 3, 0.12, 0.0),
    ("family_error_shift", 3, 0.06, 0.20),
    ("identical_unidentifiable", 0, 0.03, 0.0),
)
PARAMETERS = dict(unit_bp=80, k=5, gate_fraction=0.22, foreground_reads=(320, 180),
                  unknown_reads=100, bootstrap_replicates=40, seeds=(11, 29, 47))


def mutate(sequence: str, positions: list[int], rng: random.Random) -> str:
    result = list(sequence)
    for position in positions:
        result[position] = rng.choice(BASES.replace(result[position], ""))
    return "".join(result)


def make_case(seed: int, differences: int, error: float, family_error_shift: float
              ) -> tuple[dict[str, str], list[str], list[str]]:
    """Simulate unit-sized circular reads and record their actual source labels."""
    rng = random.Random(seed)
    founder = "".join(rng.choices(BASES, k=PARAMETERS["unit_bp"]))
    positions = rng.sample(range(len(founder)), 8)
    related = mutate(founder, positions[:differences], rng)
    decoy = mutate(founder, positions, rng)
    catalogue = {"f1": founder, "f2": related, "decoy_zero": decoy}
    reads: list[str] = []
    labels: list[str] = []
    for family, count in (("f1", 320), ("f2", 180)):
        for _ in range(count):
            unit = catalogue[family]
            # Family-specific sequencing error is deliberately absent from the model.
            per_base_error = min(0.49, error + family_error_shift * (family == "f2"))
            noisy = mutate(unit, [i for i in range(len(unit)) if rng.random() < per_base_error], rng)
            shift = rng.randrange(len(unit))
            reads.append(noisy[shift:] + noisy[:shift])
            labels.append(family)
    for _ in range(PARAMETERS["unknown_reads"]):
        reads.append("".join(rng.choices(BASES, k=len(founder))))
        labels.append("unknown")
    paired = list(zip(reads, labels))
    rng.shuffle(paired)
    return catalogue, [read for read, _ in paired], [label for _, label in paired]


def evaluate(seed: int, scenario: tuple[str, int, float, float]) -> dict:
    name, differences, error, family_error_shift = scenario
    catalogue, reads, labels = make_case(seed, differences, error, family_error_shift)
    gate = round(PARAMETERS["unit_bp"] * PARAMETERS["gate_fraction"])
    start = time.perf_counter()
    accepted, mapped, rejected, ties = align_and_gate(reads, catalogue, gate)
    mapping_seconds = time.perf_counter() - start
    start = time.perf_counter()
    estimated, unresolved, groups = shared_signature_fit(accepted, catalogue, PARAMETERS["k"])
    regression_seconds = time.perf_counter() - start
    # Reuse the exact gate decisions to count negative-source leakage.
    accepted_by_label = {"f1": 0, "f2": 0, "unknown": 0}
    accepted_counts = Counter(accepted)
    for read, label in zip(reads, labels):
        if accepted_counts[read]:
            accepted_by_label[label] += 1
            accepted_counts[read] -= 1
    rng = random.Random(seed + 100_000)
    distributions = {key: [] for key in estimated if estimated[key] is not None}
    distributions.update({key: [] for key in unresolved})
    start = time.perf_counter()
    if accepted:
        for _ in range(PARAMETERS["bootstrap_replicates"]):
            sample = rng.choices(accepted, k=len(accepted))
            fit, group_fit, _ = shared_signature_fit(sample, catalogue, PARAMETERS["k"])
            for key in distributions:
                distributions[key].append(fit[key] if key in fit else group_fit[key])
    bootstrap_seconds = time.perf_counter() - start
    intervals = {key: [float(np.percentile(values, 2.5)),
                       float(np.percentile(values, 97.5))]
                 for key, values in distributions.items() if values}
    truth = {"f1": 320, "f2": 180, "decoy_zero": 0}
    score = {}
    for method, values in (("ordinary_mapping", mapped), ("shared_signature", estimated)):
        positives = [abs(values[key] / truth[key] - 1) for key in ("f1", "f2")
                     if values[key] is not None]
        score[method] = dict(identifiable_positive_mare=float(np.mean(positives)) if positives else None,
                             all_positive_mare_refusal_penalized=float(np.mean([
                                 1.0 if values[key] is None else abs(values[key] / truth[key] - 1)
                                 for key in ("f1", "f2")])),
                             individual_refused=sum(values[key] is None for key in truth),
                             zero_decoy_assigned_reads=values["decoy_zero"],
                             positive_interval_covered={key: intervals[key][0] <= truth[key] <= intervals[key][1]
                                                        for key in ("f1", "f2") if method == "shared_signature"
                                                        and key in intervals})
    digest = lambda value: hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
    return dict(seed=seed, scenario=name, truth_reads=truth, catalogue_sha256=digest(catalogue),
                reads_sha256=digest(reads), n_reads=len(reads), gate_mismatches=gate,
                accepted_reads=len(accepted), rejected_reads=rejected, tied_mapping_reads=ties,
                accepted_by_source=accepted_by_label, identifiable_groups=groups,
                unresolved_group_estimates=unresolved, ordinary_mapping=mapped,
                shared_signature=estimated, bootstrap_95pct=intervals, scores=score,
                seconds=dict(mapping=mapping_seconds, regression=regression_seconds,
                             bootstrap=bootstrap_seconds))


def run(outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    cases = [evaluate(seed, scenario) for scenario in SCENARIOS for seed in PARAMETERS["seeds"]]
    source_sha256 = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in (Path(__file__), Path(model.__file__))}
    summary = dict(status="development_only_no_heldout", parameters=PARAMETERS, scenarios=SCENARIOS,
                   cases=cases, source_sha256=source_sha256,
                   process_max_rss_raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   resource_note="macOS ru_maxrss reports bytes; timing includes Python/Numpy only; no T7 I/O")
    (outdir / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    run(parser.parse_args().outdir)
