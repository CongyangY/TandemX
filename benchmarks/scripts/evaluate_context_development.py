"""Small mechanism experiment; all scenarios are development, not holdout."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import statistics
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from benchmarks.challenge.context_prototype import (
    ContextConfig, PeriodicContextPrototype, partition_evidence,
)
from tandemx.utils.kmers import circular_kmer_counts, iter_linear_canonical_kmers


def development_scenarios() -> tuple[str, list[tuple[str, str, int]]]:
    rng = random.Random("context-development-v1-not-holdout")
    unit = "".join(rng.choice("ACGT") for _ in range(120))
    array = unit * 6
    flank = "".join(rng.choice("ACGT") for _ in range(80))
    substituted = list(array)
    for pos in rng.sample(range(len(array)), 7):
        substituted[pos] = rng.choice([base for base in "ACGT" if base != array[pos]])
    inserted = array
    for pos in sorted(rng.sample(range(len(array)), 7), reverse=True):
        inserted = inserted[:pos] + rng.choice("ACGT") + inserted[pos:]
    # Truth is assigned by construction, before inspecting any prediction.
    # Insertions within the planted repeat interval remain part of its bp truth.
    return unit, [
        ("complete_array", flank + array + flank, len(array)),
        ("shared_24bp_decoy", (unit[:24] + flank) * 60, 0),
        ("shared_80bp_decoy", (unit[:80] + flank) * 60, 0),
        ("array_plus_80bp_decoy", flank + array + flank + (unit[:80] + flank) * 60, len(array)),
        ("substitution_7_of_720", flank + "".join(substituted) + flank, len(array)),
        ("insertion_7_of_720", flank + inserted + flank, len(inserted)),
        ("partial_single_unit", flank + unit + flank, len(unit)),
    ]


def evaluate() -> dict[str, object]:
    unit, scenarios = development_scenarios()
    config = ContextConfig()
    prototype = PeriodicContextPrototype({"family": unit}, config)
    multiplicities = circular_kmer_counts(unit, config.k)
    rows = []
    for name, sequence, truth_bp in scenarios:
        counts = Counter(iter_linear_canonical_kmers(sequence, config.k))
        raw_kmer_bp = statistics.median(
            counts[word] / multiplicity for word, multiplicity in multiplicities.items()
        ) * len(unit)
        partition = partition_evidence(prototype.intervals(sequence), len(sequence))
        predicted = partition["unique_bp"].get("family", 0)
        rows.append({
            "scenario": name, "input_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "read_bp": len(sequence), "truth_repeat_bp": truth_bp,
            "raw_diagnostic_depth_bp": raw_kmer_bp,
            "raw_depth_absolute_error_bp": abs(raw_kmer_bp - truth_bp),
            "context_qualified_bp": predicted,
            "context_absolute_error_bp": abs(predicted - truth_bp),
            "ambiguous_bp": partition["ambiguous_bp"],
            "unassigned_bp": partition["unassigned_bp"],
        })
    return {
        "status": "development_mechanism_only_not_validation",
        "config": asdict(config), "unit": unit,
        "unit_sha256": hashlib.sha256(unit.encode()).hexdigest(),
        "rows": rows,
        "limitations": [
            "One synthetic unit; deliberately constructed decoys; all outcomes are development.",
            "Raw depth diagnostic is not a full frozen production quantify benchmark.",
            "No SRF or minimap2 occupancy comparison; no superiority claim.",
            "Exact phase breaks at indels; partial units deliberately fail the two-unit gate.",
            "No physical abundance, confidence calibration or cross-species validation.",
        ],
        "source_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), Path("benchmarks/challenge/context_prototype.py"))
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
