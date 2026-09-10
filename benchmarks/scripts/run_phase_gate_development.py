"""Run all predeclared phase development conditions; never create holdout data."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import subprocess
import time
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path

from benchmarks.challenge.context_prototype import PeriodicContextPrototype, merge_intervals
from benchmarks.challenge.local_phase import LocalPhaseConfig, LocalPhasePrototype
from benchmarks.challenge.phase_gate_fixture import CONDITIONS, generate_case
from benchmarks.challenge.phase_gate_metrics import score_predictions
from tandemx.quantify.mvp import QuantifyConfig, quantify_toy_copy_number


def write_fasta(path: Path, records) -> None:
    with path.open("x") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def sequence_only(model, sequence: str) -> dict:
    hits = defaultdict(list)
    for pos in range(len(sequence) - model.config.k + 1):
        for family, _, _ in model.index.get(sequence[pos:pos + model.config.k], ()):
            hits[family].append((pos, pos + model.config.k))
    return {family: merge_intervals(rows) for family, rows in hits.items()}


def mapping_predictions(paf: Path) -> dict:
    predictions = defaultdict(lambda: defaultdict(list))
    for line in paf.read_text().splitlines():
        row = line.split("\t")
        if len(row) < 12:
            raise ValueError("Malformed PAF")
        start, end, matches, block = int(row[2]), int(row[3]), int(row[9]), int(row[10])
        if block >= 100 and matches / block >= 0.9:
            predictions[row[0]][row[5]].append((start, end))
    return predictions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--minimap2", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    configuration = LocalPhaseConfig()
    models = {
        "sequence_only": None,
        "plus_phase": replace(configuration, maximum_drift=0, phase_fraction=0, minimum_units=0),
        "plus_phase_coverage": replace(configuration, maximum_drift=0, minimum_units=0),
        "plus_recurrence_order": replace(configuration, maximum_drift=0),
        "A3": configuration,
        "A3_without_specificity_weights": replace(configuration, specificity_weights=False),
    }
    receipt = {
        "baseline": "81827c3", "split": "development", "algorithm_round": 1,
        "seeds": ["dev-phase-gate-01", "dev-phase-gate-02"],
        "conditions": CONDITIONS, "coverages": [2, 5, 10],
        "A3_config": asdict(configuration), "minimap2": str(args.minimap2),
        "minimap2_sha256": hashlib.sha256(args.minimap2.read_bytes()).hexdigest(),
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [
            Path(__file__), Path("benchmarks/challenge/local_phase.py"),
            Path("benchmarks/challenge/phase_gate_fixture.py"),
            Path("benchmarks/challenge/phase_gate_metrics.py")]},
        "endpoint": "sampled repeat bp; A0 haploid_depth=1, k=21; interval methods k=11 or native mapper",
        "classification": "synthetic source A is 0.3*truth for F1 and 1*truth for F2; R=read_estimate/observed_coverage; A/R<0.6",
        "status": "running",
    }
    (args.outdir / "protocol.json").write_text(json.dumps(receipt, indent=2) + "\n")
    rows = []
    for seed in receipt["seeds"]:
        for condition in CONDITIONS:
            for coverage in receipt["coverages"]:
                case = generate_case(seed, condition, coverage)
                root = args.outdir / f"{seed}_{condition}_{coverage}"
                root.mkdir()
                # Only observable sequences enter either external/production method.
                reads = root / "reads.fa"
                catalogue = root / "catalogue.fa"
                write_fasta(reads, case["reads"])
                write_fasta(catalogue, case["catalogue"].items())
                (root / "truth.json").write_text(json.dumps(case, indent=2) + "\n")
                start = time.perf_counter()
                quantify_toy_copy_number(QuantifyConfig(
                    reads=reads, monomers=catalogue, outdir=root / "A0", k=21,
                    haploid_depth=1.0, genome_size=case["genome_size"], kmer_backend="python"))
                a0_seconds = time.perf_counter() - start
                with (root / "A0/copy_number.tsv").open() as handle:
                    native = list(csv.DictReader(handle, delimiter="\t"))
                a0_bp = {row["family_id"]: float(row["estimated_copy_number"]) * len(case["catalogue"][row["family_id"]]) for row in native}
                target = root / "mapping_templates.fa"
                write_fasta(target, ((family, sequence * 20) for family, sequence in case["catalogue"].items()))
                command = [str(args.minimap2), "-x", "map-hifi", "-c", "-N", "50", "-f", "1000", "-r", "100,100", "-t", "1", str(target), str(reads)]
                start = time.perf_counter()
                with (root / "A1.paf").open("x") as stdout, (root / "A1.stderr").open("x") as stderr:
                    subprocess.run(command, stdout=stdout, stderr=stderr, check=True, timeout=60)
                a1_seconds = time.perf_counter() - start
                (root / "A1.command.json").write_text(json.dumps(command) + "\n")
                predictions_by_model = {"A1": (mapping_predictions(root / "A1.paf"), a1_seconds)}
                for name, cfg in {"A2": None, **models}.items():
                    start = time.perf_counter()
                    model = PeriodicContextPrototype(case["catalogue"]) if cfg is None else LocalPhasePrototype(case["catalogue"], cfg)
                    predictions = {rid: (sequence_only(model, seq) if name == "sequence_only" else model.intervals(seq)) for rid, seq in case["reads"]}
                    predictions_by_model[name] = (predictions, time.perf_counter() - start)
                lengths = {rid: len(seq) for rid, seq in case["reads"]}
                scores = {}
                for name, (predictions, seconds) in predictions_by_model.items():
                    score = score_predictions(predictions, case["truth_intervals"], lengths, set(case["catalogue"]))
                    scores[name] = score
                    for family, values in score["per_family"].items():
                        rows.append(dict(seed=seed, condition=condition, coverage=coverage, model=name, family=family,
                                         seconds=seconds, **values))
                for family, bp in a0_bp.items():
                    truth = case["read_truth_bp"].get(family, 0)
                    rows.append(dict(seed=seed, condition=condition, coverage=coverage, model="A0", family=family,
                                     seconds=a0_seconds, read_truth_bp=truth, predicted_unique_bp=bp,
                                     absolute_relative_error=abs(bp-truth)/truth if truth else None,
                                     recall=None, precision=None, background_false_bp=None))
                for row in rows[-(len(predictions_by_model) + 1)*2:]:
                    family = row["family"]
                    source_a = case["source_truth_bp"][family] * (0.3 if family == "F1" else 1.0)
                    estimate = row["predicted_unique_bp"] / (case["read_bases"] / case["genome_size"])
                    row.update(binary_truth=family == "F1", binary_call=(source_a / estimate < 0.6 if estimate else False),
                               excess_abundance_bp=max(0, row["predicted_unique_bp"] - row["read_truth_bp"]))
                (root / "scores.json").write_text(json.dumps(scores, indent=2) + "\n")
                print(f"completed {root.name}", flush=True)
    (args.outdir / "rows.json").write_text(json.dumps(rows, indent=2) + "\n")
    receipt["status"] = "complete"
    (args.outdir / "completion.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
