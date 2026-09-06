#!/usr/bin/env python3
"""Re-score archived predictions without rerunning or modifying any tool output."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import yaml

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays, score_families
from benchmarks.challenge.sequence_metrics import score_cyclic_recovery
from benchmarks.challenge.schema import ArrayRecord, digest_file, read_table, write_table
from benchmarks.challenge.run import source_manifest


def arrays(path: Path) -> list[ArrayRecord]:
    return [ArrayRecord(r["read_id"], int(r["start"]), int(r["end"]), int(r["period"]), r["sequence"], r["family_id"])
            for r in read_table(path)]


def rescore(runs: list[Path], outdir: Path) -> None:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Re-scoring requires a new empty output directory")
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    hashes = {}
    for source in runs:
        if not json.loads((source / "validation.json").read_text()).get("complete"):
            raise ValueError(f"Incomplete original run: {source}")
        raw = source / "raw_runs.tsv"
        hashes[str(raw)] = digest_file(raw)
        configuration = source / "run_config.yaml"
        config = yaml.safe_load(configuration.read_text())
        hashes[str(configuration)] = digest_file(configuration)
        for row in read_table(raw):
            # Explicit first-repetition snapshot; do not imply other repetitions agree.
            if row["repetition"] != "1":
                continue
            result = {"source_run": source.name, "scenario": row["scenario"], "seed": row["seed"],
                      "tool": row["tool"], "original_status": row["status"]}
            if row["status"] == "ok":
                dataset = source / "datasets" / row["dataset_id"]
                run = source / "runs" / row["dataset_id"] / row["tool"] / "rep1"
                truth_file, prediction_file = dataset / "truth_arrays.tsv", run / "predictions.tsv"
                read_file = dataset / "truth_reads.tsv"
                truth, prediction = arrays(truth_file), arrays(prediction_file)
                lengths = {r["read_id"]: int(r["length_bp"]) for r in read_table(read_file)}
                truth_units = {r.family_id: r.sequence for r in truth}
                files = [truth_file, prediction_file, read_file]
                if row["tool"] == "tandemx":
                    catalog = run / "discover" / "monomers.fa"
                    sequences = list(read_fasta(catalog).values())
                    files.append(catalog)
                else:
                    sequences = [r.sequence for r in prediction]
                result.update(score_arrays(prediction, truth, lengths, config["minimum_iou"])[0])
                result.update(score_families(sequences, truth_units)[0])
                result.update(score_cyclic_recovery(sequences, truth_units)[0])
                hashes.update({str(p): digest_file(p) for p in files})
            rows.append(result)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    write_table(outdir / "rescored_metrics.tsv", ({k: row.get(k, "NA") for k in fields} for row in rows), fields)
    manifest = source_manifest(Path(__file__).resolve().parents[2], outdir / "source_snapshot")
    manifest.update({"rescore_script_sha256": digest_file(Path(__file__)), "source_files": hashes,
                     "edlib_version": importlib.metadata.version("edlib"),
                     "scope": "first repetition; archived outputs unchanged; no new tool runtime measurements"})
    (outdir / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    rescore(args.runs, args.outdir)
