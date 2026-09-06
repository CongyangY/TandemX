"""Paired, alternating executions of archived and current TandemX sources."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import sys

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, write_table

OUTPUTS = ("candidate_reads.tsv", "candidate_monomers.fa", "monomers.fa", "families.tsv",
           "monomer_membership.tsv", "family_similarity.tsv")


def summarize_dataset(dataset: str, group: list[dict]) -> dict:
    valid = all(r["exit_code"] == 0 for r in group)
    same = valid and len({r["output_digest"] for r in group}) == 1
    row = {"dataset": dataset, "attempts": len(group), "all_successful": valid,
           "all_six_outputs_identical": same}
    for variant in ("baseline", "native_seed"):
        subset = [r for r in group if r["variant"] == variant]
        for metric in ("runtime_seconds", "peak_rss_mib", "cpu_user_seconds", "cpu_system_seconds"):
            row[f"{variant}_median_{metric}"] = statistics.median(r[metric] for r in subset) if valid else "NA"
    row["speedup_baseline_over_native"] = (row["baseline_median_runtime_seconds"] / row["native_seed_median_runtime_seconds"] if valid else "NA")
    row["rss_ratio_native_over_baseline"] = (row["native_seed_median_peak_rss_mib"] / row["baseline_median_peak_rss_mib"] if valid else "NA")
    return row


def compare(baseline: Path, datasets: list[Path], outdir: Path, repetitions: int) -> None:
    if repetitions < 2 or not (baseline / "tandemx").is_dir():
        raise ValueError("Need an archived TandemX source and at least two repetitions")
    outdir.mkdir(parents=True, exist_ok=False)
    provenance = source_manifest(Path(__file__).resolve().parents[2], outdir / "new_source_snapshot")
    baseline_provenance = source_manifest(baseline, outdir / "baseline_source_snapshot")
    provenance.update(baseline=baseline_provenance, script_sha256=digest_file(Path(__file__)),
                      timing_note="paired shuffled order; no profiler; same input and scientific parameters",
                      input_sha256={str(p / "reads.fa"): digest_file(p / "reads.fa") for p in datasets})
    origin = baseline.parent / "environment.json"
    if origin.is_file():
        provenance["baseline_origin_metadata"] = json.loads(origin.read_text())
    shutil.copyfile(Path(__file__), outdir / "runner_snapshot.py")
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    sources = {"baseline": outdir / "baseline_source_snapshot", "native_seed": outdir / "new_source_snapshot"}
    rows = []
    for dataset in datasets:
        for replicate in range(1, repetitions + 1):
            order = list(sources)
            random.Random(dataset.name + str(replicate)).shuffle(order)
            for label in order:
                folder = outdir / dataset.name / label / f"rep{replicate}"
                folder.mkdir(parents=True)
                command = [sys.executable, "-m", "tandemx.cli", "discover", "--reads", str(dataset / "reads.fa"),
                           "--outdir", str(folder / "discover"), "--min-period", "30", "--max-period", "1000",
                           "--min-repeat-span", "100", "--min-support-reads", "1", "--discovery-method", "elastic",
                           "--clustering-method", "sequence", "--cluster-identity", "0.95", "--kmer-backend", "rust",
                           "--threads", "1", "--no-progress"]
                (folder / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                measured = run_process(command, folder / "stdout.log", folder / "stderr.log", 120,
                                       {**os.environ, "PYTHONPATH": str(sources[label])}, sources[label])
                row = {"dataset": dataset.name, "variant": label, "replicate": replicate, **measured}
                if measured["exit_code"] == 0:
                    hashes = {name: digest_file(folder / "discover" / name) for name in OUTPUTS}
                    row["output_digest"] = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
                    (folder / "output_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
                rows.append(row)
                (folder / "receipt.json").write_text(json.dumps(row, indent=2) + "\n")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    write_table(outdir / "raw_runs.tsv", [{k: r.get(k, "NA") for k in fields} for r in rows], fields)
    summaries = [summarize_dataset(p.name, [r for r in rows if r["dataset"] == p.name]) for p in datasets]
    write_table(outdir / "summary.tsv", summaries, list(summaries[0]))
    complete = all(r["all_successful"] and r["all_six_outputs_identical"] for r in summaries)
    (outdir / "validation.json").write_text(json.dumps({"complete": True, "all_pairs_pass": complete,
                                                        "runs": len(rows), "datasets": len(datasets)}, indent=2) + "\n")
    if not complete:
        raise RuntimeError("Execution failure or changed output; see paired comparison records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--datasets", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    compare(args.baseline, args.datasets, args.outdir, args.repetitions)
