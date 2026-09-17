"""Run the frozen full-read development comparison without tuning."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

from tandemx.io.sequences import read_sequence_records_many
from tandemx.quantify.mvp import QuantifyConfig, quantify_toy_copy_number

from .full_read_research import FullReadConfig, classify_read, stream_files
from .generate_full_read_dev import PROTOCOL
from .model import align_and_gate

ROOT = Path(__file__).resolve().parents[2]
INPUT = Path(__file__).with_name("full_read_dev_v1")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_input() -> dict:
    manifest = json.loads((INPUT / "manifest.json").read_text())
    if digest(PROTOCOL) != manifest["protocol_sha256"]:
        raise ValueError("Frozen M1 full-read protocol changed")
    for name, expected in manifest["files"].items():
        if digest(INPUT / name) != expected:
            raise ValueError(f"Frozen M1 full-read input changed: {name}")
    return manifest


def ordinary_segments(sequence: str, catalogue: dict[str, str]) -> list[dict]:
    rows = []
    for start in range(0, len(sequence), 80):
        end = min(start + 80, len(sequence))
        window = sequence[start:end]
        if len(window) < 80 or set(window) - set("ACGT"):
            status, identity = "unknown", None
        else:
            accepted, mapped, rejected, tied = align_and_gate([window], catalogue, 18)
            if rejected:
                status, identity = "unknown", None
            elif tied:
                status, identity = "ambiguous", None
            else:
                assert accepted
                status, identity = "assigned", next(name for name, value in mapped.items() if value)
        rows.append(dict(start=start, end=end, status=status, identity=identity))
    return rows


def worker(method: str, outdir: Path) -> dict:
    verify_input()
    catalogue = {}
    from tandemx.quantify.mvp import read_monomer_fasta
    for item in read_monomer_fasta(INPUT / "catalogue.fa"):
        catalogue[item.family_id] = item.sequence
    reads = INPUT / "reads.fa"
    start = time.perf_counter()
    result: dict = {"method": method}
    if method == "full_read":
        settings = json.loads(PROTOCOL.read_text())["candidate_settings"]
        config = FullReadConfig(**settings)
        result["summary"] = stream_files(reads, INPUT / "catalogue.fa", config)
        result["wall_seconds"] = time.perf_counter() - start
        result["peak_process_rss_raw"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["segments"] = {
            record.id: classify_read(record.sequence, catalogue, config)
            for record in read_sequence_records_many([reads])
        }
        result["timing_scope"] = "one streaming inference pass; scoring trace generated after resource measurement"
    elif method == "ordinary_mapping":
        result["segments"] = {
            record.id: ordinary_segments(record.sequence, catalogue)
            for record in read_sequence_records_many([reads])
        }
        result["timing_scope"] = "one pass of 80-bp chunked mapping"
    elif method == "production_quantify":
        config = QuantifyConfig(reads=reads, monomers=INPUT / "catalogue.fa",
                                genome_size=20_000, outdir=outdir, k=21,
                                haploid_depth=1.0, kmer_backend="python")
        estimates = quantify_toy_copy_number(config)
        result["estimated_bp"] = {row.family_id: row.estimated_bp for row in estimates}
        result["confidence_warning"] = {
            row.family_id: {"confidence": row.confidence, "warning": row.warning}
            for row in estimates
        }
        result["timing_scope"] = "production quantify end-to-end function including output write"
    else:
        raise ValueError(method)
    if method != "full_read":
        result["wall_seconds"] = time.perf_counter() - start
        result["peak_process_rss_raw"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "worker.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def intersect(left: dict, right: dict) -> int:
    return max(0, min(left["end"], right["end"]) - max(left["start"], right["start"]))


def score_segments(segments: dict, truth: list[dict], catalogue: list[str]) -> dict:
    assigned: Counter[str] = Counter({name: 0 for name in catalogue})
    mass: Counter[str] = Counter()
    wrong = background_false = shared_ambiguous = unknown_background = unknown_positive = 0
    for read in truth:
        rows = segments[read["read_id"]]
        if rows[0]["start"] != 0 or rows[-1]["end"] != read["length_bp"]:
            raise ValueError("Uncovered read")
        if any(a["end"] != b["start"] for a, b in zip(rows, rows[1:])):
            raise ValueError("Overlapping or skipped bp")
        for row in rows:
            length = row["end"] - row["start"]
            mass[row["status"]] += length
            if row["status"] == "assigned":
                assigned[row["identity"]] += length
            for span in read["spans"]:
                bp = intersect(row, span)
                if row["status"] == "assigned" and row["identity"] != span["label"]:
                    wrong += bp
                    if span["label"] == "background":
                        background_false += bp
                if span["label"] == "twin_a" and row["status"] == "ambiguous":
                    shared_ambiguous += bp
                if row["status"] == "unknown":
                    if span["label"] == "background":
                        unknown_background += bp
                    else:
                        unknown_positive += bp
    return dict(assigned_by_family=dict(assigned), mass_by_status=dict(mass),
                wrong_family_assigned_bp=wrong, background_false_assigned_bp=background_false,
                shared_twin_ambiguous_bp=shared_ambiguous,
                unknown_background_bp=unknown_background,
                unknown_positive_bp=unknown_positive)


def run(outdir: Path) -> dict:
    manifest = verify_input()
    if outdir.exists():
        raise FileExistsError(outdir)
    outdir.mkdir(parents=True)
    methods = ("full_read", "ordinary_mapping", "production_quantify")
    outputs = {}
    for method in methods:
        target = outdir / method
        subprocess.run([sys.executable, "-m", "benchmarks.m1_shared_signature.run_full_read_dev", "--worker", method,
                        "--outdir", str(target)], cwd=ROOT, check=True)
        outputs[method] = json.loads((target / "worker.json").read_text())
    truth = json.loads((INPUT / "truth.json").read_text())
    catalogue = json.loads(PROTOCOL.read_text())["catalogue_families"]
    truth_bp: Counter[str] = Counter()
    for read in truth:
        for span in read["spans"]:
            truth_bp[span["label"]] += span["end"] - span["start"]
    positive = ("f1", "f2", "twin_a")
    rows = {}
    for method, output in outputs.items():
        scoring = (score_segments(output["segments"], truth, catalogue)
                   if "segments" in output else None)
        estimates = scoring["assigned_by_family"] if scoring else output["estimated_bp"]
        mare = (sum(abs(estimates[name] - truth_bp[name]) / truth_bp[name]
                    for name in positive) / len(positive) if scoring else None)
        rows[method] = dict(positive_family_bp_mare_with_abstentions_penalized=mare,
                            estimates_bp=estimates,
                            zero_decoy_assigned_bp=estimates["decoy_zero"],
                            wall_seconds=output["wall_seconds"],
                            peak_process_rss_raw=output["peak_process_rss_raw"],
                            segment_scoring=scoring,
                            endpoint="assigned_read_bp" if scoring else "uncalibrated_genome_scaled_estimated_bp",
                            confidence_warning=output.get("confidence_warning"),
                            timing_scope=output["timing_scope"])
    result = dict(status="development_only_no_backend_promotion", frozen_commit="9346de2",
                  input_manifest=manifest, truth_bp=dict(truth_bp), rows=rows,
                  resource_comparison_status="not_rankable_unequal_timing_scopes",
                  resource_note="Separate child processes on macOS; ru_maxrss raw bytes; equal input but unequal timed work; toy inputs only",
                  source_sha256={str(path.relative_to(ROOT)): digest(path) for path in (
                      Path(__file__), Path(__file__).with_name("full_read_research.py"),
                      Path(__file__).with_name("occupancy_research.py"),
                      Path(__file__).with_name("model.py"),
                      ROOT / "tandemx/quantify/mvp.py",
                      ROOT / "tandemx/io/sequences.py")})
    (outdir / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--worker", choices=("full_read", "ordinary_mapping", "production_quantify"))
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.outdir)
    else:
        run(args.outdir)
