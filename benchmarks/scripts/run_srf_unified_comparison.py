"""Execute the preregistered SRF unified comparator; never an algorithm tuner."""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import ArrayRecord, digest_file, read_table, write_table
from benchmarks.challenge.simulate import Scenario
from benchmarks.challenge.unified_simulate import generate_unified_dataset
from benchmarks.challenge.unified_correspondence import match_native_catalogue
from benchmarks.challenge.unified_metrics import interval_metrics
from benchmarks.challenge.unified_native_io import (
    load_tandemx_arrays, load_tandemx_catalogue, mapping_paf_to_unique_arrays,
    mapping_paf_to_all_arrays,
)
from benchmarks.scripts.run_srf_pilot import workflow, parse_bed

SEEDS = (2026091001, 2026091002, 2026091003)
BASELINE = "81827c3e3fcd04bddbeba26c8165ac680f1b232e"
PROTOCOL = "docs/srf_formal_unified_protocol_20260910.md"
METHODS = ("tandemx", "srf_k151", "srf_k101", "competitive_mapping")


def conditions() -> tuple[Scenario, ...]:
    return (Scenario("clean"), Scenario("substitution", substitution_rate=.01),
            Scenario("indel", insertion_rate=.001, deletion_rate=.001),
            Scenario("divergence", unit_divergence=.02),
            Scenario("low_abundance", positive_fraction=.1),
            Scenario("shared_fragment", period=120, negative_kind="shared_fragment"))


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def stage(folder: Path, name: str, command: list[str], stdout: Path | None = None, *, deadline: float | None = None) -> dict:
    if deadline is not None and time.perf_counter() >= deadline:
        raise TimeoutError("Global launch deadline reached; remaining stages not_run")
    save(folder / f"{name}.command.json", command)
    result = run_process(command, stdout or folder / f"{name}.stdout.log",
                         folder / f"{name}.stderr.log", 180)
    save(folder / f"{name}.receipt.json", result)
    return {"stage": name, **result}


def resources(stages: list[dict], elapsed: float, bases: int) -> dict:
    return {"workflow_wall_seconds": elapsed,
            "native_stage_wall_seconds": sum(r["runtime_seconds"] for r in stages),
            "cpu_user_seconds": sum(r["cpu_user_seconds"] for r in stages),
            "cpu_system_seconds": sum(r["cpu_system_seconds"] for r in stages),
            "peak_rss_mib": max((r["peak_rss_mib"] for r in stages), default=0),
            "throughput_bp_per_second": bases / elapsed if elapsed else None,
            "memory_method": "max_direct_child_wait4_RSS_sequential_stages_controller_excluded",
            "peak_temporary_disk": "not_measured"}


def failed(stages: list[dict]) -> bool:
    return any(r["exit_code"] != 0 or r["timed_out"] for r in stages)


def abundance_metrics(native_bp: dict[str, float], correspondence: dict, truths: list[ArrayRecord],
                      truth_catalogue: dict[str, str]) -> dict:
    actual: dict[str, int] = defaultdict(int)
    for r in truths:
        actual[r.family_id] += r.end - r.start
    estimated = dict.fromkeys(truth_catalogue, 0.0)
    unassigned = 0.0
    for family, amount in native_bp.items():
        if amount < 0:
            raise ValueError("Negative native abundance")
        match = correspondence.get(family)
        if match is None:
            raise ValueError(f"Native abundance lacks catalogue entry: {family}")
        if match.status == "unique":
            estimated[match.matches[0]] += amount
        else:
            unassigned += amount
    per = [{"family_id": f, "truth_bp": actual[f], "estimated_bp": estimated[f],
            "absolute_relative_error": abs(estimated[f] - actual[f]) / actual[f] if actual[f] else None,
            "excess_bp": max(estimated[f] - actual[f], 0)} for f in truth_catalogue]
    errors = [r["absolute_relative_error"] for r in per if r["absolute_relative_error"] is not None]
    return {"MARE": sum(errors) / len(errors) if errors else None,
            "abundance_excess_bp": sum(r["excess_bp"] for r in per),
            "unassigned_abundance_bp": unassigned, "zero_truth_family_count": sum(actual[f] == 0 for f in truth_catalogue),
            "zero_truth_family_ids": ";".join(f for f in truth_catalogue if actual[f] == 0), "per_family": per}


def score(folder: Path, catalogue: dict[str, str], arrays: list[ArrayRecord], native_bp: dict[str, float],
          dataset: Path, discovery: bool = True) -> dict:
    truth_catalogue = read_fasta(dataset / "truth_monomers.fa")
    truths = [ArrayRecord(r["read_id"], int(r["start"]), int(r["end"]), int(r["period"]),
                         r["sequence"], r["family_id"]) for r in read_table(dataset / "truth_arrays.tsv")]
    lengths = {r["read_id"]: int(r["length_bp"]) for r in read_table(dataset / "truth_reads.tsv")}
    correspondence = match_native_catalogue(catalogue, truth_catalogue, .9, 64)
    native_to_truth = {f: v.matches[0] if v.status == "unique" else None for f, v in correspondence.items()}
    intervals = interval_metrics(arrays, truths, lengths, native_to_truth)
    quant = abundance_metrics(native_bp, correspondence, truths, truth_catalogue)
    present_truth = {r.family_id for r in truths}
    recovered = len({f for f in native_to_truth.values() if f} & present_truth)
    truth_read_ids = {r.read_id for r in truths}
    negative_bp = sum(n for rid, n in lengths.items() if rid not in truth_read_ids)
    save(folder / "correspondence.json", {k: asdict(v) for k, v in correspondence.items()})
    save(folder / "interval_metrics.json", intervals)
    save(folder / "native_abundance.json", native_bp)
    write_table(folder / "family_abundance.tsv", quant.pop("per_family"),
                ["family_id", "truth_bp", "estimated_bp", "absolute_relative_error", "excess_bp"])
    write_table(folder / "native_predictions.tsv", [asdict(r) for r in arrays], list(ArrayRecord.__dataclass_fields__))
    return {**intervals["global"], **quant, "family_recovery": recovered / len(present_truth) if discovery and present_truth else None,
            "family_recovery_state": "assessed" if discovery else "N/A_shared_TandemX_catalogue",
            "matched_truth_families": recovered, "truth_families": len(present_truth), "potential_truth_catalogue_families": len(truth_catalogue),
            "native_catalogue_count": len(catalogue), "negative_read_total_bp": negative_bp,
            "negative_read_attribution_fraction": intervals["global"]["negative_read_predicted_bp"] / negative_bp if negative_bp else None}


def run_tx(executable: Path, dataset: Path, folder: Path, bases: int, deadline: float) -> dict:
    folder.mkdir(parents=True)
    start = time.perf_counter()
    discover = folder / "discover"
    command = [str(executable), "discover", "--reads", str(dataset / "reads.fa"), "--outdir", str(discover),
               "--min-period", "30", "--max-period", "1000", "--min-repeat-span", "100",
               "--min-support-reads", "1", "--min-read-length", "1", "--discovery-method", "cascade",
               "--clustering-method", "sequence", "--kmer-backend", "rust", "--threads", "1", "--no-progress"]
    stages = [stage(folder, "discover", command, deadline=deadline)]
    catalogue, arrays, amounts = {}, [], {}
    state = "failed" if failed(stages) else "ok"
    if state == "ok" and (discover / "monomers.fa").stat().st_size:
        catalogue = load_tandemx_catalogue(discover / "monomers.fa")
        arrays = load_tandemx_arrays(discover)
        quant = [str(executable), "quantify", "--reads", str(dataset / "reads.fa"), "--catalogue",
                 str(discover / "monomers.fa"), "--outdir", str(folder / "quantify"), "--genome-size", str(bases),
                 "--haploid-depth", "1", "--k", "21", "--kmer-backend", "rust", "--no-progress"]
        stages.append(stage(folder, "quantify", quant, deadline=deadline))
        state = "failed" if failed(stages) else "ok"
        if state == "ok":
            amounts = {r["family_id"]: float(r["estimated_bp"]) for r in read_table(folder / "quantify/copy_number.tsv")}
    elif state == "ok":
        state = "no_catalogue"
        arrays = load_tandemx_arrays(discover)
    elapsed = time.perf_counter() - start
    result = {"method": "tandemx", "status": state, "stages": stages, **resources(stages, elapsed, bases)}
    if state != "failed":
        result.update(score(folder, catalogue, arrays, amounts, dataset))
    return result


def run_srf(tools: dict[str, Path], dataset: Path, folder: Path, bases: int, k: int, deadline: float) -> dict:
    measured = workflow(dataset / "reads.fa", folder, tools, 20, 180, k=k, deadline=deadline)
    result = {"method": f"srf_k{k}", "status": measured["status"], "stages": measured["stages"],
              **resources(measured["stages"], measured["workflow_wall_seconds"], bases)}
    if measured["status"] in ("ok", "no_catalogue", "no_eligible_kmers"):
        cat = read_fasta(folder / "srf.fa") if (folder / "srf.fa").exists() else {}
        arrays = parse_bed(folder / "srf.bed", cat) if measured["status"] == "ok" else []
        amounts = {}
        if measured["status"] == "ok":
            for line in (folder / "srf.abundance.tsv").read_text().splitlines():
                fields = line.split("\t")
                if len(fields) != 5 or fields[0] in amounts:
                    raise ValueError("Malformed/duplicate SRF native abundance row")
                amounts[fields[0]] = float(fields[1])
            bed_amounts: dict[str, int] = defaultdict(int)
            for r in arrays:
                bed_amounts[r.family_id] += r.end - r.start
            if any(amounts.get(f, 0) != amount for f, amount in bed_amounts.items()):
                raise ValueError("SRF retained native BED bp differs from native abundance")
        result.update(score(folder, cat, arrays, amounts, dataset))
    return result


def run_mapping(tools: dict[str, Path], dataset: Path, folder: Path, txfolder: Path, txreceipt: dict, bases: int, deadline: float) -> dict:
    folder.mkdir(parents=True)
    if txreceipt["status"] not in ("ok", "no_catalogue"):
        return {"method": "competitive_mapping", "status": "not_run_dependency_failure"}
    start = time.perf_counter()
    cat = load_tandemx_catalogue(txfolder / "discover/monomers.fa") if txreceipt["status"] != "no_catalogue" else {}
    stages, amounts, all_arrays, ambiguous = [], {}, [], 0
    if cat:
        template = folder / "templates.fa"
        with template.open("w") as handle:
            for family, seq in cat.items():
                handle.write(f">{family}\n{(seq * ((10000 + len(seq) - 1) // len(seq)))[:10000]}\n")
        command = [str(tools["minimap2"]), "-x", "map-hifi", "-c", "-N1000000", "-f1000", "-r100,100", "-t1",
                   str(template), str(dataset / "reads.fa")]
        stages.append(stage(folder, "map", command, folder / "mapping.paf", deadline=deadline))
        if not failed(stages):
            unique, ambiguous = mapping_paf_to_unique_arrays(folder / "mapping.paf", cat)
            all_arrays = mapping_paf_to_all_arrays(folder / "mapping.paf", cat)
            amounts = dict.fromkeys(cat, 0)
            for r in unique:
                amounts[r.family_id] += r.end - r.start
    elapsed = time.perf_counter() - start
    result = {"method": "competitive_mapping", "status": "failed" if failed(stages) else "ok" if cat else "no_catalogue",
              "stages": stages, "ambiguous_native_bp": ambiguous, **resources(stages, elapsed, bases)}
    shared = txreceipt["stages"][0]
    result.update(shared_discovery_plus_mapping_wall_seconds=shared["runtime_seconds"] + elapsed,
                  shared_discovery_plus_mapping_peak_rss_mib=max(shared["peak_rss_mib"], result["peak_rss_mib"]),
                  shared_discovery_plus_mapping_cpu_seconds=shared["cpu_user_seconds"] + shared["cpu_system_seconds"] + result["cpu_user_seconds"] + result["cpu_system_seconds"])
    if result["status"] != "failed":
        result.update(score(folder, cat, all_arrays, amounts, dataset, discovery=False))
    return result


def fingerprint(root: Path, tools: dict[str, Path], executable: Path) -> dict:
    files = [p for tree in ("tandemx", "benchmarks/challenge", "rust-core/src") for p in (root / tree).rglob("*")
             if p.is_file() and p.suffix in (".py", ".so", ".rs")]
    files += [root / "rust-core/Cargo.toml", root / "rust-core/Cargo.lock", root / PROTOCOL, Path(__file__).resolve(), root / "benchmarks/scripts/run_srf_pilot.py", executable]
    files += list(tools.values())
    return {str(p): digest_file(p) for p in sorted(set(files))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tandemx", type=Path, default=Path(sys.executable).with_name("tandemx"))
    parser.add_argument("--tools-root", type=Path, default=Path("/Volumes/T7/Codex/TandemX/tools/src"))
    parser.add_argument("--tool-build-receipt", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print(json.dumps({"seeds": SEEDS, "conditions": [asdict(c) for c in conditions()]}, indent=2))
        return
    root = Path(__file__).resolve().parents[2]
    tools = {"kmc": args.tools_root / "KMC/bin/kmc", "dump": args.tools_root / "KMC/bin/kmc_dump",
             "srf": args.tools_root / "srf/srf", "k8": args.tools_root / "k8-1.2/k8-arm64-Darwin",
             "minimap2": args.tools_root / "minimap2/minimap2", "utils": args.tools_root / "srf/srfutils.js"}
    if args.tool_build_receipt is None or not args.tool_build_receipt.is_file():
        raise ValueError("--tool-build-receipt is required for formal execution")
    signature = fingerprint(root, tools, args.tandemx)
    signature[str(args.tool_build_receipt)] = digest_file(args.tool_build_receipt)
    old_build = root / "paper/evidence/comparator_builds/provenance.json"
    signature[str(old_build)] = digest_file(old_build)
    if subprocess.check_output(["git", "diff", BASELINE, "--", "tandemx", "rust-core"], cwd=root):
        raise RuntimeError("Production differs from frozen baseline")
    if args.resume:
        env = json.loads((args.outdir / "environment.json").read_text())
        if env["fingerprint"] != signature:
            raise ValueError("Resume source/protocol/binary fingerprint mismatch")
    else:
        args.outdir.mkdir(parents=True, exist_ok=False)
        provenance = source_manifest(root, args.outdir / "source_snapshot")
        shutil.copyfile(args.tool_build_receipt, args.outdir / "clean_kmc_build_receipt.json")
        shutil.copyfile(old_build, args.outdir / "earlier_comparator_build_provenance.json")
        for name in (PROTOCOL, "benchmarks/scripts/run_srf_unified_comparison.py", "benchmarks/scripts/run_srf_pilot.py"):
            destination = args.outdir / "source_snapshot" / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / name, destination)
        save(args.outdir / "environment.json", {"fingerprint": signature, "provenance": provenance,
             "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
             "platform": platform.platform(), "cpu_count": os.cpu_count(), "load_average": os.getloadavg(),
             "python": sys.executable, "baseline": BASELINE, "started_unix": time.time(),
             "timing_scope": "serial_native_methods_three_simulation_seeds_no_technical_repetitions"})
    start = time.perf_counter()
    rows = []
    stop = False
    for condition in conditions():
        for seed in SEEDS:
            key = f"{condition.name}_s{seed}"
            dataset = args.outdir / "datasets" / key
            if not dataset.exists():
                if time.perf_counter() - start > 2700:
                    stop = True
                    break
                generate_unified_dataset(condition, seed, dataset)
            manifest = json.loads((dataset / "manifest.json").read_text())
            for name, item in manifest["files"].items():
                if digest_file(dataset / name) != item["sha256"]:
                    raise ValueError("Dataset checksum mismatch")
            bases = manifest["total_bases"]
            txreceipt = None
            for method in METHODS:
                folder = args.outdir / "runs" / key / method
                receipt = folder / "result.json"
                if receipt.exists():
                    result = json.loads(receipt.read_text())
                    for name, value in result["output_sha256"].items():
                        if digest_file(folder / name) != value:
                            raise ValueError("Existing run output checksum mismatch")
                else:
                    if time.perf_counter() - start > 2700:
                        stop = True
                        break
                    if folder.exists():
                        raise RuntimeError(f"Incomplete attempt retained; do not overwrite: {folder}")
                    print(f"START {key} {method}", flush=True)
                    try:
                        if method == "tandemx":
                            result = run_tx(args.tandemx, dataset, folder, bases, start + 2700)
                        elif method.startswith("srf_"):
                            result = run_srf(tools, dataset, folder, bases, int(method.split("k")[1]), start + 2700)
                        else:
                            result = run_mapping(tools, dataset, folder, folder.with_name("tandemx"), txreceipt, bases, start + 2700)
                    except TimeoutError as error:
                        result = {"method": method, "status": "not_run_deadline", "reason": str(error)}
                    except Exception as error:
                        folder.mkdir(parents=True, exist_ok=True)
                        save(folder / "adapter_error.json", {"error": repr(error), "status": "scoring_or_adapter_error"})
                        raise
                    result.update(condition=condition.name, seed=seed, input_sha256=digest_file(dataset / "reads.fa"))
                    result["retained_output_bytes"] = sum(p.stat().st_size for p in folder.rglob("*") if p.is_file())
                    result["output_sha256"] = {str(p.relative_to(folder)): digest_file(p) for p in folder.rglob("*") if p.is_file()}
                    save(receipt, result)
                    print(f"DONE {key} {method} {result['status']}", flush=True)
                if method == "tandemx":
                    txreceipt = result
                rows.append({k: v for k, v in result.items() if not isinstance(v, (dict, list))})
                fields = list(dict.fromkeys(k for row in rows for k in row))
                write_table(args.outdir / "summary.tsv", [{k: row.get(k) for k in fields} for row in rows], fields)
            if stop:
                break
        if stop:
            break
    completed = len(rows)
    seen = {(r["condition"], r["seed"], r["method"]) for r in rows}
    rows.extend({"condition": c.name, "seed": seed, "method": m, "status": "not_run"}
                for c in conditions() for seed in SEEDS for m in METHODS if (c.name, seed, m) not in seen)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    write_table(args.outdir / "summary.tsv", [{k: row.get(k) for k in fields} for row in rows], fields)
    save(args.outdir / "completion.json", {"completed_cells": completed, "expected_cells": 72,
         "unrun_cells": 72 - completed, "time_limit_reached": stop,
         "failures": sum(r["status"] == "failed" for r in rows), "elapsed_seconds": time.perf_counter() - start})


if __name__ == "__main__":
    main()
