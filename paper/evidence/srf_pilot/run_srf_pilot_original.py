"""Run SRF's documented KMC -> motif -> mapping -> abundance pipeline.

This bounded diagnostic keeps native outputs and all preparation costs. It is a
catalogue/abundance workflow, not a directly equivalent per-read finder. No truth
sequence or truth label is supplied to an external command.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays
from benchmarks.challenge.run import json_safe, run_process, source_manifest
from benchmarks.challenge.schema import ArrayRecord, digest_file, read_table, write_table
from benchmarks.challenge.sequence_metrics import score_cyclic_recovery


def parse_bed(path: Path, motifs: dict[str, str], keep_only: bool = True) -> list[ArrayRecord]:
    arrays = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 8 or fields[3] not in motifs or int(fields[7]) not in (0, 1, 2):
            raise ValueError(f"Malformed SRF native BED row: {line}")
        if keep_only and int(fields[7]) == 0:
            continue
        period = int(fields[6])
        if period != len(motifs[fields[3]]):
            raise ValueError("Native BED motif length differs from SRF catalogue")
        arrays.append(ArrayRecord(fields[0], int(fields[1]), int(fields[2]), period,
                                  motifs[fields[3]], fields[3]))
    return arrays


def workflow(reads: Path, outdir: Path, tools: dict[str, Path], minimum_count: int, timeout: float) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    temporary = outdir / "tmp"
    temporary.mkdir()
    start = time.perf_counter()
    with reads.open() as handle:
        total_bases = sum(len(line.strip()) for line in handle if not line.startswith(">"))
    if total_bases < 1:
        raise ValueError("Empty FASTA")
    db, dumped, catalog = outdir / "counts", outdir / "counts.txt", outdir / "srf.fa"
    elongated, paf, bed, abundance = [outdir / name for name in ("srf.enlong.fa", "srf.paf", "srf.bed", "srf.abundance.tsv")]
    kmc, dump, srf, k8, mm2, utils = [str(tools[name]) for name in ("kmc", "dump", "srf", "k8", "minimap2", "utils")]
    commands = [("count", [kmc, "-fm", "-k151", "-t1", "-m2", "-sm", f"-ci{minimum_count}", "-cs100000", str(reads), str(db), str(temporary)], None),
                ("dump", [dump, str(db), str(dumped)], None),
                ("assemble", [srf, "-p", "srf", str(dumped)], catalog),
                ("elongate", [k8, utils, "enlong", str(catalog)], elongated),
                ("map", [mm2, "-c", "-N1000000", "-f1000", "-r100,100", "-t1", str(elongated), str(reads)], paf),
                ("filter", [k8, utils, "paf2bed", str(paf)], bed),
                ("abundance", [k8, utils, "bed2abun", "-g", str(total_bases), str(bed)], abundance)]
    records = []
    status = "ok"
    for name, command, output in commands:
        if name == "elongate" and not read_fasta(catalog):
            status = "no_catalogue"
            break  # Native discovery completed; downstream stages explicitly not run.
        (outdir / f"{name}.command.json").write_text(json.dumps(command, indent=2) + "\n")
        measured = run_process(command, output or outdir / f"{name}.stdout.log", outdir / f"{name}.stderr.log", timeout)
        records.append({"stage": name, **measured})
        if measured["exit_code"] or measured["timed_out"]:
            status = "failed"
            break
    elapsed = time.perf_counter() - start
    record = {"status": status, "workflow_wall_seconds": elapsed,
              "external_stage_wall_seconds": sum(r["runtime_seconds"] for r in records),
              "maximum_external_stage_peak_rss_mib": max(r["peak_rss_mib"] for r in records),
              "memory_method": "max direct-child wait4 RSS across sequential native stages; benchmark controller excluded",
              "input_sha256": digest_file(reads), "total_input_bases": total_bases,
              "k": 151, "minimum_count": minimum_count, "threads_parameter": 1,
              "temporary_and_output_bytes_after_run": sum(p.stat().st_size for p in outdir.rglob("*") if p.is_file()),
              "stages": records, "skipped_stages": [name for name, _, _ in commands[len(records):]],
              "warning": "illustrative_readme_k151_count_preset;library_coverage_not_known;no_default_optimality_claim"}
    (outdir / "workflow_receipt.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def run(datasets: list[Path], tools_root: Path, outdir: Path, counts: list[int]) -> None:
    if any(c < 1 for c in counts):
        raise ValueError("Positive KMC count thresholds required")
    outdir.mkdir(parents=True, exist_ok=False)
    tools = {"kmc": tools_root / "KMC/bin/kmc", "dump": tools_root / "KMC/bin/kmc_dump",
             "srf": tools_root / "srf/srf", "k8": tools_root / "k8-1.2/k8-arm64-Darwin",
             "minimap2": tools_root / "minimap2/minimap2", "utils": tools_root / "srf/srfutils.js"}
    for path in tools.values():
        if not path.is_file():
            raise ValueError(f"Missing SRF dependency: {path}")
    provenance = source_manifest(Path(__file__).resolve().parents[2], outdir / "source_snapshot")
    provenance.update(executable_sha256={str(p): digest_file(p) for p in tools.values()},
                      script_sha256=digest_file(Path(__file__)))
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    rows = []
    for dataset in datasets:
        reads = dataset / "reads.fa"
        for count in counts:
            folder = outdir / dataset.name / f"ci{count}"
            measured = workflow(reads, folder, tools, count, 120)
            row = {"dataset": dataset.name, "minimum_count": count,
                   **{k: v for k, v in measured.items() if k not in ("stages", "skipped_stages")}}
            if measured["status"] in ("ok", "no_catalogue"):
                catalog = read_fasta(folder / "srf.fa")
                predictions = parse_bed(folder / "srf.bed", catalog) if measured["status"] == "ok" else []
                truths = [ArrayRecord(r["read_id"], int(r["start"]), int(r["end"]), int(r["period"]), r["sequence"], r["family_id"])
                          for r in read_table(dataset / "truth_arrays.tsv")]
                lengths = {r["read_id"]: int(r["length_bp"]) for r in read_table(dataset / "truth_reads.tsv")}
                metrics, _ = score_arrays(predictions, truths, lengths)
                row.update(metrics)
                row.update(score_cyclic_recovery(list(catalog.values()), {r.family_id: r.sequence for r in truths})[0])
                row["native_catalogue_size"] = len(catalog)
                row["native_motif_lengths"] = ";".join(str(len(s)) for s in catalog.values())
                write_table(folder / "predictions.tsv", [asdict(r) for r in predictions], list(ArrayRecord.__dataclass_fields__))
                row["truth_arrays_sha256"] = digest_file(dataset / "truth_arrays.tsv")
                row["truth_reads_sha256"] = digest_file(dataset / "truth_reads.tsv")
            rows.append(row)
            (folder / "metrics.json").write_text(json.dumps(json_safe(row), indent=2, allow_nan=False) + "\n")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    write_table(outdir / "summary.tsv", [{k: r.get(k, "NA") for k in fields} for r in rows], fields)
    receipt = {"complete": True, "runs": len(rows), "failed_runs": sum(r["status"] == "failed" for r in rows),
               "evidence_scope": "bounded SRF workflow pilot; no technical repetitions; no calibrated copy-number claim"}
    (outdir / "validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if receipt["failed_runs"]:
        raise RuntimeError("SRF pilot contains failures; see preserved logs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, nargs="+", required=True)
    parser.add_argument("--tools-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-counts", type=int, nargs="+", default=[100])
    args = parser.parse_args()
    run(args.datasets, args.tools_root, args.outdir, args.minimum_counts)
