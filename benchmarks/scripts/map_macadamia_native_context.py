"""Map frozen Macadamia original reads to a preselected assembly-only context."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from benchmarks.scripts.score_rice_native_context import parse_paf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_binding(config: dict, repo: Path) -> tuple[Path, Path, list[Path]]:
    context = repo / config["context_fasta"]
    if sha256(context) != config["context_fasta_sha256"]:
        raise ValueError("selected context changed after freeze")
    if context.open().readline().split()[0] != ">" + config["context_name"]:
        raise ValueError("selected context name changed")
    audit = json.loads((repo / config["input_readback"]).read_text())
    if audit["status"] != "pass":
        raise ValueError("Macadamia source readback has not passed")
    by_name = {row["name"]: row for row in audit["results"]}
    sources = []
    for item in config["read_inputs"]:
        name = "macadamia_original_hifi_" + item["run"]
        if by_name[name]["status"] != "pass" or by_name[name]["sha256"] != item["sha256"]:
            raise ValueError(f"source readback changed: {name}")
        path = Path(item["path"])
        if by_name[name]["path"] != str(path) or not path.is_file():
            raise ValueError(f"source path changed: {name}")
        sources.append(path)
    executable = Path(config["minimap2_executable"])
    if sha256(executable) != config["minimap2_sha256"]:
        raise ValueError("minimap2 executable hash changed")
    version = subprocess.check_output([str(executable), "--version"], text=True).strip()
    if version != config["minimap2_version"]:
        raise ValueError("minimap2 version changed")
    return context, executable, sources


def score(paf: Path, config: dict, outdir: Path) -> dict:
    start_required = config["array_start_0_in_context"] - config["minimum_natural_flank_bp_each_side"]
    end_required = config["array_end_0_in_context"] + config["minimum_natural_flank_bp_each_side"]
    min_identity = config["minimum_nmatch_over_alignment_columns"]
    all_count = primary_count = geometric_count = 0
    good: dict[str, dict] = {}
    with gzip.open(paf, "rt") as stream:
        for line in stream:
            row = parse_paf(line)
            all_count += 1
            if not row["primary"]:
                continue
            primary_count += 1
            if (row["target"] != config["context_name"] or
                    row["target_start_0"] > start_required or
                    row["target_end_0"] < end_required):
                continue
            geometric_count += 1
            if row["identity"] < min_identity:
                continue
            prior = good.get(row["read_id"])
            if prior is None or (row["identity"], row["mapq"]) > (prior["identity"], prior["mapq"]):
                good[row["read_id"]] = row
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "qualifying_record_alignments.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["read_id", "read_length", "query_start_0", "query_end_0",
                         "strand", "target_start_0", "target_end_0", "matching_bases",
                         "alignment_columns", "identity", "mapq"])
        for row in sorted(good.values(), key=lambda value: value["read_id"]):
            writer.writerow([row[k] for k in ("read_id", "read_length", "query_start_0",
                                              "query_end_0", "strand", "target_start_0",
                                              "target_end_0", "matching_bases", "alignment_columns")]
                            + [f"{row['identity']:.9f}", row["mapq"]])
    zmws = []
    for read_id in good:
        match = re.match(r"^[^/]+/(\d+)/(?:ccs|\d+_\d+)$", read_id)
        if match is None:
            zmws = []
            break
        zmws.append(match.group(1))
    distinct_zmw = len(set(zmws)) if zmws and len(zmws) == len(good) else None
    threshold = config["minimum_distinct_original_records"]
    result = {
        "reported_alignment_count": all_count,
        "reported_primary_alignment_count": primary_count,
        "primary_geometric_spanner_alignment_count": geometric_count,
        "qualifying_distinct_record_count": len(good),
        "required_distinct_record_count": threshold,
        "distinct_zmw_count": distinct_zmw,
        "record_gate_pass": len(good) >= threshold,
        "molecule_gate_pass": distinct_zmw is not None and distinct_zmw >= threshold,
        "molecule_identity_boundary": None if distinct_zmw is not None else "original_ZMW_unverifiable_from_FASTQ_ids",
        "minimum_identity": min_identity,
        "minimum_natural_flank_bp_each_side": config["minimum_natural_flank_bp_each_side"],
        "paf_gzip_sha256": sha256(paf),
        "qualifying_tsv_sha256": sha256(outdir / "qualifying_record_alignments.tsv"),
    }
    (outdir / "score.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def run(protocol: Path, outdir: Path, repo: Path) -> dict:
    config = json.loads(protocol.read_text())
    context, executable, sources = verify_binding(config, repo)
    outdir.mkdir(parents=True, exist_ok=True)
    command = [str(executable), *config["minimap2_flags"], str(context), *(str(path) for path in sources)]
    paf = outdir / "context_alignments.paf.gz"
    partial = outdir / "context_alignments.paf.gz.partial"
    began = time.monotonic()
    with (outdir / "minimap2.stderr.txt").open("wb") as error_stream:
        with partial.open("wb") as raw_output:
            with gzip.GzipFile(fileobj=raw_output, mode="wb", mtime=0) as compressed:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error_stream)
                assert process.stdout is not None
                for block in iter(lambda: process.stdout.read(4 * 1024 * 1024), b""):
                    compressed.write(block)
                process.stdout.close()
                returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"minimap2 exited {returncode}; partial PAF retained")
    partial.replace(paf)
    result = score(paf, config, outdir)
    result.update({"command": command, "elapsed_seconds": round(time.monotonic() - began, 3),
                   "input_readback_sha256": sha256(repo / config["input_readback"]),
                   "protocol_sha256": sha256(protocol), "status": "completed"})
    (outdir / "run_receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.protocol, args.outdir, Path.cwd()), sort_keys=True))


if __name__ == "__main__":
    main()
