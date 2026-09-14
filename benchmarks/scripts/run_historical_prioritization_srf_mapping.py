#!/usr/bin/env python3
"""Run the preregistered historical SRF/mapping endpoint without tuning."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.run import run_process
from benchmarks.challenge.schema import digest_file
from benchmarks.challenge.unified_correspondence import match_native_catalogue
from benchmarks.challenge.unified_native_io import load_tandemx_catalogue, mapping_paf_to_unique_arrays


ACK = "HISTORICAL_PROTOCOL_V1_REVIEWED"
METHODS = ("srf_k151", "srf_k101", "competitive_mapping")
PRIORITIZATION_FIELDS = (
    "family_id", "method", "reference_proxy_positive", "old_assembly_bp", "new_assembly_bp",
    "state", "native_retained_read_bp", "read_estimated_bp", "old_read_ratio",
    "predicted_proxy_positive", "outcome",
)
SNAPSHOT_FILES = (
    "docs/historical_prioritization_srf_mapping_protocol_v1_20260913.md",
    "benchmarks/scripts/run_historical_prioritization_srf_mapping.py",
    "benchmarks/challenge/unified_correspondence.py",
    "benchmarks/challenge/unified_native_io.py",
)


class CellFailure(RuntimeError):
    """A native stage or declared resource check failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def save_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def tree_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def parse_df_pk_free_bytes(output: str) -> int:
    """Parse POSIX ``df -Pk`` available KiB without using exFAT statvfs."""
    lines = [line for line in output.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("df -Pk returned no filesystem row")
    fields = lines[-1].split()
    if len(fields) < 6:
        raise ValueError("Malformed df -Pk filesystem row")
    available_kib = int(fields[3])
    if available_kib < 0:
        raise ValueError("df -Pk returned negative available space")
    return available_kib * 1024


def free_bytes(path: Path) -> tuple[int, str]:
    """Measure free bytes; macOS ``statvfs`` is unreliable on this exFAT volume."""
    if sys.platform == "darwin":
        output = subprocess.check_output(["df", "-Pk", str(path)], text=True)
        return parse_df_pk_free_bytes(output), "darwin_df_Pk_available_kib_x_1024"
    return shutil.disk_usage(path).free, "shutil_disk_usage_free"


def load_config(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError("Configuration must be a JSON object")
    return value


def _positive_int(value: object, label: str) -> int:
    if type(value) is not int or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def validate_config(config: Mapping[str, Any], root: Path, *, verify_files: bool) -> dict[str, Any]:
    """Validate frozen semantics; optionally perform the expensive file hashes."""
    if config.get("schema_version") != 1:
        raise ValueError("Unsupported schema_version")
    if config.get("status") != "frozen_before_formal_historical_native_output":
        raise ValueError("Configuration is not marked frozen before native output")
    guard = config.get("preregistration_guard", {})
    if guard.get("native_output_inspected_for_parameter_selection") is not False:
        raise ValueError("Native-output selection guard must be false")
    if guard.get("accept_preexisting_native_outputs") is not False:
        raise ValueError("Pre-existing native output import must be disabled")
    if guard.get("required_execution_ack") != ACK:
        raise ValueError("Unexpected execution acknowledgement")
    if tuple(config.get("methods", ())) != METHODS:
        raise ValueError(f"Methods/order must be {METHODS}")

    srf = config.get("srf", {})
    if srf.get("k_values") != [151, 101] or srf.get("minimum_count") != 20:
        raise ValueError("SRF k151/k101 and ci20 are frozen")
    if srf.get("threads") != 1 or srf.get("kmc_memory_gib") != 2 or srf.get("strict_memory") is not True:
        raise ValueError("SRF resource arguments differ from the frozen protocol")
    mapping = config.get("competitive_mapping", {})
    expected_mapping = {"template_bp": 10000, "preset": "map-hifi", "max_secondary": 1000000,
                        "minimizer_occurrence_floor": 1000, "chaining_bandwidth": "100,100",
                        "threads": 1, "minimum_alignment_block_bp": 100, "minimum_identity": 0.9,
                        "retain_primary_and_secondary": True, "cross_family_overlap": "exclude_and_report_once"}
    if mapping != expected_mapping:
        raise ValueError("Competitive-mapping rule differs from the frozen simulation rule")
    correspondence = config.get("correspondence", {})
    if correspondence != {"identity": 0.9, "max_integer_multiple": 64,
                           "length_tolerance_fraction": 0.1,
                           "unique_only": True, "redistribute_unassigned": False}:
        raise ValueError("Correspondence rule differs from the frozen rule")
    normalization = config.get("normalization", {})
    if normalization.get("formula") != "native_retained_read_bp/total_library_bases*genome_size_bp":
        raise ValueError("Unexpected normalization")
    evaluation = config.get("evaluation", {})
    if evaluation.get("collapse_threshold") != 0.6 or evaluation.get("minimum_new_assembly_bp") != 15000:
        raise ValueError("Historical endpoint threshold changed")
    if evaluation.get("zero_estimate_prediction") != "N/A_undefined_denominator":
        raise ValueError("Zero-estimate semantics changed")

    resources = config.get("resources", {})
    if resources.get("serial_cells") is not True or resources.get("minimum_free_bytes") != 250_000_000_000:
        raise ValueError("Unexpected cell scheduling or free-disk stop rule")
    _positive_int(resources.get("maximum_cell_bytes"), "maximum_cell_bytes")
    _positive_int(resources.get("global_launch_seconds"), "global_launch_seconds")
    for stage, seconds in resources.get("stage_timeout_seconds", {}).items():
        _positive_int(seconds, f"stage timeout {stage}")

    species = config.get("species")
    if not isinstance(species, list) or [row.get("id") for row in species] != ["ey15_2", "macadamia_jansenii"]:
        raise ValueError("Species/order must be ey15_2 then macadamia_jansenii")
    manifest: dict[str, Any] = {"verified": verify_files, "files": []}
    for row in species:
        reads = row.get("reads", [])
        expected_read_files = 1 if row["id"] == "ey15_2" else 2
        if len(reads) != expected_read_files:
            raise ValueError(f"{row['id']} requires {expected_read_files} FASTQ file(s)")
        if row["id"] == "macadamia_jansenii" and [r["run_accession"] for r in reads] != ["SRR13557763", "SRR13557762"]:
            raise ValueError("Macadamia FASTQ order changed")
        if sum(_positive_int(r["read_count"], "read_count") for r in reads) != row["total_read_count"]:
            raise ValueError(f"Read-count total mismatch for {row['id']}")
        if sum(_positive_int(r["total_bases"], "total_bases") for r in reads) != row["total_library_bases"]:
            raise ValueError(f"Base total mismatch for {row['id']}")
        _positive_int(row["genome_size_bp"], "genome_size_bp")
        for descriptor in [*reads, row["catalogue"], row["reference_metrics"], row["source_config"]]:
            path = Path(descriptor["path"])
            if not path.is_absolute():
                path = root / path
            item = {"path": str(path), "expected_bytes": descriptor["bytes"],
                    "expected_sha256": descriptor["sha256"]}
            if verify_files:
                if not path.is_file():
                    raise FileNotFoundError(path)
                item.update(observed_bytes=path.stat().st_size, observed_sha256=digest_file(path))
                if item["observed_bytes"] != item["expected_bytes"] or item["observed_sha256"] != item["expected_sha256"]:
                    raise ValueError(f"Input binding mismatch: {path}")
            manifest["files"].append(item)
        if verify_files:
            catalogue = load_tandemx_catalogue(_resolve(root, row["catalogue"]["path"]))
            if len(catalogue) != row["catalogue"]["records"]:
                raise ValueError(f"Catalogue record count mismatch for {row['id']}")
            reference = read_tsv(_resolve(root, row["reference_metrics"]["path"]))
            validate_reference(reference, row)

    protocol = config.get("protocol", {})
    if protocol.get("path") not in SNAPSHOT_FILES:
        raise ValueError("Unexpected protocol path")
    if verify_files:
        protocol_path = root / protocol["path"]
        protocol_item = {"path": str(protocol_path), "expected_bytes": protocol["bytes"],
                         "expected_sha256": protocol["sha256"], "observed_bytes": protocol_path.stat().st_size,
                         "observed_sha256": digest_file(protocol_path)}
        if (protocol_item["observed_bytes"] != protocol_item["expected_bytes"] or
                protocol_item["observed_sha256"] != protocol_item["expected_sha256"]):
            raise ValueError("Protocol SHA-256 mismatch")
        manifest["files"].append(protocol_item)
        upstream = config["upstream_simulation_protocol"]
        upstream_path = root / upstream["path"]
        upstream_item = {"path": str(upstream_path), "expected_sha256": upstream["sha256"],
                         "observed_bytes": upstream_path.stat().st_size,
                         "observed_sha256": digest_file(upstream_path)}
        if upstream_item["observed_sha256"] != upstream_item["expected_sha256"]:
            raise ValueError("Upstream simulation protocol SHA-256 mismatch")
        manifest["files"].append(upstream_item)
        for descriptor in config.get("implementation_sources", []):
            path = root / descriptor["path"]
            item = {"path": str(path), "expected_sha256": descriptor["sha256"],
                    "observed_bytes": path.stat().st_size, "observed_sha256": digest_file(path)}
            if item["observed_sha256"] != item["expected_sha256"]:
                raise ValueError(f"Implementation-source binding mismatch: {path}")
            manifest["files"].append(item)
        for name, tool in config.get("tools", {}).items():
            path = Path(tool["path"])
            item = {"name": name, "path": str(path), "expected_bytes": tool["bytes"],
                    "expected_sha256": tool["sha256"], "observed_bytes": path.stat().st_size if path.is_file() else None,
                    "observed_sha256": digest_file(path) if path.is_file() else None}
            if (not path.is_file() or item["observed_bytes"] != item["expected_bytes"] or
                    item["observed_sha256"] != item["expected_sha256"]):
                raise ValueError(f"Tool binding mismatch: {name}")
            manifest["files"].append(item)
    return manifest


def _resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def validate_reference(rows: list[dict[str, str]], species: Mapping[str, Any]) -> list[dict[str, str]]:
    required = {"family_id", "old_assembly_bp", "new_assembly_bp", "old_new_ratio",
                "reference_state", "eligibility"}
    if not rows or required - rows[0].keys():
        raise ValueError("Reference metrics lacks required fields")
    eligible = []
    seen: set[str] = set()
    for row in rows:
        family = row["family_id"]
        if family in seen:
            raise ValueError(f"Duplicate reference family: {family}")
        seen.add(family)
        old, new = float(row["old_assembly_bp"]), float(row["new_assembly_bp"])
        expected_eligible = new >= 15000
        if (row["eligibility"] == "eligible") != expected_eligible:
            raise ValueError(f"Eligibility mismatch for {family}")
        if expected_eligible:
            expected_state = "reference_collapse" if old / new < 0.6 else "reference_retained"
            if row["reference_state"] != expected_state:
                raise ValueError(f"Reference state mismatch for {family}")
            if abs(float(row["old_new_ratio"]) - old / new) > 1e-12:
                raise ValueError(f"Old/new ratio mismatch for {family}")
            eligible.append(row)
    if len(eligible) != species["reference_metrics"]["eligible_records"]:
        raise ValueError(f"Eligible-family count mismatch for {species['id']}")
    return eligible


def plan_rows(config: Mapping[str, Any]) -> list[dict[str, object]]:
    return [{"cell_index": index, "species": species["id"], "method": method,
             "fastq_count": len(species["reads"]), "total_library_bases": species["total_library_bases"]}
            for index, (species, method) in enumerate(
                ((s, m) for s in config["species"] for m in METHODS), 1)]


def kmc_command(config: Mapping[str, Any], reads_list: Path, cell: Path, k: int) -> list[str]:
    srf = config["srf"]
    return [config["tools"]["kmc"]["path"], "-fq", f"-k{k}", "-t1", "-m2", "-sm",
            f"-ci{srf['minimum_count']}", f"-cs{srf['maximum_count']}", f"@{reads_list}",
            str(cell / "counts"), str(cell / "tmp")]


def srf_commands(config: Mapping[str, Any], species: Mapping[str, Any], cell: Path, k: int) -> list[tuple[str, list[str], Path | None]]:
    tools = config["tools"]
    reads = [item["path"] for item in species["reads"]]
    return [
        ("count", kmc_command(config, cell / "reads.list", cell, k), None),
        ("dump", [tools["kmc_dump"]["path"], str(cell / "counts"), str(cell / "counts.txt")], None),
        ("assemble", [tools["srf"]["path"], "-p", "srf", str(cell / "counts.txt")], cell / "srf.fa"),
        ("elongate", [tools["k8"]["path"], tools["srfutils"]["path"], "enlong", str(cell / "srf.fa")], cell / "srf.enlong.fa"),
        ("map", [tools["minimap2"]["path"], "-c", "-N1000000", "-f1000", "-r100,100", "-t1",
                 str(cell / "srf.enlong.fa"), *reads], cell / "srf.paf"),
        ("filter", [tools["k8"]["path"], tools["srfutils"]["path"], "paf2bed", str(cell / "srf.paf")], cell / "srf.bed"),
        ("abundance", [tools["k8"]["path"], tools["srfutils"]["path"], "bed2abun", "-g",
                       str(species["total_library_bases"]), str(cell / "srf.bed")], cell / "srf.abundance.tsv"),
    ]


def mapping_command(config: Mapping[str, Any], species: Mapping[str, Any], cell: Path) -> list[str]:
    reads = [item["path"] for item in species["reads"]]
    return [config["tools"]["minimap2"]["path"], "-x", "map-hifi", "-c", "-N1000000",
            "-f1000", "-r100,100", "-t1", str(cell / "templates.fa"), *reads]


def make_templates(catalogue: Mapping[str, str], path: Path, length: int = 10000) -> None:
    with path.open("w") as handle:
        for family, sequence in catalogue.items():
            template = (sequence * ((length + len(sequence) - 1) // len(sequence)))[:length]
            handle.write(f">{family}\n{template}\n")


def _stage_outputs(cell: Path, stage_name: str, stdout: Path | None) -> list[Path]:
    if stage_name == "count":
        return [cell / "counts.kmc_pre", cell / "counts.kmc_suf"]
    if stage_name == "dump":
        return [cell / "counts.txt"]
    return [stdout] if stdout else []


def execute_stage(config: Mapping[str, Any], cell: Path, name: str, command: list[str], stdout: Path | None,
                  deadline: float) -> dict[str, object]:
    timeout = config["resources"]["stage_timeout_seconds"][name]
    remaining = deadline - time.perf_counter()
    if remaining <= 0:
        raise CellFailure("global_launch_deadline")
    timeout = min(timeout, remaining)
    save_json(cell / f"{name}.command.json", command)
    started = utc_now()
    measured = run_process(command, stdout or cell / f"{name}.stdout.log", cell / f"{name}.stderr.log", timeout)
    outputs = []
    for path in _stage_outputs(cell, name, stdout):
        outputs.append({"path": str(path), "exists": path.is_file(),
                        "bytes": path.stat().st_size if path.is_file() else None,
                        "sha256": digest_file(path) if path.is_file() else None})
    usage = tree_bytes(cell)
    free, free_method = free_bytes(cell)
    receipt = {"stage": name, "command": command, "started_utc": started, "ended_utc": utc_now(),
               **measured, "outputs": outputs, "cell_bytes_after_stage": usage,
               "free_bytes_after_stage": free, "free_bytes_method": free_method,
               "disk_measurement": "post_stage_not_continuous_peak"}
    save_json(cell / f"{name}.receipt.json", receipt)
    if measured["exit_code"] != 0 or measured["timed_out"]:
        raise CellFailure(f"stage_failed:{name}")
    if usage > config["resources"]["maximum_cell_bytes"]:
        raise CellFailure(f"maximum_cell_bytes_after:{name}")
    if free < config["resources"]["minimum_free_bytes"]:
        raise CellFailure(f"minimum_free_bytes_after:{name}")
    return receipt


def parse_srf_abundance(path: Path, bed: Path, catalogue: Mapping[str, str]) -> dict[str, int]:
    amounts: dict[str, int] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 5 or fields[0] not in catalogue or fields[0] in amounts:
            raise ValueError(f"Malformed/duplicate SRF abundance row {number}")
        value = int(fields[1])
        if value < 0:
            raise ValueError("Negative SRF abundance")
        amounts[fields[0]] = value
    bed_amounts: dict[str, int] = defaultdict(int)
    for number, line in enumerate(bed.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 8 or fields[3] not in catalogue or int(fields[7]) not in (0, 1, 2):
            raise ValueError(f"Malformed SRF BED row {number}")
        if int(fields[7]) > 0:
            bed_amounts[fields[3]] += int(fields[2]) - int(fields[1])
    if any(amounts.get(family, 0) != value for family, value in bed_amounts.items()):
        raise ValueError("SRF abundance differs from positive-flag BED bases")
    return amounts


def assign_srf_abundance(native: Mapping[str, int], correspondence: Mapping[str, Any],
                         catalogue: Mapping[str, str]) -> tuple[dict[str, int], list[dict[str, object]]]:
    assigned = dict.fromkeys(catalogue, 0)
    unassigned = []
    for native_id, amount in native.items():
        match = correspondence[native_id]
        if match.status == "unique":
            assigned[match.matches[0]] += amount
        else:
            unassigned.append({"native_id": native_id, "status": match.status,
                               "native_retained_read_bp": amount,
                               "candidate_matches": ";".join(match.matches)})
    return assigned, unassigned


def normalize_amounts(raw: Mapping[str, int], total_bases: int, genome_size: int) -> dict[str, float]:
    if total_bases <= 0 or genome_size <= 0:
        raise ValueError("Normalization denominators must be positive")
    return {family: amount / total_bases * genome_size for family, amount in raw.items()}


def score_estimates(reference: list[dict[str, str]], estimates: Mapping[str, float] | None,
                    method: str, technical_state: str = "ok") -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in reference:
        if row["eligibility"] != "eligible":
            continue
        truth = row["reference_state"] == "reference_collapse"
        base = {"family_id": row["family_id"], "method": method,
                "reference_proxy_positive": truth, "old_assembly_bp": float(row["old_assembly_bp"]),
                "new_assembly_bp": float(row["new_assembly_bp"])}
        if estimates is None:
            rows.append({**base, "state": technical_state, "native_retained_read_bp": "N/A",
                         "read_estimated_bp": "N/A", "old_read_ratio": "N/A",
                         "predicted_proxy_positive": "N/A", "outcome": "N/A"})
            continue
        estimate = estimates.get(row["family_id"], 0.0)
        if estimate < 0:
            raise ValueError("Negative normalized estimate")
        if estimate == 0:
            rows.append({**base, "state": "zero_native_support", "read_estimated_bp": 0.0,
                         "old_read_ratio": "N/A", "predicted_proxy_positive": "N/A", "outcome": "N/A"})
            continue
        ratio = float(row["old_assembly_bp"]) / estimate
        call = ratio < 0.6
        outcome = "TP" if truth and call else "FN" if truth else "FP" if call else "TN"
        rows.append({**base, "state": "ok", "read_estimated_bp": estimate,
                     "old_read_ratio": ratio, "predicted_proxy_positive": call, "outcome": outcome})
    counts = Counter(str(row["outcome"]) for row in rows)
    return rows, {"eligible_families": len(rows), "available_families": len(rows) - counts["N/A"],
                  "unavailable_families": counts["N/A"],
                  "zero_supported_eligible_families": sum(row["state"] == "zero_native_support" for row in rows),
                  "TP": counts["TP"], "FN": counts["FN"], "FP": counts["FP"], "TN": counts["TN"],
                  "interpretation_boundary": "agreement_with_retrospective_newer_assembly_proxy_not_physical_accuracy"}


def _run_srf(config: Mapping[str, Any], species: Mapping[str, Any], cell: Path, catalogue: dict[str, str],
             deadline: float, k: int) -> tuple[str, dict[str, float], dict[str, object]]:
    (cell / "tmp").mkdir()
    (cell / "reads.list").write_text("".join(f"{item['path']}\n" for item in species["reads"]))
    stages = []
    commands = srf_commands(config, species, cell, k)
    for name, command, output in commands:
        if name == "assemble" and (cell / "counts.txt").stat().st_size == 0:
            return "no_eligible_kmers", dict.fromkeys(catalogue, 0.0), {
                "stages": stages, "unassigned": [], "assigned_raw": dict.fromkeys(catalogue, 0)}
        if name == "elongate" and not read_fasta(cell / "srf.fa"):
            return "no_catalogue", dict.fromkeys(catalogue, 0.0), {
                "stages": stages, "unassigned": [], "assigned_raw": dict.fromkeys(catalogue, 0)}
        stages.append(execute_stage(config, cell, name, command, output, deadline))
    native_catalogue = read_fasta(cell / "srf.fa")
    native = parse_srf_abundance(cell / "srf.abundance.tsv", cell / "srf.bed", native_catalogue)
    correspondence = match_native_catalogue(native_catalogue, catalogue,
                                             config["correspondence"]["identity"],
                                             config["correspondence"]["max_integer_multiple"])
    save_json(cell / "correspondence.json", {key: asdict(value) for key, value in correspondence.items()})
    assigned_raw, unassigned = assign_srf_abundance(native, correspondence, catalogue)
    estimates = normalize_amounts(assigned_raw, species["total_library_bases"], species["genome_size_bp"])
    for row in unassigned:
        row["normalized_estimated_bp"] = row["native_retained_read_bp"] / species["total_library_bases"] * species["genome_size_bp"]
    return "ok", estimates, {"stages": stages, "unassigned": unassigned, "assigned_raw": assigned_raw,
                             "native_catalogue_count": len(native_catalogue),
                             "unassigned_native_read_bp": sum(int(row["native_retained_read_bp"]) for row in unassigned)}


def _run_mapping(config: Mapping[str, Any], species: Mapping[str, Any], cell: Path, catalogue: dict[str, str],
                 deadline: float) -> tuple[str, dict[str, float], dict[str, object]]:
    make_templates(catalogue, cell / "templates.fa", config["competitive_mapping"]["template_bp"])
    save_json(cell / "correspondence.json", {"state": "direct_frozen_TandemX_family_ids",
                                              "catalogue_records": len(catalogue)})
    command = mapping_command(config, species, cell)
    stage = execute_stage(config, cell, "map", command, cell / "mapping.paf", deadline)
    arrays, ambiguous = mapping_paf_to_unique_arrays(
        cell / "mapping.paf", catalogue,
        config["competitive_mapping"]["minimum_alignment_block_bp"],
        config["competitive_mapping"]["minimum_identity"])
    raw = dict.fromkeys(catalogue, 0)
    for record in arrays:
        raw[record.family_id] += record.end - record.start
    estimates = normalize_amounts(raw, species["total_library_bases"], species["genome_size_bp"])
    return "ok", estimates, {"stages": [stage], "ambiguous_native_read_bp": ambiguous,
                             "unassigned_native_read_bp": ambiguous, "assigned_raw": raw,
                             "unassigned": [{"native_id": "cross_family_overlap", "status": "ambiguous",
                                             "native_retained_read_bp": ambiguous,
                                             "normalized_estimated_bp": ambiguous / species["total_library_bases"] * species["genome_size_bp"],
                                             "candidate_matches": "multiple_TandemX_families"}] if ambiguous else []}


def run_cell(config: Mapping[str, Any], species: Mapping[str, Any], method: str, folder: Path,
             root: Path, deadline: float) -> dict[str, object]:
    folder.mkdir(parents=True, exist_ok=False)
    catalogue = load_tandemx_catalogue(_resolve(root, species["catalogue"]["path"]))
    reference = read_tsv(_resolve(root, species["reference_metrics"]["path"]))
    eligible = validate_reference(reference, species)
    started = time.perf_counter()
    status, estimates, detail = "technical_failure", None, {"stages": [], "unassigned": []}
    failure = None
    try:
        if method.startswith("srf_k"):
            status, estimates, detail = _run_srf(config, species, folder, catalogue, deadline, int(method[5:]))
        else:
            status, estimates, detail = _run_mapping(config, species, folder, catalogue, deadline)
    except (CellFailure, OSError, ValueError, subprocess.SubprocessError) as error:
        failure = f"{type(error).__name__}:{error}"
        status = "technical_failure"
    scored, summary = score_estimates(reference, estimates, method, failure or status)
    write_tsv(folder / "family_prioritization.tsv", scored, list(PRIORITIZATION_FIELDS))
    assigned_raw = detail.get("assigned_raw", {})
    estimate_rows = ([{"family_id": family, "native_retained_read_bp": assigned_raw[family],
                       "read_estimated_bp": value} for family, value in estimates.items()]
                     if estimates is not None else [])
    write_tsv(folder / "family_estimates.tsv", estimate_rows,
              ["family_id", "native_retained_read_bp", "read_estimated_bp"])
    unassigned = detail.get("unassigned", [])
    write_tsv(folder / "unassigned_abundance.tsv", unassigned,
              ["native_id", "status", "native_retained_read_bp", "normalized_estimated_bp", "candidate_matches"])
    stages = detail.get("stages", [])
    result = {"species": species["id"], "method": method, "status": status, "failure": failure,
              **summary, "catalogue_records": len(catalogue), "eligible_records": len(eligible),
              "native_catalogue_count": detail.get("native_catalogue_count", "N/A_shared_catalogue"),
              "ambiguous_native_read_bp": detail.get("ambiguous_native_read_bp", 0),
              "unassigned_native_read_bp": detail.get("unassigned_native_read_bp", 0),
              "runtime_seconds": time.perf_counter() - started,
              "stage_runtime_seconds": sum(float(row["runtime_seconds"]) for row in stages),
              "peak_stage_rss_mib": max((float(row["peak_rss_mib"]) for row in stages), default=0),
              "cell_bytes": tree_bytes(folder), "completed_utc": utc_now(),
              "native_output_reuse": False}
    save_json(folder / "result.json", result)
    return result


def snapshot_sources(root: Path, outdir: Path) -> dict[str, str]:
    snapshot = outdir / "source_snapshot"
    hashes = {}
    for name in SNAPSHOT_FILES:
        source, destination = root / name, snapshot / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        hashes[name] = digest_file(source)
        if digest_file(destination) != hashes[name]:
            raise ValueError(f"Source changed while snapshotting: {name}")
    return hashes


def execute(config_path: Path, config: dict[str, Any], outdir: Path, root: Path) -> None:
    manifest = validate_config(config, root, verify_files=True)
    if outdir.exists():
        raise FileExistsError(f"Formal output root must not exist: {outdir}")
    launch_free, launch_free_method = free_bytes(outdir.parent)
    if launch_free < config["resources"]["minimum_free_bytes"]:
        raise CellFailure("minimum_free_bytes_before_launch")
    outdir.mkdir(parents=True)
    shutil.copyfile(config_path, outdir / "run_config.json")
    config_hash = digest_file(config_path)
    source_hashes = snapshot_sources(root, outdir)
    manifest.update({"config_sha256": config_hash, "protocol_sha256": config["protocol"]["sha256"],
                     "validated_utc": utc_now()})
    save_json(outdir / "input_manifest.json", manifest)
    git = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=root,
                         text=True, capture_output=True, check=False)
    save_json(outdir / "environment.json", {"python": sys.version, "platform": platform.platform(),
              "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
              "git_status": git.stdout, "config_sha256": config_hash, "source_sha256": source_hashes,
              "tools": config["tools"], "threads": 1, "started_utc": utc_now(),
              "memory_method": "max_direct_child_wait4_RSS_per_sequential_stage_controller_excluded",
              "disk_method": "post_stage_cell_bytes_and_free_bytes_not_continuous_peak",
              "launch_free_bytes": launch_free, "launch_free_bytes_method": launch_free_method})
    plan = plan_rows(config)
    write_tsv(outdir / "plan.tsv", plan, list(plan[0]))
    deadline = time.perf_counter() + config["resources"]["global_launch_seconds"]
    summaries = []
    stop_reason = None
    for item in plan:
        species = next(row for row in config["species"] if row["id"] == item["species"])
        if stop_reason:
            summaries.append({**item, "status": "not_run_prior_failure", "failure": stop_reason})
            continue
        folder = outdir / species["id"] / str(item["method"])
        result = run_cell(config, species, str(item["method"]), folder, root, deadline)
        summaries.append({**item, **result})
        if result["status"] == "technical_failure":
            stop_reason = str(result["failure"])
    fields = list(dict.fromkeys(key for row in summaries for key in row))
    write_tsv(outdir / "summary.tsv", [{key: row.get(key, "N/A") for key in fields} for row in summaries], fields)
    complete = stop_reason is None and len(summaries) == 6
    save_json(outdir / "completion.json", {"complete": complete, "cells_planned": 6,
              "cells_completed": sum(row.get("status") not in {"not_run_prior_failure", "technical_failure"} for row in summaries),
              "failed_cells": sum(row.get("status") == "technical_failure" for row in summaries),
              "not_run_cells": sum(row.get("status") == "not_run_prior_failure" for row in summaries),
              "stop_reason": stop_reason, "completed_utc": utc_now(),
              "interpretation_boundary": "retrospective_proxy_agreement_not_physical_accuracy"})
    files = sorted(path for path in outdir.rglob("*") if path.is_file() and path.name != "archive_manifest.json")
    save_json(outdir / "archive_manifest.json", {"files": [{"path": str(path.relative_to(outdir)),
              "bytes": path.stat().st_size, "sha256": digest_file(path)} for path in files]})
    if not complete:
        raise RuntimeError(f"Formal comparator incomplete: {stop_reason}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true", help="Validate frozen structure only; do not hash large inputs or run tools")
    mode.add_argument("--audit-bindings", action="store_true", help="Recompute every input/tool/protocol hash; do not run tools")
    mode.add_argument("--execute", action="store_true", help="Run all six native cells after review")
    parser.add_argument("--ack")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    config = load_config(args.config)
    if args.plan:
        validate_config(config, root, verify_files=False)
        print(json.dumps({"status": "structurally_valid_unexecuted", "cells": plan_rows(config)}, indent=2))
        return
    if args.audit_bindings:
        print(json.dumps(validate_config(config, root, verify_files=True), indent=2))
        return
    if args.ack != ACK:
        raise ValueError(f"Formal execution requires --ack {ACK}")
    if args.outdir is None:
        raise ValueError("--outdir is required with --execute")
    execute(args.config.resolve(), config, args.outdir.resolve(), root)


if __name__ == "__main__":
    main()
