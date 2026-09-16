"""Validate preregistered benchmark inputs without changing source files.

Usage: python benchmarks/scripts/validate_input_manifest.py MANIFEST.json
The JSON report is written to stdout; exit 1 means an input is ineligible.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, BinaryIO


SPLITS = {"development", "validation", "final-heldout"}
TRUTH_TYPES = {"none", "simulated-exact", "edited-delta", "physical", "assembly-proxy"}
STATUSES = {"eligible", "pending", "invalid", "storage-failure", "source-unresolved"}
FILE_KINDS = {"fasta", "fastq"}
FILE_ROLES = {"reads", "assembly", "catalogue"}
TASK_ROLES = {
    "input-integrity-only": set(),
    "raw-read-discovery": {"reads"},
    "assembly-structure": {"assembly"},
    "fixed-catalogue-abundance": {"reads", "catalogue"},
    "read-assembly-deficit": {"reads", "assembly"},
    "read-assembly-structure": {"reads", "assembly"},
}
PAIRING_STATES = {"verified", "unresolved", "not-applicable"}
PROVENANCE_STATES = {"verified", "unresolved"}
REQUIRED = {
    "dataset", "sample", "species", "truth_type", "family", "array_coordinates",
    "coverage", "read_source", "assembly_source", "split", "allowed_tuning",
    "metrics", "competitors", "software_versions", "donor_id", "accessions",
    "assembly_lineage", "family_homology_group", "status", "files",
    "benchmark_task", "source_pairing_status", "source_pairing_evidence",
    "source_provenance_status", "source_provenance_evidence",
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
BASES = frozenset(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn")


def _require_text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label}: nonempty string required")
    return value.strip()


def _open(path: Path) -> BinaryIO:
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    if path.suffix == ".gz":
        raise ValueError(f"{path}: .gz suffix without gzip magic")
    return path.open("rb")


def _records(handle: BinaryIO, kind: str) -> tuple[int, int]:
    count = bases = 0
    if kind == "fastq":
        while header := handle.readline():
            sequence = handle.readline().rstrip(b"\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip(b"\r\n")
            if not header.startswith(b"@") or not sequence or not plus.startswith(b"+"):
                raise ValueError(f"malformed FASTQ record {count + 1}")
            if len(sequence) != len(quality) or not set(sequence) <= BASES:
                raise ValueError(f"invalid FASTQ bases/quality at record {count + 1}")
            count += 1
            bases += len(sequence)
    else:
        seen = False
        current_len = 0
        for line in handle:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if len(line) == 1 or (seen and current_len == 0):
                    raise ValueError("empty FASTA header or sequence")
                count += 1
                seen = True
                current_len = 0
            else:
                if not seen or not line or not set(line) <= BASES:
                    raise ValueError(f"invalid FASTA sequence at record {count}")
                current_len += len(line)
                bases += len(line)
        if seen and current_len == 0:
            raise ValueError("empty final FASTA sequence")
    if count == 0:
        raise ValueError("empty sequence file")
    return count, bases


def check_file(spec: dict[str, Any]) -> dict[str, Any]:
    path = Path(_require_text(spec.get("path"), "file.path"))
    kind = spec.get("kind")
    role = spec.get("role")
    expected_sha = spec.get("sha256")
    expected_size = spec.get("bytes")
    if kind not in FILE_KINDS or not isinstance(expected_sha, str) or not HEX64.fullmatch(expected_sha):
        raise ValueError(f"{path}: kind or sha256 invalid")
    if role not in FILE_ROLES or (role in {"assembly", "catalogue"} and kind != "fasta"):
        raise ValueError(f"{path}: invalid file role/kind")
    if not path.is_absolute():
        raise ValueError(f"{path}: absolute path required")
    if not isinstance(expected_size, int) or expected_size <= 0:
        raise ValueError(f"{path}: positive byte size required")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    actual_size = path.stat().st_size
    if actual_size != expected_size or digest.hexdigest() != expected_sha:
        raise ValueError(f"{path}: size/SHA-256 mismatch")
    with _open(path) as handle:
        records, bases = _records(handle, kind)
    return {"path": str(path.resolve()), "kind": kind, "role": role, "bytes": actual_size,
            "sha256": digest.hexdigest(), "records": records, "bases": bases,
            "structure": "valid"}


def validate(manifest: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("schema_version") != 2 or not isinstance(manifest.get("datasets"), list):
        raise ValueError("schema_version 2 and datasets list required")
    errors: list[str] = []
    checked: list[dict[str, Any]] = []
    groups: dict[tuple[str, str], set[str]] = defaultdict(set)
    ids: set[str] = set()
    for entry in manifest["datasets"]:
        name = entry.get("dataset", "<unnamed>")
        try:
            missing = REQUIRED - entry.keys()
            if missing:
                raise ValueError(f"missing fields: {', '.join(sorted(missing))}")
            name = _require_text(name, "dataset")
            if name in ids:
                raise ValueError("duplicate dataset")
            ids.add(name)
            if entry["split"] not in SPLITS or entry["truth_type"] not in TRUTH_TYPES:
                raise ValueError("invalid split or truth_type")
            if entry["status"] not in STATUSES:
                raise ValueError("invalid status")
            task = entry["benchmark_task"]
            if task not in TASK_ROLES:
                raise ValueError("invalid benchmark_task")
            if entry["source_pairing_status"] not in PAIRING_STATES:
                raise ValueError("invalid source_pairing_status")
            if entry["source_provenance_status"] not in PROVENANCE_STATES:
                raise ValueError("invalid source_provenance_status")
            if not isinstance(entry["source_pairing_evidence"], list):
                raise ValueError("source_pairing_evidence must be list")
            if not isinstance(entry["source_provenance_evidence"], list):
                raise ValueError("source_provenance_evidence must be list")
            for field in ("source_pairing_evidence", "source_provenance_evidence"):
                for evidence in entry[field]:
                    _require_text(evidence, field)
            for field in ("sample", "species", "donor_id", "family_homology_group", "assembly_lineage", "read_source", "assembly_source"):
                _require_text(entry[field], field)
            if entry["coverage"] is not None and (not isinstance(entry["coverage"], (int, float)) or entry["coverage"] < 0):
                raise ValueError("coverage must be nonnegative numeric or null")
            if not isinstance(entry["family"], list) or not isinstance(entry["array_coordinates"], list):
                raise ValueError("family and array_coordinates must be lists")
            if not isinstance(entry["accessions"], list) or not entry["accessions"]:
                raise ValueError("accessions must be nonempty list")
            for accession in entry["accessions"]:
                _require_text(accession, "accession")
            if not isinstance(entry["allowed_tuning"], bool):
                raise ValueError("allowed_tuning must be boolean")
            if entry["split"] != "development" and entry["allowed_tuning"]:
                raise ValueError("tuning is permitted only in development")
            if not all(isinstance(entry[f], list) for f in ("metrics", "competitors", "files")):
                raise ValueError("metrics, competitors and files must be lists")
            if not isinstance(entry["software_versions"], dict):
                raise ValueError("software_versions must be object")
            for key in ("species", "donor_id", "sample", "assembly_lineage", "family_homology_group"):
                groups[(key, entry[key])].add(entry["split"])
            for accession in entry["accessions"]:
                groups[("accession", accession)].add(entry["split"])
            files = [check_file(file) for file in entry["files"]]
            for file in files:
                groups[("file_path", file["path"])].add(entry["split"])
                groups[("file_sha256", file["sha256"])].add(entry["split"])
            if entry["status"] == "eligible":
                if task == "input-integrity-only":
                    raise ValueError("integrity-only inventory cannot be benchmark-eligible")
                missing_roles = TASK_ROLES[task] - {file["role"] for file in files}
                if missing_roles:
                    raise ValueError(f"eligible {task} missing file roles: {sorted(missing_roles)}")
                if entry["source_provenance_status"] != "verified" or not entry["source_provenance_evidence"]:
                    raise ValueError("eligible dataset needs verified source provenance and evidence")
                if {"reads", "assembly"} <= TASK_ROLES[task]:
                    if entry["source_pairing_status"] != "verified" or not entry["source_pairing_evidence"]:
                        raise ValueError("eligible read+assembly task needs verified pairing and evidence")
                elif entry["source_pairing_status"] == "unresolved":
                    raise ValueError("eligible dataset cannot have unresolved source pairing")
            checked.append({"dataset": name, "status": entry["status"], "files": files})
        except (ValueError, OSError, EOFError) as exc:
            errors.append(f"{name}: {exc}")
    for (key, value), splits in sorted(groups.items()):
        if len(splits) > 1:
            errors.append(f"split leakage: {key}={value} in {sorted(splits)}")
    return {"valid": not errors, "datasets_checked": checked, "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        result = validate(json.loads(args.manifest.read_text()))
    except (ValueError, OSError) as exc:
        result = {"valid": False, "datasets_checked": [], "errors": [str(exc)]}
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
