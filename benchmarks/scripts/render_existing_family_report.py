#!/usr/bin/env python3
"""Render a TandemX family report from compact, already-produced outputs.

This utility copies catalogue/abundance/array evidence only.  It neither opens
raw reads nor reruns discovery, quantification, or localization.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from tandemx.compare.mvp import (
    DEFAULT_COLLAPSE_THRESHOLD,
    DEFAULT_OVEREXPANSION_THRESHOLD,
    CompareConfig,
    compare_toy_abundance,
)
from tandemx.report_html import write_html_report


REQUIRED_DISCOVERY = ("families.tsv", "monomers.fa", "candidate_reads.tsv")
OPTIONAL_DISCOVERY = ("monomer_membership.tsv", "family_hierarchy.tsv", "family_audit_summary.json", "discovery_summary.json")
RECOVERY_LOCKED = ("recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy(source: Path, destination: Path, manifest: list[dict[str, object]]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    source_hash, destination_hash = sha256(source), sha256(destination)
    if source_hash != destination_hash:
        raise RuntimeError(f"Copied file hash mismatch: {source}")
    manifest.append({"source_path": str(source.resolve()), "destination_path": str(destination),
                     "source_sha256": source_hash, "destination_sha256": destination_hash,
                     "bytes": source.stat().st_size})


def _copy_recovery(source: Path, destination: Path, manifest: list[dict[str, object]]) -> bool:
    lock_path = source / "candidate_lock.json"
    if not lock_path.is_file():
        return False
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("status") != "candidate_generation_complete":
        raise ValueError(f"Recovery lock is not complete: {lock_path}")
    outputs = lock.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError(f"Recovery lock lacks output hashes: {lock_path}")
    for name in RECOVERY_LOCKED:
        path, expected = source / name, outputs.get(name)
        if not isinstance(expected, str) or not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Recovery lock hash mismatch for {path}")
    for name in (*RECOVERY_LOCKED, "candidate_lock.json", "recovery_report.html", "flank_audit.tsv",
                 "recovery_assessment.json", "recovery_validation_report.html"):
        path = source / name
        if path.is_file():
            _copy(path, destination / name, manifest)
    return True


def render_existing_family_report(*, discover: Path, copy_number: Path, arrays: Path | None,
                                  recovery: Path | None, outdir: Path) -> Path:
    """Create a new self-contained report directory from hash-recorded products."""
    discover, copy_number, outdir = Path(discover), Path(copy_number), Path(outdir)
    if outdir.exists():
        raise ValueError(f"Output directory must not already exist: {outdir}")
    missing = [str(discover / name) for name in REQUIRED_DISCOVERY if not (discover / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing required discovery products: " + ", ".join(missing))
    if not copy_number.is_file():
        raise FileNotFoundError(f"Missing copy-number table: {copy_number}")
    if arrays is not None and not Path(arrays).is_file():
        raise FileNotFoundError(f"Missing arrays BED: {arrays}")
    outdir.mkdir(parents=True)
    manifest: list[dict[str, object]] = []
    for name in REQUIRED_DISCOVERY + OPTIONAL_DISCOVERY:
        if (discover / name).is_file():
            _copy(discover / name, outdir / "discover" / name, manifest)
    _copy(copy_number, outdir / "quantify" / "copy_number.tsv", manifest)
    quantify_root = copy_number.parent.parent
    for name in ("run_config.yaml",):
        if (copy_number.parent / name).is_file():
            _copy(copy_number.parent / name, outdir / "quantify" / name, manifest)
    if (quantify_root / "automatic_defaults.json").is_file():
        _copy(quantify_root / "automatic_defaults.json", outdir / "automatic_defaults.json", manifest)
    if arrays is not None:
        _copy(Path(arrays), outdir / "locate" / "arrays.bed", manifest)
        compare_toy_abundance(CompareConfig(copy_number=outdir / "quantify" / "copy_number.tsv",
                                            arrays=outdir / "locate" / "arrays.bed", outdir=outdir / "compare",
                                            collapse_threshold=DEFAULT_COLLAPSE_THRESHOLD,
                                            overexpansion_threshold=DEFAULT_OVEREXPANSION_THRESHOLD))
        (outdir / "comparison_reuse_config.json").write_text(json.dumps({
            "command": "existing CompareConfig.compare_toy_abundance",
            "collapse_threshold": DEFAULT_COLLAPSE_THRESHOLD,
            "overexpansion_threshold": DEFAULT_OVEREXPANSION_THRESHOLD,
            "interpretation": "recomputed from frozen copied copy-number and array tables; not a new benchmark",
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    recovery_copied = _copy_recovery(Path(recovery), outdir / "recovery", manifest) if recovery else False
    records = [
        {"step": "discover", "exit_status": 0, "output_validated": True, "notes": "source_reuse_hash_verified"},
        {"step": "quantify", "exit_status": 0, "output_validated": True, "notes": "source_reuse_hash_verified"},
    ]
    if arrays is not None:
        records.extend({"step": step, "exit_status": 0, "output_validated": True, "notes": "source_reuse_hash_verified"} for step in ("locate", "compare"))
    if recovery_copied:
        records.append({"step": "recovery", "exit_status": 0, "output_validated": True, "notes": "source_reuse_hash_verified"})
    (outdir / "source_reuse_manifest.json").write_text(json.dumps({
        "schema_version": 1, "source_reuse_only": True, "raw_reads_accessed": False,
        "discovery_rerun": False, "quantification_rerun": False, "localization_rerun": False,
        "compare_recomputed_from_frozen_tables": arrays is not None, "recovery_copied": recovery_copied,
        "files": manifest, "records": records,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_html_report(outdir, records)
    return outdir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discover", required=True, type=Path)
    parser.add_argument("--copy-number", required=True, type=Path)
    parser.add_argument("--arrays", type=Path)
    parser.add_argument("--recovery", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    render_existing_family_report(discover=args.discover, copy_number=args.copy_number,
                                  arrays=args.arrays, recovery=args.recovery, outdir=args.outdir)


if __name__ == "__main__":
    main()
