#!/usr/bin/env python3
"""Archive compact unitFinder source-enrollment evidence without genome payloads."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

from benchmarks.scripts.download_unitfinder_zh13_assembly import file_hashes


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate(
    source: Path, config_path: Path
) -> tuple[dict[str, Any], dict[str, Any], Path, dict[str, Any]]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt_path = source / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("complete") is True and receipt.get("fate") == "source_enrollment_passed":
        record = receipt["compressed_file"]
        payload = source / record["file"]
        observed = file_hashes(payload)
        if observed != {key: record[key] for key in ("bytes", "md5", "sha256")}:
            raise ValueError("enrolled source payload changed")
        if (
            observed["bytes"] != config["assembly"]["expected_bytes"]
            or observed["md5"] != config["assembly"]["expected_md5"]
            or receipt.get("fasta_qc", {}).get("record_count", 0) <= 0
        ):
            raise ValueError("enrolled source does not satisfy the frozen source contract")
    elif receipt.get("complete") is False and receipt.get("fate") in {
        "source_enrollment_failure",
        "range_resume_failure",
    }:
        record = receipt["partial_file"]
        payload = source / record["file"]
        observed = file_hashes(payload)
        if observed != {key: record[key] for key in ("bytes", "md5", "sha256")}:
            raise ValueError("retained partial source changed")
        diagnosis_path = source / "transport_diagnosis.json"
        if diagnosis_path.is_file():
            diagnosis = json.loads(diagnosis_path.read_text(encoding="utf-8"))
            if diagnosis.get("run_receipt_sha256") != digest(receipt_path):
                raise ValueError("source failure diagnosis differs from receipt")
    else:
        raise ValueError("source enrollment receipt lacks a recognized explicit fate")
    return config, receipt, payload, observed


def archive(source: Path, config_path: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite archive: {outdir}")
    source = source.resolve()
    config, receipt, payload, payload_hashes = validate(source, config_path)
    outdir.mkdir(parents=True)
    copied: list[dict[str, Any]] = []

    def copy(path: Path, target: Path) -> None:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        if path.stat().st_size != target.stat().st_size or digest(path) != digest(target):
            raise OSError(f"unitFinder source archive copy differs: {path}")
        copied.append(
            {
                "file": target.relative_to(outdir).as_posix(),
                "source": str(path),
                "bytes": target.stat().st_size,
                "sha256": digest(target),
            }
        )

    copy(config_path, outdir / "preregistration.json")
    for name in ("run_receipt.json", "transport_diagnosis.json"):
        path = source / name
        if path.is_file():
            copy(path, outdir / "results" / name)
    snapshot = source / "source_snapshot"
    if snapshot.is_dir():
        for path in sorted(item for item in snapshot.rglob("*") if item.is_file()):
            copy(path, outdir / "results/source_snapshot" / path.relative_to(snapshot))
    for name in ("independent_verification.json", "sequence_lengths.tsv"):
        path = source / name
        if path.is_file():
            copy(path, outdir / "results" / name)

    external = {
        "path": str(payload),
        "retained_outside_git": True,
        **payload_hashes,
    }
    headline = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": receipt["complete"],
        "experiment_fate": receipt["fate"],
        "accession": config["assembly"]["accession"],
        "external_payload": external,
        "fasta_qc": receipt.get("fasta_qc"),
        "independent_verification": (
            json.loads(
                (source / "independent_verification.json").read_text(encoding="utf-8")
            )
            if (source / "independent_verification.json").is_file()
            else None
        ),
        "boundary": config["boundary"],
    }
    headline_path = outdir / "headline_summary.json"
    headline_path.write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if receipt["complete"]:
        message = (
            "The exact official assembly file passed frozen compressed-size, publisher "
            "MD5 and streaming FASTA checks. The large payload remains on T7 and is "
            "identified here by path and hashes. This is input identity, not unitFinder "
            "reproducibility or biological accuracy."
        )
    else:
        message = (
            "The frozen official-assembly acquisition failed and its partial payload is "
            "retained on T7 by exact path and hashes. This is a source/transport failure, "
            "not a unitFinder execution or accuracy result."
        )
    readme = outdir / "README.md"
    readme.write_text(f"# {config['experiment_id']}\n\n{message}\n", encoding="utf-8")
    for path in (headline_path, readme):
        copied.append(
            {
                "file": path.name,
                "source": "generated_by_archive_unitfinder_source_enrollment",
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )
    manifest = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": receipt["complete"],
        "file_count": len(copied),
        "files": copied,
    }
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
