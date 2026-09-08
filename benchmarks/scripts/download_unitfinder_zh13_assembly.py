#!/usr/bin/env python3
"""Download and stream-QC the frozen official ZH13-T2T assembly input."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any, BinaryIO
from urllib.request import urlopen


def file_hashes(path: Path) -> dict[str, Any]:
    md5 = hashlib.md5()  # noqa: S324 - required to verify the publisher checksum
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            size += len(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    return {"bytes": size, "md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def stream_download(source: BinaryIO, target: Path) -> int:
    size = 0
    with target.open("xb") as output:
        while chunk := source.read(4 * 1024 * 1024):
            output.write(chunk)
            size += len(chunk)
    return size


def fasta_qc(path: Path) -> dict[str, Any]:
    records = bases = acgt_bases = other_bases = 0
    current_bases = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ValueError(f"empty FASTA header at line {line_number}")
                if records and current_bases == 0:
                    raise ValueError(f"empty FASTA record before line {line_number}")
                records += 1
                current_bases = 0
                continue
            if records == 0:
                raise ValueError(f"sequence before FASTA header at line {line_number}")
            upper = line.upper()
            line_acgt = sum(upper.count(base) for base in "ACGT")
            bases += len(line)
            acgt_bases += line_acgt
            other_bases += len(line) - line_acgt
            current_bases += len(line)
    if records == 0 or current_bases == 0 or bases == 0:
        raise ValueError("empty or incomplete FASTA")
    return {
        "record_count": records,
        "base_count": bases,
        "acgt_bases": acgt_bases,
        "other_bases": other_bases,
        "other_fraction": other_bases / bases,
    }


def write_receipt(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def download(config_path: Path, outdir: Path, timeout: float = 60.0) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"choose a new output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    assembly = config["assembly"]
    outdir.mkdir(parents=True)
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    shutil.copyfile(config_path, snapshot / config_path.name)
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    receipt_path = outdir / "run_receipt.json"
    receipt: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "complete": False,
        "fate": "download_in_progress",
        "assembly": {
            "accession": assembly["accession"],
            "filename": assembly["filename"],
            "url": assembly["url"],
        },
        "boundary": config["boundary"],
    }
    write_receipt(receipt_path, receipt)
    partial = outdir / f"{assembly['filename']}.partial"
    final = outdir / assembly["filename"]
    try:
        with urlopen(assembly["url"], timeout=timeout) as source:
            response = {
                "status": getattr(source, "status", None),
                "content_length": source.headers.get("Content-Length"),
                "last_modified": source.headers.get("Last-Modified"),
                "content_type": source.headers.get("Content-Type"),
            }
            transferred = stream_download(source, partial)
        hashes = file_hashes(partial)
        if transferred != assembly["expected_bytes"] or hashes["bytes"] != assembly["expected_bytes"]:
            raise ValueError("downloaded compressed byte count differs from frozen source")
        if hashes["md5"] != assembly["expected_md5"]:
            raise ValueError("downloaded MD5 differs from frozen publisher checksum")
        qc = fasta_qc(partial)
        partial.replace(final)
        receipt.update(
            {
                "complete": True,
                "fate": "source_enrollment_passed",
                "response": response,
                "compressed_file": {"file": final.name, **hashes},
                "fasta_qc": qc,
            }
        )
    except Exception as error:
        receipt.update(
            {
                "fate": "source_enrollment_failure",
                "error_type": type(error).__name__,
                "error": str(error),
                "partial_file": (
                    {"file": partial.name, **file_hashes(partial)}
                    if partial.is_file()
                    else None
                ),
            }
        )
    write_receipt(receipt_path, receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()
    receipt = download(args.config, args.outdir, args.timeout)
    return 0 if receipt["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
