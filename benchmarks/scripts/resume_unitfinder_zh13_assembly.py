#!/usr/bin/env python3
"""Resume the frozen truncated ZH13 source transfer without changing its evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any
from urllib.request import Request, urlopen

from benchmarks.scripts.download_unitfinder_zh13_assembly import fasta_qc, file_hashes


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate_parent(config: dict[str, Any]) -> tuple[Path, Path]:
    record = config["parent_failure"]
    parent = Path(record["directory"]).resolve()
    for name, expected in record["artifacts"].items():
        path = parent / name
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"ZH13 parent failure artifact changed: {name}")
    receipt = json.loads((parent / "run_receipt.json").read_text(encoding="utf-8"))
    diagnosis = json.loads(
        (parent / "transport_diagnosis.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("complete") is not False
        or receipt.get("fate") != record["fate"]
        or diagnosis.get("classification") != record["classification"]
        or diagnosis.get("safe_range_resume_offset") != record["resume_offset"]
    ):
        raise ValueError("ZH13 parent failure receipt or diagnosis changed")
    partial = parent / receipt["partial_file"]["file"]
    if partial.stat().st_size != record["resume_offset"]:
        raise ValueError("ZH13 parent partial size differs from resume offset")
    return parent, partial


def write_receipt(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def resume(config_path: Path, outdir: Path, timeout: float = 60.0) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"choose a new output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    parent, parent_partial = validate_parent(config)
    assembly = config["assembly"]
    offset = config["parent_failure"]["resume_offset"]
    outdir.mkdir(parents=True)
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    shutil.copyfile(config_path, snapshot / config_path.name)
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    partial = outdir / f"{assembly['filename']}.partial"
    final = outdir / assembly["filename"]
    receipt_path = outdir / "run_receipt.json"
    receipt: dict[str, Any] = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "complete": False,
        "fate": "range_resume_in_progress",
        "parent_failure_directory": str(parent),
        "resume_offset": offset,
        "boundary": config["boundary"],
    }
    write_receipt(receipt_path, receipt)
    try:
        shutil.copyfile(parent_partial, partial)
        if partial.stat().st_size != offset or digest(partial) != digest(parent_partial):
            raise ValueError("copied parent partial differs before range request")
        request = Request(assembly["url"], headers={"Range": f"bytes={offset}-"})
        with urlopen(request, timeout=timeout) as source:
            response = {
                "status": source.status,
                "content_range": source.headers.get("Content-Range"),
                "content_length": source.headers.get("Content-Length"),
                "last_modified": source.headers.get("Last-Modified"),
                "etag": source.headers.get("ETag"),
                "accept_ranges": source.headers.get("Accept-Ranges"),
            }
            if source.status != config["range_acceptance"]["status"]:
                raise ValueError("range server did not return frozen HTTP 206")
            if response["content_range"] != config["range_acceptance"]["content_range"]:
                raise ValueError("range server Content-Range differs from frozen response")
            appended = 0
            with partial.open("ab") as output:
                while chunk := source.read(4 * 1024 * 1024):
                    output.write(chunk)
                    appended += len(chunk)
        if appended != config["range_acceptance"]["remaining_bytes"]:
            raise ValueError("range response ended before the frozen remaining byte count")
        hashes = file_hashes(partial)
        if hashes["bytes"] != assembly["expected_bytes"]:
            raise ValueError("resumed compressed byte count differs from frozen source")
        if hashes["md5"] != assembly["expected_md5"]:
            raise ValueError("resumed MD5 differs from frozen publisher checksum")
        qc = fasta_qc(partial)
        partial.replace(final)
        receipt.update(
            {
                "complete": True,
                "fate": "source_enrollment_passed",
                "response": response,
                "appended_bytes": appended,
                "compressed_file": {"file": final.name, **hashes},
                "fasta_qc": qc,
            }
        )
    except Exception as error:
        receipt.update(
            {
                "fate": "range_resume_failure",
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
    receipt = resume(args.config, args.outdir, args.timeout)
    return 0 if receipt["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
