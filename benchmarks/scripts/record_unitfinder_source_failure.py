#!/usr/bin/env python3
"""Diagnose and retain a failed frozen unitFinder source enrollment."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Callable
from urllib.request import Request, urlopen


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def head_metadata(url: str, timeout: float) -> dict[str, Any]:
    request = Request(url, method="HEAD")
    with urlopen(request, timeout=timeout) as response:
        return {
            "status": response.status,
            "content_length": response.headers.get("Content-Length"),
            "last_modified": response.headers.get("Last-Modified"),
            "etag": response.headers.get("ETag"),
            "accept_ranges": response.headers.get("Accept-Ranges"),
            "content_type": response.headers.get("Content-Type"),
        }


def classify(
    partial_bytes: int,
    expected_bytes: int,
    gzip_exit_code: int,
    head_content_length: str | None,
) -> str:
    if (
        partial_bytes < expected_bytes
        and gzip_exit_code != 0
        and head_content_length is not None
        and int(head_content_length) == expected_bytes
    ):
        return "truncated_transfer_before_expected_content_length"
    return "unresolved_source_enrollment_failure"


def diagnose(
    result_dir: Path,
    config_path: Path,
    timeout: float = 60.0,
    head_reader: Callable[[str, float], dict[str, Any]] = head_metadata,
) -> dict[str, Any]:
    output = result_dir / "transport_diagnosis.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite source diagnosis: {output}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt_path = result_dir / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("complete") is not False or receipt.get("fate") != "source_enrollment_failure":
        raise ValueError("source receipt is not an explicit failed enrollment")
    partial_record = receipt.get("partial_file")
    if not isinstance(partial_record, dict):
        raise ValueError("failed source receipt lacks a retained partial file")
    partial = result_dir / partial_record["file"]
    if (
        not partial.is_file()
        or partial.stat().st_size != partial_record["bytes"]
        or digest(partial) != partial_record["sha256"]
    ):
        raise ValueError("retained partial source changed")
    gzip_check = subprocess.run(
        ["gzip", "-t", str(partial)], capture_output=True, text=True, check=False
    )
    head = head_reader(config["assembly"]["url"], timeout)
    classification = classify(
        partial.stat().st_size,
        config["assembly"]["expected_bytes"],
        gzip_check.returncode,
        head.get("content_length"),
    )
    result = {
        "schema_version": 1,
        "complete": True,
        "classification": classification,
        "run_receipt_sha256": digest(receipt_path),
        "config_sha256": digest(config_path),
        "partial_file": partial_record,
        "expected_bytes": config["assembly"]["expected_bytes"],
        "expected_md5": config["assembly"]["expected_md5"],
        "head_response": head,
        "gzip_test": {
            "command": ["gzip", "-t", str(partial)],
            "exit_code": gzip_check.returncode,
            "stderr": gzip_check.stderr.strip(),
        },
        "safe_range_resume_offset": (
            partial.stat().st_size
            if classification == "truncated_transfer_before_expected_content_length"
            else None
        ),
        "boundary": (
            "This classifies transport completeness only. It does not validate the "
            "assembly, unitFinder output, Table S1 agreement or biological accuracy."
        ),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()
    diagnose(args.result_dir, args.config, args.timeout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
