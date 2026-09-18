"""Bounded readback and disposable write/readback probe for a mounted T7."""

from __future__ import annotations

import argparse
import errno
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Any


T7_MOUNT = Path("/Volumes/T7")
CHUNK_BYTES = 4 * 1024 * 1024


def mount_identity() -> tuple[int, bool]:
    if not T7_MOUNT.is_mount():
        raise RuntimeError("T7 mount is unavailable")
    return T7_MOUNT.stat().st_dev, True


def stream_sha256(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    count = 0
    with path.open("rb", buffering=0) as handle:
        if hasattr(fcntl, "F_NOCACHE"):
            fcntl.fcntl(handle.fileno(), fcntl.F_NOCACHE, 1)
        while block := handle.read(CHUNK_BYTES):
            digest.update(block)
            count += len(block)
    return digest.hexdigest(), count


def check_inputs(manifest: Path) -> dict[str, Any]:
    expected = json.loads(manifest.read_text())
    device, _ = mount_identity()
    results: list[dict[str, Any]] = []
    for item in expected["inputs"]:
        path = Path(item["path"])
        if not path.is_relative_to(T7_MOUNT):
            raise ValueError(f"input outside T7: {path}")
        began = time.monotonic()
        actual_size = path.stat().st_size
        required_size = item["expected_bytes"]
        if required_size is not None and actual_size != required_size:
            raise ValueError(f"{item['name']}: size mismatch {actual_size} != {required_size}")
        actual_hash, read_bytes = stream_sha256(path)
        if actual_hash != item["expected_sha256"] or read_bytes != actual_size:
            raise ValueError(f"{item['name']}: hash/read-size mismatch")
        expanded_bytes = None
        if path.suffix == ".gz":
            expanded_bytes = 0
            with gzip.open(path, "rb") as handle:
                while block := handle.read(CHUNK_BYTES):
                    expanded_bytes += len(block)
        if mount_identity()[0] != device:
            raise RuntimeError("T7 mount device changed during input audit")
        results.append({
            "name": item["name"],
            "path": str(path),
            "bytes": actual_size,
            "sha256": actual_hash,
            "gzip_complete": True if expanded_bytes is not None else None,
            "decompressed_bytes": expanded_bytes,
            "elapsed_seconds": round(time.monotonic() - began, 3),
            "status": "pass",
        })
        print(f"PASS {item['name']} {actual_size} bytes", flush=True)
    return {"mode": "inputs", "mount_device": device, "results": results, "status": "pass"}


def write_probe(directory: Path, size_gib: int) -> dict[str, Any]:
    if size_gib < 2 or size_gib > 8:
        raise ValueError("write probe must be between 2 and 8 GiB")
    if not directory.is_relative_to(T7_MOUNT):
        raise ValueError("probe directory outside T7")
    device, _ = mount_identity()
    if not directory.is_dir():
        raise FileNotFoundError(directory)
    requested = size_gib * 1024**3
    free = shutil.disk_usage(directory).free
    if free < requested * 2:
        raise RuntimeError(f"insufficient free space: {free} < {requested * 2}")
    path: Path | None = None
    result: dict[str, Any] = {"mode": "write_probe", "bytes": requested,
                              "mount_device": device, "status": "failed"}
    try:
        fd, filename = tempfile.mkstemp(prefix="tandemx_t7_probe_", suffix=".bin", dir=directory)
        path = Path(filename)
        writer_hash = hashlib.sha256()
        began = time.monotonic()
        with os.fdopen(fd, "wb", buffering=0) as handle:
            for _ in range(requested // CHUNK_BYTES):
                block = os.urandom(CHUNK_BYTES)
                handle.write(block)
                writer_hash.update(block)
            handle.flush()
            os.fsync(handle.fileno())
            result["fsync_completed"] = True
            if hasattr(fcntl, "F_FULLFSYNC"):
                try:
                    fcntl.fcntl(handle.fileno(), fcntl.F_FULLFSYNC)
                    result["fullfsync_completed"] = True
                except OSError as exc:
                    if exc.errno not in {errno.EINVAL, errno.ENOTSUP, errno.ENOTTY}:
                        raise
                    result["fullfsync_completed"] = False
        result["write_fsync_seconds"] = round(time.monotonic() - began, 3)
        if path.stat().st_size != requested:
            raise IOError("written file has unexpected size")
        began = time.monotonic()
        read_hash, read_bytes = stream_sha256(path)
        result["readback_seconds"] = round(time.monotonic() - began, 3)
        result["write_sha256"] = writer_hash.hexdigest()
        result["read_sha256"] = read_hash
        if read_bytes != requested or read_hash != writer_hash.hexdigest():
            raise IOError("write/readback hash or size mismatch")
        if mount_identity()[0] != device:
            raise RuntimeError("T7 mount device changed during write probe")
        result["status"] = "pass"
        return result
    finally:
        if path is not None and path.exists():
            path.unlink()
            result["temporary_file_deleted"] = True
        else:
            result["temporary_file_deleted"] = path is None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["inputs", "write_probe"])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--probe-dir", type=Path)
    parser.add_argument("--probe-gib", type=int, default=3)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    try:
        if args.mode == "inputs":
            if args.manifest is None:
                parser.error("--manifest is required for inputs")
            result = check_inputs(args.manifest)
        else:
            if args.probe_dir is None:
                parser.error("--probe-dir is required for write_probe")
            result = write_probe(args.probe_dir, args.probe_gib)
        code = 0
    except (OSError, RuntimeError, ValueError, EOFError) as exc:
        result = {"mode": args.mode, "status": "failed", "error": repr(exc)}
        code = 1
    result["started_at"] = started
    result["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "output": str(args.output)}), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
