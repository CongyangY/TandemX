"""Download one immutable public file with parallel HTTP ranges and full hashes."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import time
import urllib.request
from urllib.parse import urlsplit


def planned_ranges(start: int, total: int, chunk_bytes: int) -> list[tuple[int, int]]:
    """Return inclusive, non-overlapping ranges covering ``start:total``."""
    if start < 0 or total < 1 or start > total or chunk_bytes < 1:
        raise ValueError("invalid range plan")
    return [
        (offset, min(total - 1, offset + chunk_bytes - 1))
        for offset in range(start, total, chunk_bytes)
    ]


def digest_file(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def _download_range(
    url: str,
    path: Path,
    start: int,
    end: int,
    total: int,
    timeout: float,
    attempts: int = 5,
) -> dict[str, int | str]:
    expected = end - start + 1
    partial = path.with_name(path.name + ".partial")
    if path.exists():
        if path.stat().st_size != expected:
            raise ValueError(f"completed range has wrong size: {path}")
        return {"path": str(path), "start": start, "end": end, "bytes": expected}
    for attempt in range(1, attempts + 1):
        have = partial.stat().st_size if partial.exists() else 0
        if have > expected:
            raise ValueError(f"partial range exceeds expected size: {partial}")
        if have == expected:
            partial.rename(path)
            break
        request_start = start + have
        request = urllib.request.Request(
            url,
            headers={"Range": f"bytes={request_start}-{end}", "Accept-Encoding": "identity"},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                content_range = response.headers.get("Content-Range", "")
                expected_header = f"bytes {request_start}-{end}/{total}"
                if response.status != 206 or content_range != expected_header:
                    raise ValueError(
                        f"server returned {response.status} {content_range!r}; "
                        f"expected 206 {expected_header!r}"
                    )
                with partial.open("ab") as handle:
                    while block := response.read(1024 * 1024):
                        have += len(block)
                        if have > expected:
                            raise ValueError(f"range transfer exceeded its budget: {partial}")
                        handle.write(block)
            if have != expected:
                raise OSError(f"short range transfer: {have} of {expected} bytes")
            partial.rename(path)
            break
        except (OSError, TimeoutError, http.client.HTTPException) as error:
            if attempt == attempts:
                raise
            print(
                f"retry range={start}-{end} attempt={attempt} "
                f"bytes={have} error={error}",
                flush=True,
            )
            time.sleep(2)
    return {"path": str(path), "start": start, "end": end, "bytes": expected}


def assemble_and_verify(
    output: Path,
    prefix: Path | None,
    parts: list[Path],
    expected_bytes: int,
    expected_md5: str,
) -> dict[str, int | str | bool]:
    assembling = output.with_name(output.name + ".assembling")
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    written = 0
    with assembling.open("wb") as destination:
        for source in ([prefix] if prefix is not None else []) + parts:
            if source is None:
                continue
            with source.open("rb") as handle:
                while block := handle.read(8 * 1024 * 1024):
                    destination.write(block)
                    md5.update(block)
                    sha256.update(block)
                    written += len(block)
    if written != expected_bytes or md5.hexdigest() != expected_md5:
        raise ValueError(
            "assembled download failed size or publisher MD5 verification; "
            "all parts were retained"
        )
    os.replace(assembling, output)
    return {
        "bytes": written,
        "md5": md5.hexdigest(),
        "sha256": sha256.hexdigest(),
        "publisher_md5_verified": True,
    }


def download(
    url: str,
    output: Path,
    expected_bytes: int,
    expected_md5: str,
    workers: int,
    chunk_mib: int,
    existing_prefix: Path | None = None,
    timeout: float = 90.0,
) -> dict:
    parsed = urlsplit(url)
    if parsed.scheme != "https" or not parsed.hostname or parsed.query:
        raise ValueError("download URL must be an immutable HTTPS object without a query")
    if expected_bytes < 1 or not re.fullmatch(r"[0-9a-f]{32}", expected_md5):
        raise ValueError("expected size and MD5 are required")
    if not 1 <= workers <= 32 or chunk_mib < 1:
        raise ValueError("workers or chunk size outside the supported range")
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    if output.exists() and output.stat().st_size == expected_bytes:
        md5, sha256 = digest_file(output)
        if md5 != expected_md5:
            raise ValueError("existing full-sized output has the wrong MD5")
        return {
            "url": url,
            "output": str(output),
            "bytes": expected_bytes,
            "md5": md5,
            "sha256": sha256,
            "publisher_md5_verified": True,
            "transfer": "existing_complete_reverified",
            "elapsed_seconds": time.time() - started,
        }

    prefix = existing_prefix
    if prefix is not None and prefix.resolve() == output.resolve() and not prefix.exists():
        prior_backup = output.with_name(output.name + ".prefix")
        if prior_backup.exists():
            prefix = prior_backup
    if prefix is not None and prefix.resolve() == output.resolve() and prefix.exists():
        backup = output.with_name(output.name + ".prefix")
        if backup.exists():
            if backup.stat().st_size != output.stat().st_size:
                raise ValueError("existing prefix backup conflicts with current partial output")
            output.unlink()
        else:
            output.rename(backup)
        prefix = backup
    prefix_bytes = prefix.stat().st_size if prefix is not None and prefix.exists() else 0
    if prefix_bytes > expected_bytes:
        raise ValueError("existing prefix exceeds expected object size")

    parts_dir = output.with_name(output.name + ".parts")
    parts_dir.mkdir(exist_ok=True)
    ranges = planned_ranges(prefix_bytes, expected_bytes, chunk_mib * 1024 * 1024)
    part_specs = [
        (parts_dir / f"part_{start:020d}_{end:020d}", start, end)
        for start, end in ranges
    ]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _download_range, url, path, start, end, expected_bytes, timeout
            ): (start, end)
            for path, start, end in part_specs
        }
        completed = 0
        for future in as_completed(futures):
            future.result()
            completed += 1
            if completed == len(futures) or completed % max(1, workers) == 0:
                print(f"completed_ranges={completed}/{len(futures)}", flush=True)

    verification = assemble_and_verify(
        output,
        prefix if prefix_bytes else None,
        [path for path, _, _ in part_specs],
        expected_bytes,
        expected_md5,
    )
    receipt = {
        "url": url,
        "output": str(output),
        "expected_bytes": expected_bytes,
        "expected_md5": expected_md5,
        "workers": workers,
        "chunk_mib": chunk_mib,
        "existing_prefix": str(prefix) if prefix_bytes else None,
        "existing_prefix_bytes": prefix_bytes,
        "range_count": len(part_specs),
        **verification,
        "transfer": "parallel_ranges_complete",
        "elapsed_seconds": time.time() - started,
    }
    receipt_path = output.with_name(output.name + ".download_receipt.json")
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--expected-md5", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-mib", type=int, default=128)
    parser.add_argument("--existing-prefix", type=Path)
    parser.add_argument("--timeout", type=float, default=90.0)
    args = parser.parse_args()
    receipt = download(
        args.url,
        args.output,
        args.expected_bytes,
        args.expected_md5,
        args.workers,
        args.chunk_mib,
        args.existing_prefix,
        args.timeout,
    )
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
