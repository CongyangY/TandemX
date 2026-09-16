"""Copy one independently verified SRA FASTQ chunk to a new output directory.

This never overwrites or deletes a prior attempt. The source must already have
a complete converter receipt; promotion rechecks its bytes and FASTQ records.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_fastq(path: Path, accession: str, start: int, end: int) -> tuple[int, int]:
    records = bases = 0
    with gzip.open(path, "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip(b"\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip(b"\r\n")
            if header.rstrip(b"\r\n") != f"@{accession}.{start + records}".encode():
                raise ValueError("FASTQ spot identifiers are not contiguous")
            if plus.rstrip(b"\r\n") != b"+" or not sequence or len(sequence) != len(quality):
                raise ValueError("FASTQ structure or quality length is invalid")
            records += 1
            bases += len(sequence)
    if records != end - start + 1:
        raise ValueError(f"FASTQ has {records} records; expected {end - start + 1}")
    return records, bases


def promote(source: Path, source_receipt: Path, target_dir: Path) -> dict[str, object]:
    receipt_bytes = source_receipt.read_bytes()
    receipt = json.loads(receipt_bytes)
    if receipt.get("status") != "complete" or receipt.get("vdb_dump_exit_code") != 0:
        raise ValueError("source receipt is not a completed SRA extraction")
    accession = receipt.get("run_accession")
    span = receipt.get("spot_range")
    if not isinstance(accession, str) or not isinstance(span, dict):
        raise ValueError("source accession or spot range is missing")
    start, end = span.get("start"), span.get("end")
    if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
        raise ValueError("source spot range is invalid")
    expected_name = f"spot_{start}_{end}.fastq.gz"
    if source.name != expected_name:
        raise ValueError("source FASTQ name does not match receipt range")
    if source.stat().st_size != receipt.get("fastq_bytes"):
        raise ValueError("source FASTQ size does not match receipt")
    source_sha = sha256_file(source)
    if source_sha != receipt.get("fastq_sha256"):
        raise ValueError("source FASTQ hash does not match receipt")
    records, bases = audit_fastq(source, accession, start, end)
    if records != receipt.get("record_count") or bases != receipt.get("total_bases"):
        raise ValueError("source FASTQ totals do not match receipt")

    target_dir.mkdir(parents=True, exist_ok=False)
    output = target_dir / expected_name
    partial = target_dir / f"{expected_name}.partial"
    target_receipt = target_dir / f"spot_{start}_{end}.receipt.json"
    try:
        with source.open("rb") as src, partial.open("xb") as dst:
            shutil.copyfileobj(src, dst, length=8 << 20)
            dst.flush()
            os.fsync(dst.fileno())
        if sha256_file(partial) != source_sha:
            raise ValueError("target readback hash mismatch")
        os.replace(partial, output)
        with target_receipt.open("xb") as handle:
            handle.write(receipt_bytes)
            handle.flush()
            os.fsync(handle.fileno())
        if sha256_file(output) != source_sha:
            raise ValueError("published target readback hash mismatch")
    except BaseException:
        partial.unlink(missing_ok=True)
        output.unlink(missing_ok=True)
        target_receipt.unlink(missing_ok=True)
        raise
    return {
        "status": "promoted_and_readback_verified",
        "source_fastq": str(source),
        "source_receipt": str(source_receipt),
        "target_fastq": str(output),
        "target_receipt": str(target_receipt),
        "fastq_sha256": source_sha,
        "receipt_sha256": hashlib.sha256(receipt_bytes).hexdigest(),
        "record_count": records,
        "total_bases": bases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--source-receipt", required=True, type=Path)
    parser.add_argument("--target-dir", required=True, type=Path)
    parser.add_argument("--promotion-receipt", required=True, type=Path)
    args = parser.parse_args()
    result = promote(args.source, args.source_receipt, args.target_dir)
    args.promotion_receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
