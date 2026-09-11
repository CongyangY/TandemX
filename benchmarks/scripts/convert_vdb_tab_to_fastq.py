"""Convert a bounded ``vdb-dump`` READ/QUALITY tab stream into FASTQ.

This is a benchmark-acquisition utility, not a TandemX production input
sampler.  It accepts only a declared contiguous SRA spot range and records that
range in its receipt; a range must never be described as a random sample.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import BinaryIO


BIOLOGICAL = "SRA_READ_TYPE_BIOLOGICAL"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _quality_to_fastq(value: str, expected_length: int) -> bytes:
    fields = value.split(",")
    if len(fields) != expected_length:
        raise ValueError(
            f"QUALITY has {len(fields)} values but READ_LEN is {expected_length}"
        )
    qualities = bytearray()
    for token in fields:
        token = token.strip()
        if not re.fullmatch(r"\d+", token):
            raise ValueError(f"QUALITY value is not an unsigned integer: {token!r}")
        score = int(token)
        if not 0 <= score <= 93:
            raise ValueError(f"QUALITY value outside [0, 93]: {score}")
        qualities.append(score + 33)
    return bytes(qualities)


def convert(
    source: BinaryIO,
    output: Path,
    receipt: Path,
    *,
    accession: str,
    start_spot: int,
    end_spot: int,
) -> dict[str, object]:
    """Stream strict four-column VDB tab rows to a deterministic gzip FASTQ.

    The required tab columns, in order, are READ, QUALITY, READ_LEN and
    READ_TYPE.  The record identifier is derived from the declared accession
    and absolute spot number, which is unique within one contiguous range.
    """
    if not re.fullmatch(r"[ESD]RR\d+", accession):
        raise ValueError("accession must be a public ENA/SRA run accession")
    if start_spot < 1 or end_spot < start_spot:
        raise ValueError("spot range is invalid")
    if output.exists() or receipt.exists():
        raise ValueError("refusing to overwrite output or receipt")
    partial = output.with_name(output.name + ".partial")
    if partial.exists():
        raise ValueError(f"partial output already exists: {partial}")
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    source_digest = hashlib.sha256()
    records = bases = 0
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "failed_partial_retained",
        "run_accession": accession,
        "spot_range": {"start": start_spot, "end": end_spot},
        "input_format": "vdb_dump_SEQUENCE_READ_QUALITY_READ_LEN_READ_TYPE_tab",
        "selection": "fixed_contiguous_spot_range_not_random_sampling",
        "quality_conversion": "integer_0_to_93_to_fastq_phred_plus_33",
        "only_read_type": BIOLOGICAL,
        "converter_sha256": _sha256_file(Path(__file__)),
    }
    try:
        with partial.open("wb") as raw, gzip.GzipFile(
            filename="", fileobj=raw, mode="wb", compresslevel=1, mtime=0
        ) as compressed:
            for line_number, raw_line in enumerate(source, 1):
                source_digest.update(raw_line)
                if raw_line in (b"\n", b"\r\n"):
                    continue
                try:
                    line = raw_line.decode("ascii").rstrip("\r\n")
                except UnicodeDecodeError as exc:
                    raise ValueError(f"row {line_number}: tab data is not ASCII") from exc
                columns = line.split("\t")
                if len(columns) != 4:
                    raise ValueError(f"row {line_number}: expected exactly 4 tab columns")
                read, quality, length_text, read_type = columns
                if not length_text.isdigit() or int(length_text) <= 0:
                    raise ValueError(f"row {line_number}: READ_LEN is not a positive integer")
                read_length = int(length_text)
                if len(read) != read_length:
                    raise ValueError(
                        f"row {line_number}: READ has {len(read)} bases but READ_LEN is {read_length}"
                    )
                if read_type != BIOLOGICAL:
                    raise ValueError(f"row {line_number}: disallowed READ_TYPE {read_type!r}")
                fastq_quality = _quality_to_fastq(quality, read_length)
                spot = start_spot + records
                if spot > end_spot:
                    raise ValueError("source contains more rows than the declared spot range")
                identifier = f"@{accession}.{spot}".encode("ascii")
                compressed.write(identifier + b"\n" + read.encode("ascii") + b"\n+\n" + fastq_quality + b"\n")
                records += 1
                bases += read_length
        expected_records = end_spot - start_spot + 1
        if records != expected_records:
            raise ValueError(f"received {records} rows but declared range requires {expected_records}")
        partial.rename(output)
        result.update(
            status="complete",
            record_count=records,
            total_bases=bases,
            source_tab_sha256=source_digest.hexdigest(),
            fastq_sha256=_sha256_file(output),
            fastq_bytes=output.stat().st_size,
            ids_are_unique=True,
        )
        receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result
    except Exception as exc:
        result.update(
            error=str(exc),
            completed_records=records,
            completed_bases=bases,
            source_tab_sha256_so_far=source_digest.hexdigest(),
            partial_fastq_path=str(partial),
            partial_fastq_bytes=partial.stat().st_size if partial.exists() else 0,
        )
        receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="VDB tab file; omit to read stdin")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--start-spot", required=True, type=int)
    parser.add_argument("--end-spot", required=True, type=int)
    args = parser.parse_args()
    if args.input is None:
        convert(sys.stdin.buffer, args.output, args.receipt, accession=args.accession,
                start_spot=args.start_spot, end_spot=args.end_spot)
    else:
        with args.input.open("rb") as source:
            convert(source, args.output, args.receipt, accession=args.accession,
                    start_spot=args.start_spot, end_spot=args.end_spot)


if __name__ == "__main__":
    main()
