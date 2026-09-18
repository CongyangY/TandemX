"""Combine exact extracted Macadamia FASTQ records without changing record bytes."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.map_macadamia_native_context import sha256


def combine(extraction_receipt: Path, output: Path) -> dict:
    source = json.loads(extraction_receipt.read_text())
    if source["status"] != "complete" or source["distinct_zmw_verified"]:
        raise ValueError("Expected complete provisional source records")
    expected = {item["id"]: item["record_sha256"]
                for run in source["runs"] for item in run["reads"]}
    if len(expected) != source["qualifying_record_count"]:
        raise ValueError("Extraction receipt has duplicate or missing IDs")
    records = {}
    for run in source["runs"]:
        path = Path(run["selected_fastq"])
        if sha256(path) != run["selected_fastq_sha256"]:
            raise ValueError("Extracted source FASTQ changed")
        with gzip.open(path, "rb") as stream:
            while header := stream.readline():
                record = header + stream.readline() + stream.readline() + stream.readline()
                read_id = header[1:].split(None, 1)[0].decode("ascii")
                if read_id not in expected or read_id in records:
                    raise ValueError("Unexpected or repeated original read")
                if hashlib.sha256(record).hexdigest() != expected[read_id]:
                    raise ValueError("Original FASTQ record bytes changed")
                records[read_id] = record
    if set(records) != set(expected):
        raise ValueError("Missing original read")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as stream:
            for read_id in sorted(records):
                stream.write(records[read_id])
    result = {"status": "combined_exact_original_records", "record_count": len(records),
              "distinct_zmw_verified": False,
              "extraction_receipt_sha256": sha256(extraction_receipt),
              "combined_fastq_sha256": sha256(output),
              "per_record_sha256": expected}
    (output.parent / "combination_receipt.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extraction-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(combine(args.extraction_receipt, args.output), sort_keys=True))


if __name__ == "__main__":
    main()
