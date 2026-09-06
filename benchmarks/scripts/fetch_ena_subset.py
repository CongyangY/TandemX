#!/usr/bin/env python3
"""Retrieve a bounded FASTQ prefix from a public ENA run with provenance.

The complete remote FASTQ checksum is NOT verified for a prefix. This limitation
is explicit in the receipt. Prefix samples are engineering pilots, not unbiased
whole-library abundance samples. Only HTTPS ENA-hosted FASTQ URLs are accepted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


class BoundedReader:
    def __init__(self, stream, limit: int):
        self.stream, self.limit, self.bytes_read = stream, limit, 0

    def read(self, size: int = -1) -> bytes:
        remaining = self.limit - self.bytes_read
        if remaining <= 0:
            raise ValueError("Compressed download budget exceeded before completing the requested reads")
        data = self.stream.read(min(remaining, size) if size >= 0 else remaining)
        self.bytes_read += len(data)
        return data


def extract_fastq_prefix(source, destination, read_count: int) -> tuple[int, int, list[str]]:
    if read_count < 1:
        raise ValueError("Read count must be positive")
    bases = 0
    first_headers: list[str] = []
    seen: set[bytes] = set()
    for index in range(read_count):
        lines = [source.readline() for _ in range(4)]
        header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
        if not header.startswith(b"@") or not separator.startswith(b"+"):
            raise ValueError(f"Malformed or incomplete FASTQ record {index + 1}")
        if not sequence or len(sequence) != len(quality) or set(sequence.upper()) - set(b"ACGTN"):
            raise ValueError(f"Invalid FASTQ sequence/quality at record {index + 1}")
        identifier = header[1:].split()[0]
        if identifier in seen:
            raise ValueError(f"Duplicate read identifier: {identifier.decode()}")
        seen.add(identifier)
        if len(first_headers) < 3:
            first_headers.append(header.decode())
        destination.write(b">" + header[1:] + b"\n" + sequence.upper() + b"\n")
        bases += len(sequence)
    return read_count, bases, first_headers


def fetch_subset(accession: str, read_count: int, outdir: Path, max_download_bytes: int) -> dict:
    if not re.fullmatch(r"[ESD]RR\d+", accession):
        raise ValueError("Use a run accession (ERR/SRR/DRR followed by digits)")
    if read_count < 1 or max_download_bytes < 1:
        raise ValueError("Read count and download budget must be positive")
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    query = urllib.parse.urlencode({"accession": accession, "result": "read_run", "format": "tsv",
                                   "fields": "run_accession,sample_accession,scientific_name,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"})
    api = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query
    with urllib.request.urlopen(api, timeout=30) as response:
        metadata_text = response.read(1024 * 1024).decode()
    (outdir / "ena_metadata.tsv").write_text(metadata_text)
    rows = list(csv.DictReader(io.StringIO(metadata_text), delimiter="\t"))
    if len(rows) != 1 or rows[0]["run_accession"] != accession:
        raise ValueError("ENA did not return exactly the requested run")
    row = rows[0]
    paths = row["fastq_ftp"].split(";")
    if len(paths) != 1:
        raise ValueError("This long-read prefix pilot expects one FASTQ file; paired/multipart runs need a different sampler")
    url = "https://" + paths[0]
    if urllib.parse.urlsplit(url).hostname != "ftp.sra.ebi.ac.uk":
        raise ValueError(f"Unexpected ENA FASTQ hostname: {url}")
    if read_count > int(row["read_count"]):
        raise ValueError("Requested more reads than recorded by ENA")
    for label, identifier in (("run", accession), ("sample", row["sample_accession"])):
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/ena/browser/api/xml/{identifier}", timeout=30) as response:
            (outdir / f"ena_{label}.xml").write_bytes(response.read(2 * 1024 * 1024))
    temp = outdir / "reads.fa.partial"
    with urllib.request.urlopen(url, timeout=60) as response:
        bounded = BoundedReader(response, max_download_bytes)
        with gzip.GzipFile(fileobj=bounded) as fastq, temp.open("wb") as fasta:
            count, bases, headers = extract_fastq_prefix(fastq, fasta, read_count)
    output = outdir / "reads.fa"
    temp.rename(output)
    with output.open("rb") as handle:
        sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    receipt = {"accession": accession, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
               "metadata_url": api, "source_url": url, "metadata": row,
               "sampling": "first_N_records_in_archive_order", "requested_reads": read_count,
               "observed_reads": count, "observed_bases": bases, "compressed_bytes_read": bounded.bytes_read,
               "max_compressed_bytes": max_download_bytes, "output_sha256": sha256,
               "output_bytes": output.stat().st_size, "first_headers": headers,
               "source_full_file_md5": row["fastq_md5"], "source_full_file_md5_verified": False,
               "warning": "prefix_is_not_a_random_sample;remote_gzip_trailer_and_full_file_checksum_not_checked"}
    (outdir / "subset_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--read-count", type=int, default=1000)
    parser.add_argument("--max-download-mib", type=int, default=256)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = fetch_subset(args.accession, args.read_count, args.outdir, args.max_download_mib * 1024 * 1024)
    except (ValueError, OSError, EOFError) as error:
        parser.exit(2, f"error: {error}\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
