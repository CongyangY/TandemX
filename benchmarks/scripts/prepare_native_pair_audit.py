"""Copy small archived native-read subsets with closed-loop integrity checks.

This is an input audit, not a read-to-assembly mapping or accuracy benchmark.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fastq_stats(path: Path) -> tuple[int, int]:
    count = bases = 0
    with gzip.open(path, "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip(b"\r\n")
            separator = handle.readline()
            quality = handle.readline().rstrip(b"\r\n")
            if (not header.startswith(b"@") or not separator.startswith(b"+")
                    or not sequence or len(sequence) != len(quality)):
                raise ValueError(f"invalid FASTQ record {count + 1}: {path}")
            count += 1
            bases += len(sequence)
    if count == 0:
        raise ValueError(f"empty FASTQ: {path}")
    return count, bases


def prepare(config: Path, outdir: Path, reserve_bytes: int = 10 * 1024**3) -> dict:
    spec = json.loads(config.read_text())
    if spec.get("schema_version") != 1 or not spec.get("subsets"):
        raise ValueError("invalid audit config")
    outdir.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(outdir).free
    total = sum(row["bytes"] for row in spec["subsets"])
    if free - total < reserve_bytes:
        raise ValueError("internal-disk reserve would be breached")
    results = []
    for row in spec["subsets"]:
        source = Path(row["source"])
        if source.stat().st_size != row["bytes"]:
            raise ValueError(f"source size mismatch: {source}")
        target = outdir / f"{row['dataset']}.fastq.gz"
        partial = target.with_suffix(target.suffix + ".partial")
        digest = hashlib.sha256()
        try:
            with source.open("rb") as inp, partial.open("wb") as out:
                for chunk in iter(lambda: inp.read(1024 * 1024), b""):
                    digest.update(chunk)
                    out.write(chunk)
                out.flush()
                os.fsync(out.fileno())
            if digest.hexdigest() != row["sha256"]:
                raise ValueError(f"source SHA-256 mismatch: {source}")
            if sha256_file(partial) != row["sha256"]:
                raise ValueError(f"copy readback mismatch: {partial}")
            count, bases = fastq_stats(partial)
            if (count, bases) != (row["read_count"], row["total_bases"]):
                raise ValueError(f"FASTQ totals mismatch: {source}")
            partial.replace(target)
        finally:
            partial.unlink(missing_ok=True)
        results.append({"dataset": row["dataset"], "source": str(source),
                        "copy": str(target), "bytes": row["bytes"],
                        "sha256": row["sha256"], "read_count": count,
                        "total_bases": bases, "gzip_trailer_checked": True,
                        "source_and_copy_hash_verified": True})
    receipt = {"schema_version": 1, "scope": "bounded_native_read_subset_integrity_only",
               "config_sha256": sha256_file(config), "subsets": results,
               "full_source_current_readback": "not_performed",
               "native_interval_support": "not_measured"}
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def extract_selected(source: Path, expected_sha256: str, read_ids: set[str],
                     out_fastq_gz: Path) -> dict:
    """Preserve exact original FASTQ record bytes for named archive reads."""
    if not read_ids:
        raise ValueError("no read IDs requested")
    if sha256_file(source) != expected_sha256:
        raise ValueError("source SHA-256 mismatch")
    selected: dict[str, bytes] = {}
    total_reads = 0
    with gzip.open(source, "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if (not header.startswith(b"@") or not separator.startswith(b"+")
                    or not sequence or len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n"))):
                raise ValueError(f"invalid FASTQ record {total_reads + 1}")
            total_reads += 1
            read_id = header[1:].split(None, 1)[0].decode("ascii")
            if read_id in read_ids:
                if read_id in selected:
                    raise ValueError(f"duplicate selected read ID: {read_id}")
                selected[read_id] = header + sequence + separator + quality
    if set(selected) != read_ids:
        raise ValueError(f"missing requested reads: {sorted(read_ids - set(selected))}")
    out_fastq_gz.parent.mkdir(parents=True, exist_ok=True)
    with out_fastq_gz.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as out:
            for read_id in sorted(selected):
                out.write(selected[read_id])
    if fastq_stats(out_fastq_gz)[0] != len(read_ids):
        raise ValueError("selected FASTQ readback failed")
    return {"source": str(source), "source_sha256": expected_sha256,
            "source_gzip_trailer_checked": True, "source_read_count": total_reads,
            "selected_fastq": str(out_fastq_gz), "selected_fastq_sha256": sha256_file(out_fastq_gz),
            "reads": [{"id": read_id, "record_sha256": hashlib.sha256(selected[read_id]).hexdigest(),
                       "bases": len(selected[read_id].splitlines()[1])} for read_id in sorted(selected)]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.config, args.outdir), indent=2))


if __name__ == "__main__":
    main()
