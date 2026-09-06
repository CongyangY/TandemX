"""Archive compact evidence from one completed full-file FASTQ QC run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import string

from benchmarks.challenge.schema import digest_file, iter_table


FILES = ("qc.json", "length_histogram.tsv", "joint_distribution.tsv", "base_quality_histogram.tsv")


def _sum_column(path: Path, column: str, required: set[str]) -> int:
    total = 0
    for row in iter_table(path, required | {column}):
        try:
            value = int(row[column])
        except ValueError as exc:
            raise ValueError(f"Invalid {column} in {path}") from exc
        if value <= 0:
            raise ValueError(f"Nonpositive {column} in {path}")
        total += value
    return total


def archive(source: Path, outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    paths = [source / name for name in FILES]
    if any(not path.is_file() for path in paths):
        raise ValueError("Complete QC receipt and three distribution tables are required")
    receipt = json.loads(paths[0].read_text())
    try:
        read_count = int(receipt["read_count"])
        total_bases = int(receipt["total_bases"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("QC receipt has invalid denominators") from exc
    input_hash = receipt.get("input_sha256", "")
    if (
        receipt.get("complete") is not True
        or receipt.get("fastq_records_valid") is not True
        or receipt.get("gzip_trailer_checked") is not True
        or receipt.get("exact_duplicate_read_ids") != 0
        or read_count <= 0
        or total_bases <= 0
        or len(input_hash) != 64
        or any(character not in string.hexdigits for character in input_hash)
    ):
        raise ValueError("Require complete duplicate-free full-file FASTQ QC")
    if _sum_column(paths[1], "read_count", {"length_bp"}) != read_count:
        raise ValueError("Length-histogram denominator differs from QC receipt")
    if _sum_column(
        paths[2], "read_count", {"length_bin_kb", "gc_bin_percent", "mean_quality_bin_phred"}
    ) != read_count:
        raise ValueError("Joint-distribution denominator differs from QC receipt")
    if _sum_column(paths[3], "base_count", {"phred"}) != total_bases:
        raise ValueError("Base-quality denominator differs from QC receipt")

    outdir.mkdir(parents=True)
    manifest = []
    for path in paths:
        target = outdir / path.name
        shutil.copyfile(path, target)
        source_hash = digest_file(path)
        if target.stat().st_size != path.stat().st_size or digest_file(target) != source_hash:
            raise OSError(f"Archive copy differs: {path}")
        manifest.append({
            "file": target.name,
            "source": str(path.resolve()),
            "sha256": source_hash,
            "bytes": target.stat().st_size,
        })
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
