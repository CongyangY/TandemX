"""Archive compact, validated evidence from a completed whole-file sampling run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, iter_table

COMPLETED_SAMPLE_STATUSES = {"ok", "empty_random_sample"}


def _count_rows(path: Path, required: set[str]) -> int:
    total = 0
    for row in iter_table(path, required):
        try:
            count = int(row["read_count"])
        except ValueError as exc:
            raise ValueError(f"Invalid read count in {path}") from exc
        if count <= 0:
            raise ValueError(f"Nonpositive read count in {path}")
        total += count
    return total


def archive(source: Path, outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    plan_path = source / "sampling_plan.json"
    receipt_path = source / "sampling_receipt.json"
    if not plan_path.is_file() or not receipt_path.is_file():
        raise ValueError("Sampling plan and receipt are required")
    receipt = json.loads(receipt_path.read_text())
    samples = receipt.get("samples", [])
    sample_ids = [row.get("sample_id") for row in samples]
    if (
        receipt.get("complete") is not True
        or not samples
        or len(sample_ids) != len(set(sample_ids))
        or any(
            not sample_id or row.get("status") not in COMPLETED_SAMPLE_STATUSES
            for sample_id, row in zip(sample_ids, samples)
        )
        or receipt.get("plan_sha256") != digest_file(plan_path)
    ):
        raise ValueError("Require a complete, internally consistent sampling receipt")

    paths = [plan_path, receipt_path]
    for sample in samples:
        sample_id = sample["sample_id"]
        length_path = source / f"{sample_id}.length_histogram.tsv"
        joint_path = source / f"{sample_id}.joint_distribution.tsv"
        try:
            expected_reads = int(sample["read_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid receipt read count for {sample_id}") from exc
        if expected_reads < 0:
            raise ValueError(f"Invalid receipt read count for {sample_id}")
        if _count_rows(length_path, {"length_bp", "read_count"}) != expected_reads:
            raise ValueError(f"Length histogram total differs for {sample_id}")
        if _count_rows(
            joint_path,
            {"length_bin_kb", "gc_bin_percent", "mean_quality_bin_phred", "read_count"},
        ) != expected_reads:
            raise ValueError(f"Joint histogram total differs for {sample_id}")
        paths.extend([length_path, joint_path])

    outdir.mkdir(parents=True)
    manifest = []
    for path in paths:
        target = outdir / path.name
        shutil.copyfile(path, target)
        if target.read_bytes() != path.read_bytes():
            raise OSError(f"Archive copy differs: {path}")
        manifest.append(
            {
                "file": target.name,
                "source": str(path.resolve()),
                "sha256": digest_file(target),
                "bytes": target.stat().st_size,
            }
        )
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.outdir)
