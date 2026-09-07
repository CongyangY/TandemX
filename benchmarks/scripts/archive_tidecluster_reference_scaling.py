#!/usr/bin/env python3
"""Validate and archive compact TideCluster real-reference scaling evidence."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, write_table


RUN_FILES = (
    "tc_tidehunter.gff3",
    "tc_clustering.gff3_1.gff3",
    "tc_clustering.gff3",
    "tc_cmd_args.json",
    "tc_consensus/consensus_sequences_all.fasta",
    "result.json",
    "finalization_environment.json",
    "finalization_receipt.json",
    "profile/receipt.json",
    "profile/stages.tsv",
    "tidehunter.gnu_time.txt",
    "clustering.gnu_time.txt",
    "normalized/summary.json",
    "normalized/normalized_arrays.tsv",
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _verified_finalization(run_dir: Path) -> tuple[dict[str, object], list[dict[str, str]]]:
    receipt = json.loads((run_dir / "finalization_receipt.json").read_text())
    if receipt.get("complete") is not True or receipt.get("external_stages_reused") is not True:
        raise ValueError(f"Run lacks a completed no-rerun finalization: {run_dir}")
    tracked = {str(item["file"]): item for item in receipt.get("files", [])}
    expected_tracked = set(RUN_FILES) - {"finalization_receipt.json"}
    if not expected_tracked.issubset(tracked):
        raise ValueError(f"Finalization receipt lacks required files: {run_dir}")
    for item in tracked.values():
        path = run_dir / str(item["file"])
        if (
            not path.is_file()
            or path.stat().st_size != item["bytes"]
            or digest_file(path) != item["sha256"]
        ):
            raise ValueError(f"Finalized TideCluster file changed: {path}")
    missing = [name for name in RUN_FILES if not (run_dir / name).is_file()]
    if missing:
        raise ValueError(f"Compact TideCluster evidence is incomplete: {missing}")
    result = json.loads((run_dir / "result.json").read_text())
    if (
        result.get("complete") is not True
        or result.get("accuracy")
        != "not_assessed_without_independent_real_array_and_family_truth"
    ):
        raise ValueError(f"TideCluster result lacks its real-data evidence boundary: {run_dir}")
    rows = _read_tsv(run_dir / "normalized/normalized_arrays.tsv")
    if len(rows) != result["summary"]["predicted_array_count"]:
        raise ValueError(f"Normalized row count differs from result summary: {run_dir}")
    return result, rows


def _provenance_counts(
    run_dir: Path, rows: list[dict[str, str]]
) -> Counter[str]:
    if rows and "representative_selection_source" in rows[0]:
        return Counter(row["representative_selection_source"] for row in rows)
    final_rows = _read_gff_keys(run_dir / "tc_clustering.gff3")
    intermediate_rows = _read_gff_keys(run_dir / "tc_clustering.gff3_1.gff3")
    if final_rows != intermediate_rows:
        raise ValueError(
            "Legacy normalized table lacks selection provenance and GFF intervals differ"
        )
    return Counter({"exact_intermediate_interval": len(rows)})


def _read_gff_keys(path: Path) -> set[tuple[str, int, int]]:
    keys = set()
    for line in path.open(encoding="utf-8"):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 9:
            raise ValueError(f"Malformed GFF3 row: {path}")
        key = (fields[0], int(fields[3]), int(fields[4]))
        if key in keys:
            raise ValueError(f"Duplicate GFF3 interval: {path}:{key}")
        keys.add(key)
    return keys


def archive(
    sampling_receipt: Path,
    run_dirs: list[Path],
    outdir: Path,
) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    sampling = json.loads(sampling_receipt.read_text())
    if sampling.get("complete") is not True:
        raise ValueError("Reference-window sampling receipt is incomplete")
    samples = sampling.get("samples", {})
    if len(run_dirs) < 2:
        raise ValueError("Scaling archive requires at least two completed sizes")

    outdir.mkdir(parents=True)
    shutil.copyfile(sampling_receipt, outdir / "sampling_receipt.json")
    scaling_rows = []
    source_entries = []
    seen_samples = set()
    for run_dir in run_dirs:
        run_dir = run_dir.resolve()
        result, rows = _verified_finalization(run_dir)
        sample_id = str(result["sample_id"])
        if sample_id in seen_samples or sample_id not in samples:
            raise ValueError(f"Duplicate or unknown sample in scaling archive: {sample_id}")
        seen_samples.add(sample_id)
        expected_bases = samples[sample_id]["total_bases"]
        finalization = json.loads((run_dir / "finalization_environment.json").read_text())
        if finalization.get("sample_sha256") != samples[sample_id]["sha256"]:
            raise ValueError(f"Run and sampling receipt disagree: {sample_id}")
        target_root = outdir / "runs" / sample_id
        for name in RUN_FILES:
            source = run_dir / name
            target = target_root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            if source.stat().st_size != target.stat().st_size or digest_file(source) != digest_file(target):
                raise OSError(f"Archive copy differs: {source}")
            source_entries.append(
                {
                    "file": target.relative_to(outdir).as_posix(),
                    "source": str(source),
                    "bytes": target.stat().st_size,
                    "sha256": digest_file(target),
                }
            )
        provenance = _provenance_counts(run_dir, rows)
        copy_sources = Counter(row["copy_number_source"] for row in rows)
        summary = result["summary"]
        resources = result["internal_gnu_time"]
        for stage in ("tidehunter", "clustering"):
            scaling_rows.append(
                {
                    "sample_id": sample_id,
                    "input_bases": expected_bases,
                    "input_mb": expected_bases / 1_000_000,
                    "stage": stage,
                    "wall_seconds": resources[stage]["wall_seconds"],
                    "maximum_rss_kb": resources[stage]["maximum_rss_kb"],
                    "predicted_array_count": summary["predicted_array_count"],
                    "predicted_family_count": summary["predicted_family_count"],
                    "predicted_positive_sequence_count": summary[
                        "predicted_positive_sequence_count"
                    ],
                    "predicted_union_bp": summary["predicted_union_bp"],
                    "predicted_union_base_fraction": summary[
                        "predicted_union_base_fraction"
                    ],
                    "exact_intermediate_interval_count": provenance[
                        "exact_intermediate_interval"
                    ],
                    "clipped_interval_count": provenance[
                        "unique_family_consistent_overlap_after_clipping"
                    ],
                    "merged_interval_count": provenance[
                        "longest_family_consistent_overlap_after_merge"
                    ],
                    "exact_copy_number_count": copy_sources[
                        "exact_tidehunter_interval"
                    ],
                    "unavailable_copy_number_count": copy_sources[
                        "unavailable_after_interval_merge_or_resolution"
                    ],
                    "accuracy": result["accuracy"],
                    "warning": result["warning"],
                }
            )
    scaling_rows.sort(key=lambda row: (row["input_bases"], row["stage"]))
    write_table(outdir / "scaling_summary.tsv", scaling_rows, list(scaling_rows[0]))
    summary = {
        "schema_version": 1,
        "complete": True,
        "sample_ids": sorted(seen_samples),
        "run_count": len(seen_samples),
        "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
        "warning": "descriptive_nested_real_reference_windows_not_whole_genome_accuracy_truth",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (outdir / "README.md").write_text(
        "# TideCluster MorexV3 reference-window scaling\n\n"
        "This compact archive contains completed TideCluster 1.21.2 runs on "
        "deterministic nested real-reference windows. Each run retains its "
        "external outputs, normalized calls, two-stage resource measurements, "
        "finalization receipt and source-environment hashes.\n\n"
        "The calls are descriptive. The sampled reference has no independent "
        "array/family truth, the windows do not preserve whole-chromosome "
        "context, and the rows therefore do not measure accuracy. Copy number "
        "is unavailable when TideCluster merged or resolved a TideHunter "
        "interval. `scaling_summary.tsv` is the figure/table source.\n",
        encoding="utf-8",
    )
    generated = (
        "sampling_receipt.json",
        "scaling_summary.tsv",
        "summary.json",
        "README.md",
    )
    manifest = {
        "complete": True,
        "files": [
            *source_entries,
            *[
                {
                    "file": name,
                    "source": "generated_by_archive_tidecluster_reference_scaling",
                    "bytes": (outdir / name).stat().st_size,
                    "sha256": digest_file(outdir / name),
                }
                for name in generated
            ],
        ],
    }
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampling-receipt", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, action="append", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.sampling_receipt, args.run_dir, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
