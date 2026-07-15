#!/usr/bin/env python3
"""Validate and summarize an external-tool benchmark for publication review."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path
from typing import Iterable


PAIRWISE_FIELDS = [
    "dataset_id",
    "total_bases",
    "monomer_length",
    "error_rate",
    "competitor",
    "tandemx_runtime_seconds",
    "competitor_runtime_seconds",
    "tandemx_speedup",
    "tandemx_peak_memory_mb",
    "competitor_peak_memory_mb",
    "tandemx_memory_reduction_fraction",
    "memory_winner",
]

MACRO_FIELDS = [
    "tool",
    "macro_recall",
    "macro_precision",
    "macro_false_positive_rate",
    "macro_monomer_length_mae_bp",
]

CHART_FIELDS = ["tool", "runtime_seconds", "peak_memory_mb"]
DATASET_TABLE_FIELDS = [
    "dataset_id",
    "observed_read_count",
    "expected_read_length",
    "positive_reads",
    "truth_id_coverage",
    "source_mode",
]
PAIRWISE_SQL = """\
WITH tandemx AS (
    SELECT * FROM benchmark_summary WHERE tool = 'tandemx'
)
SELECT
    tandemx.dataset_id,
    tandemx.total_bases,
    tandemx.monomer_length,
    tandemx.error_rate,
    competitor.tool AS competitor,
    tandemx.median_runtime_seconds AS tandemx_runtime_seconds,
    competitor.median_runtime_seconds AS competitor_runtime_seconds,
    competitor.median_runtime_seconds / tandemx.median_runtime_seconds AS tandemx_speedup,
    tandemx.median_peak_memory_mb AS tandemx_peak_memory_mb,
    competitor.median_peak_memory_mb AS competitor_peak_memory_mb,
    1.0 - tandemx.median_peak_memory_mb / competitor.median_peak_memory_mb
        AS tandemx_memory_reduction_fraction,
    CASE
        WHEN tandemx.median_peak_memory_mb < competitor.median_peak_memory_mb THEN 'tandemx'
        ELSE competitor.tool
    END AS memory_winner
FROM tandemx
JOIN benchmark_summary AS competitor USING (dataset_id)
WHERE competitor.tool IN ('trf', 'tidehunter')
ORDER BY tandemx.total_bases ASC,
         CASE competitor.tool WHEN 'trf' THEN 1 ELSE 2 END
"""
MACRO_SQL = """\
SELECT
    tool,
    AVG(recall) AS macro_recall,
    AVG(precision) AS macro_precision,
    AVG(false_positive_rate) AS macro_false_positive_rate,
    AVG(monomer_length_mae_bp) AS macro_monomer_length_mae_bp
FROM benchmark_summary
GROUP BY tool
ORDER BY CASE tool WHEN 'tandemx' THEN 1 WHEN 'trf' THEN 2 ELSE 3 END
"""
DATASET_TABLE_SQL = """\
SELECT
    dataset_id,
    observed_read_count,
    expected_read_length,
    positive_reads,
    truth_id_coverage,
    source_mode
FROM dataset_manifest
ORDER BY observed_read_count ASC, expected_read_length ASC
"""
CHART_SQL = """\
SELECT
    CASE tool
        WHEN 'tandemx' THEN 'TandemX'
        WHEN 'trf' THEN 'TRF'
        WHEN 'tidehunter' THEN 'TideHunter'
    END AS tool,
    median_runtime_seconds AS runtime_seconds,
    median_peak_memory_mb AS peak_memory_mb
FROM benchmark_summary
WHERE dataset_id = 'long_noisy_421'
ORDER BY median_runtime_seconds ASC
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--dataset-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_inputs(
    summary_rows: list[dict[str, str]], manifest_rows: list[dict[str, str]]
) -> dict[str, object]:
    required_tools = {"tandemx", "trf", "tidehunter"}
    datasets = {row["dataset_id"] for row in manifest_rows}
    observed_pairs = {(row["dataset_id"], row["tool"]) for row in summary_rows}
    expected_pairs = {(dataset, tool) for dataset in datasets for tool in required_tools}
    issues: list[str] = []
    if not summary_rows:
        issues.append("summary_is_empty")
    if not manifest_rows:
        issues.append("dataset_manifest_is_empty")
    if observed_pairs != expected_pairs:
        issues.append("dataset_tool_matrix_incomplete")
    if any(row["source_mode"] != "reused_existing_files" for row in manifest_rows):
        issues.append("not_all_inputs_reused")
    if any(float(row["truth_id_coverage"]) != 1.0 for row in manifest_rows):
        issues.append("truth_id_coverage_incomplete")
    if any(row["observed_read_count"] != row["truth_row_count"] for row in manifest_rows):
        issues.append("read_truth_row_count_mismatch")
    if any(int(row["successful_runs"]) != 3 for row in summary_rows):
        issues.append("not_all_three_repetitions_succeeded")
    if any(row["deterministic_predictions"] != "true" for row in summary_rows):
        issues.append("prediction_digest_mismatch")
    maximum_runtime_cv = max(
        (float(row["runtime_cv"]) for row in summary_rows), default=0.0
    )
    return {
        "status": "passed" if not issues else "failed",
        "issues": issues,
        "dataset_count": len(datasets),
        "dataset_tool_combinations": len(summary_rows),
        "successful_runs": sum(int(row["successful_runs"]) for row in summary_rows),
        "all_prediction_digests_deterministic": not any(
            row["deterministic_predictions"] != "true" for row in summary_rows
        ),
        "maximum_runtime_cv": maximum_runtime_cv,
        "runtime_variability_note": (
            "The maximum is dominated by process cold-start on the 0.4 Mb dataset; "
            "publication runtime claims should emphasize the larger workloads."
            if maximum_runtime_cv > 0.1
            else "All dataset-tool runtime CV values are at most 10%."
        ),
    }


def summary_database(summary_rows: list[dict[str, str]]) -> sqlite3.Connection:
    database = sqlite3.connect(":memory:")
    database.execute(
        "CREATE TABLE benchmark_summary ("
        "dataset_id TEXT, tool TEXT, total_bases INTEGER, monomer_length INTEGER, error_rate REAL, "
        "median_runtime_seconds REAL, median_peak_memory_mb REAL, recall REAL, precision REAL, "
        "false_positive_rate REAL, monomer_length_mae_bp REAL)"
    )
    database.executemany(
        "INSERT INTO benchmark_summary VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                row["dataset_id"], row["tool"], int(row["total_bases"]),
                int(row["monomer_length"]), float(row["error_rate"]),
                float(row["median_runtime_seconds"]), float(row["median_peak_memory_mb"]),
                float(row["recall"]), float(row["precision"]),
                float(row["false_positive_rate"]), float(row["monomer_length_mae_bp"]),
            )
            for row in summary_rows
        ],
    )
    return database


def build_pairwise_rows(summary_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    database = summary_database(summary_rows)
    try:
        rows = [dict(zip(PAIRWISE_FIELDS, row, strict=True)) for row in database.execute(PAIRWISE_SQL)]
    finally:
        database.close()
    for row in rows:
        row["tandemx_runtime_seconds"] = f"{float(row['tandemx_runtime_seconds']):.6f}"
        row["competitor_runtime_seconds"] = f"{float(row['competitor_runtime_seconds']):.6f}"
        row["tandemx_speedup"] = f"{float(row['tandemx_speedup']):.6f}"
        row["tandemx_peak_memory_mb"] = f"{float(row['tandemx_peak_memory_mb']):.3f}"
        row["competitor_peak_memory_mb"] = f"{float(row['competitor_peak_memory_mb']):.3f}"
        row["tandemx_memory_reduction_fraction"] = (
            f"{float(row['tandemx_memory_reduction_fraction']):.6f}"
        )
    return rows


def build_macro_rows(summary_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    database = summary_database(summary_rows)
    try:
        output = [dict(zip(MACRO_FIELDS, row, strict=True)) for row in database.execute(MACRO_SQL)]
    finally:
        database.close()
    for row in output:
        for field in MACRO_FIELDS[1:]:
            row[field] = f"{float(row[field]):.6f}"
    return output


def build_chart_rows(summary_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    database = summary_database(summary_rows)
    try:
        return [dict(zip(CHART_FIELDS, row, strict=True)) for row in database.execute(CHART_SQL)]
    finally:
        database.close()


def build_dataset_table_rows(manifest_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    database = sqlite3.connect(":memory:")
    try:
        database.execute(
            "CREATE TABLE dataset_manifest (dataset_id TEXT, observed_read_count INTEGER, "
            "expected_read_length INTEGER, positive_reads INTEGER, truth_id_coverage REAL, "
            "source_mode TEXT)"
        )
        database.executemany(
            "INSERT INTO dataset_manifest VALUES (?, ?, ?, ?, ?, ?)",
            [
                (
                    row["dataset_id"], int(row["observed_read_count"]),
                    int(row["expected_read_length"]), int(row["positive_reads"]),
                    float(row["truth_id_coverage"]), row["source_mode"],
                )
                for row in manifest_rows
            ],
        )
        return [
            dict(zip(DATASET_TABLE_FIELDS, row, strict=True))
            for row in database.execute(DATASET_TABLE_SQL)
        ]
    finally:
        database.close()


def main() -> int:
    args = parse_args()
    summary_rows = read_tsv(args.summary)
    manifest_rows = read_tsv(args.dataset_manifest)
    validation = validate_inputs(summary_rows, manifest_rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "pairwise_comparison.tsv", build_pairwise_rows(summary_rows), PAIRWISE_FIELDS)
    write_tsv(args.outdir / "macro_accuracy.tsv", build_macro_rows(summary_rows), MACRO_FIELDS)
    write_tsv(args.outdir / "long_noisy_runtime.tsv", build_chart_rows(summary_rows), CHART_FIELDS)
    write_tsv(
        args.outdir / "dataset_table.tsv",
        build_dataset_table_rows(manifest_rows),
        DATASET_TABLE_FIELDS,
    )
    (args.outdir / "analysis_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    return 0 if validation["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
