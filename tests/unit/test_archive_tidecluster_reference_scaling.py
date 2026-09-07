from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_tidecluster_reference_scaling import RUN_FILES, archive


def _write_run(root: Path, sample_id: str, input_bases: int) -> Path:
    root.mkdir(parents=True)
    for name in RUN_FILES:
        if name == "finalization_receipt.json":
            continue
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder\n")
    normalized = root / "normalized/normalized_arrays.tsv"
    with normalized.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "family_id",
                "representative_selection_source",
                "copy_number_source",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "family_id": "TRC_1",
                "representative_selection_source": "exact_intermediate_interval",
                "copy_number_source": "exact_tidehunter_interval",
            }
        )
    result = {
        "complete": True,
        "sample_id": sample_id,
        "summary": {
            "predicted_array_count": 1,
            "predicted_family_count": 1,
            "predicted_positive_sequence_count": 1,
            "predicted_union_bp": 100,
            "predicted_union_base_fraction": 100 / input_bases,
        },
        "internal_gnu_time": {
            "tidehunter": {"wall_seconds": 2, "maximum_rss_kb": 100},
            "clustering": {"wall_seconds": 3, "maximum_rss_kb": 200},
        },
        "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
        "warning": "sampled_real_reference_no_whole_genome_context",
    }
    (root / "result.json").write_text(json.dumps(result))
    (root / "finalization_environment.json").write_text(
        json.dumps({"sample_sha256": f"sha-{sample_id}"})
    )
    tracked = []
    for name in RUN_FILES:
        if name == "finalization_receipt.json":
            continue
        path = root / name
        tracked.append(
            {"file": name, "bytes": path.stat().st_size, "sha256": digest_file(path)}
        )
    (root / "finalization_receipt.json").write_text(
        json.dumps({"complete": True, "external_stages_reused": True, "files": tracked})
    )
    return root


def test_archive_tidecluster_reference_scaling_copies_verified_runs(tmp_path: Path) -> None:
    samples = {
        "ten": {"total_bases": 10_000_000, "sha256": "sha-ten"},
        "hundred": {"total_bases": 100_000_000, "sha256": "sha-hundred"},
    }
    sampling = tmp_path / "sampling.json"
    sampling.write_text(json.dumps({"complete": True, "samples": samples}))
    runs = [
        _write_run(tmp_path / "ten", "ten", 10_000_000),
        _write_run(tmp_path / "hundred", "hundred", 100_000_000),
    ]
    result = archive(sampling, runs, tmp_path / "archive")
    assert result["complete"] is True
    assert len(result["files"]) == len(RUN_FILES) * 2 + 4
    with (tmp_path / "archive/scaling_summary.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4
    assert {row["stage"] for row in rows} == {"tidehunter", "clustering"}
    assert {row["accuracy"] for row in rows} == {
        "not_assessed_without_independent_real_array_and_family_truth"
    }


def test_archive_tidecluster_reference_scaling_rejects_changed_run_file(
    tmp_path: Path,
) -> None:
    samples = {
        "ten": {"total_bases": 10_000_000, "sha256": "sha-ten"},
        "hundred": {"total_bases": 100_000_000, "sha256": "sha-hundred"},
    }
    sampling = tmp_path / "sampling.json"
    sampling.write_text(json.dumps({"complete": True, "samples": samples}))
    runs = [
        _write_run(tmp_path / "ten", "ten", 10_000_000),
        _write_run(tmp_path / "hundred", "hundred", 100_000_000),
    ]
    (runs[0] / "result.json").write_text("changed\n")
    with pytest.raises(ValueError, match="changed"):
        archive(sampling, runs, tmp_path / "archive")
