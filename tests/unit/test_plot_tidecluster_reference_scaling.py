from __future__ import annotations

import csv
import json
from pathlib import Path

from benchmarks.scripts.plot_tidecluster_reference_scaling import plot


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_plot_tidecluster_reference_scaling_writes_editable_six_panel_figure(
    tmp_path: Path,
) -> None:
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    rows = []
    for sample_id, input_mb, arrays, families in (
        ("ten", 10, 8, 3),
        ("hundred", 100, 80, 12),
    ):
        for stage, seconds, rss in (
            ("tidehunter", input_mb * 2, input_mb * 1000),
            ("clustering", input_mb, input_mb * 2000),
        ):
            rows.append(
                {
                    "sample_id": sample_id,
                    "input_bases": input_mb * 1_000_000,
                    "input_mb": input_mb,
                    "stage": stage,
                    "wall_seconds": seconds,
                    "maximum_rss_kb": rss,
                    "predicted_array_count": arrays,
                    "predicted_family_count": families,
                    "predicted_positive_sequence_count": input_mb // 10,
                    "predicted_union_bp": input_mb * 20_000,
                    "predicted_union_base_fraction": 0.02,
                    "exact_intermediate_interval_count": arrays - 2,
                    "clipped_interval_count": 1,
                    "merged_interval_count": 1,
                    "exact_copy_number_count": arrays - 1,
                    "unavailable_copy_number_count": 1,
                    "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
                    "warning": "sampled_real_reference_no_whole_genome_context",
                }
            )
        normalized = [
            {"period": 40 + index * 10} for index in range(arrays)
        ]
        _write_tsv(
            evidence / f"runs/{sample_id}/normalized/normalized_arrays.tsv",
            normalized,
        )
    _write_tsv(evidence / "scaling_summary.tsv", rows)
    (evidence / "summary.json").write_text(
        json.dumps(
            {
                "complete": True,
                "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
            }
        )
    )
    (evidence / "archive_manifest.json").write_text(
        json.dumps({"complete": True, "files": []})
    )
    receipt = plot(evidence, tmp_path / "figures")
    assert receipt["complete"] is True
    assert receipt["panel_count"] == 6
    assert receipt["svg_text_nodes"] > 20
    assert receipt["svg_image_nodes"] == 0
    assert (tmp_path / "figures/tidecluster_morex_reference_scaling.pdf").stat().st_size > 1000
