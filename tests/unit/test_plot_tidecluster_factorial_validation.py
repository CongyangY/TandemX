import csv
import json

from benchmarks.scripts.plot_tidecluster_factorial_validation import (
    optional_float,
    plot,
)


def test_optional_float_preserves_missing_measurements() -> None:
    assert optional_float("") is None
    assert optional_float("None") is None
    assert optional_float("NA") is None
    assert optional_float("nan") is None
    assert optional_float("0") == 0.0
    assert optional_float("1.25") == 1.25


def test_plot_renders_all_six_frozen_fates_and_preserves_unavailable(tmp_path) -> None:
    result_dir = tmp_path / "result"
    result_dir.mkdir()
    rows = []
    for seed in (6401, 6402, 6403):
        for setting in ("default_primary", "matched_period_sensitivity"):
            rows.append(
                {
                    "seed": seed,
                    "setting": setting,
                    "status": "ok",
                    "array_recall": 0.8,
                    "array_precision": 0.9,
                    "base_union_recall": 0.75,
                    "base_union_precision": 0.85,
                    "matched_boundary_mae_bp": 2.5,
                    "matched_period_mae_bp": 1.0,
                    "cyclic_monomer_recall": 0.7,
                    "homologous_consensus_fraction": 0.6,
                    "tidehunter_wall_seconds": 10,
                    "tidehunter_maximum_rss_kb": 100000,
                    "clustering_wall_seconds": 20,
                    "clustering_maximum_rss_kb": 200000,
                }
            )
    rows[0]["status"] = "external_resource_failure_parent_v1"
    for name in (
        "array_recall",
        "array_precision",
        "base_union_recall",
        "base_union_precision",
        "matched_boundary_mae_bp",
        "matched_period_mae_bp",
        "cyclic_monomer_recall",
        "homologous_consensus_fraction",
    ):
        rows[0][name] = ""
    with (result_dir / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (result_dir / "run_receipt.json").write_text(json.dumps({"complete": True}))
    (result_dir / "independent_verification.json").write_text(
        json.dumps({"verification_passed": True, "failures": []})
    )
    (result_dir / "environment.json").write_text(json.dumps({"complete": True}))
    with (result_dir / "cell_fates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "seed",
                "setting",
                "status",
                "tidehunter_status",
                "clustering_status",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "seed": row["seed"],
                    "setting": row["setting"],
                    "status": row["status"],
                    "tidehunter_status": (
                        "imported_failed_exit_137" if row is rows[0] else "ok"
                    ),
                    "clustering_status": (
                        "not_started_dependency_failure" if row is rows[0] else "ok"
                    ),
                }
            )
    receipt = plot(result_dir, tmp_path / "figure")
    assert receipt["panel_count"] == 6
    assert receipt["frozen_run_count"] == 6
    assert receipt["successful_accuracy_run_count"] == 5
    assert receipt["unavailable_accuracy_run_count"] == 1
    assert receipt["svg_raster_image_element_count"] == 0
    assert receipt["svg_text_element_count"] > 0
