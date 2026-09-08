import hashlib
import json

import pytest

from benchmarks.scripts.archive_tidecluster_factorial import (
    aggregate_settings,
    archive_description,
    resource_context,
    validate_figure,
)


def test_aggregate_settings_preserves_seed_range_and_resource_values() -> None:
    rows = [
        {
            "seed": "1",
            "setting": "default",
            "array_recall": "0.8",
            "array_precision": "0.5",
            "base_union_recall": "0.7",
            "base_union_precision": "0.6",
            "matched_boundary_mae_bp": "2",
            "matched_period_mae_bp": "1",
            "cyclic_monomer_recall": "0.4",
            "homologous_consensus_fraction": "0.3",
            "tidehunter_wall_seconds": "10",
            "tidehunter_maximum_rss_kb": "100",
            "clustering_wall_seconds": "20",
            "clustering_maximum_rss_kb": "200",
        },
        {
            "seed": "2",
            "setting": "default",
            "array_recall": "1.0",
            "array_precision": "0.7",
            "base_union_recall": "0.9",
            "base_union_precision": "0.8",
            "matched_boundary_mae_bp": "4",
            "matched_period_mae_bp": "3",
            "cyclic_monomer_recall": "0.6",
            "homologous_consensus_fraction": "0.5",
            "tidehunter_wall_seconds": "12",
            "tidehunter_maximum_rss_kb": "120",
            "clustering_wall_seconds": "22",
            "clustering_maximum_rss_kb": "220",
        },
    ]
    result = aggregate_settings(rows)["default"]
    assert result["run_count"] == 2
    assert result["seeds"] == [1, 2]
    assert result["metrics"]["array_recall"] == {
        "mean": 0.9,
        "minimum": 0.8,
        "maximum": 1.0,
    }
    assert result["metrics"]["clustering_maximum_rss_kb"]["mean"] == 210


def test_validate_figure_rejects_raster_svg(tmp_path) -> None:
    directory = tmp_path / "figures_v1"
    directory.mkdir()
    svg = directory / "figure.svg"
    svg.write_text("<svg><text>test</text></svg>")
    record = {
        "bytes": svg.stat().st_size,
        "sha256": hashlib.sha256(svg.read_bytes()).hexdigest(),
    }
    provenance = {
        "complete": True,
        "panel_count": 6,
        "frozen_run_count": 6,
        "svg_raster_image_element_count": 0,
        "svg_text_element_count": 1,
        "outputs": {svg.name: record},
    }
    (directory / "figure_provenance.json").write_text(json.dumps(provenance))
    assert validate_figure(tmp_path)["panel_count"] == 6
    provenance["svg_raster_image_element_count"] = 1
    (directory / "figure_provenance.json").write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match="raster image"):
        validate_figure(tmp_path)


def test_validate_figure_prefers_v2_and_requires_identical_panel_source(tmp_path) -> None:
    for version, png_text in (("figures_v1", "failed"), ("figures_v2", "accepted")):
        directory = tmp_path / version
        directory.mkdir()
        outputs = {}
        for name, text in (
            ("figure.svg", "<svg><text>test</text></svg>"),
            ("tidecluster_factorial_validation.png", png_text),
            ("panel_source.tsv", "panel\tvalue\nA\t1\n"),
        ):
            path = directory / name
            path.write_text(text)
            outputs[name] = {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        (directory / "figure_provenance.json").write_text(
            json.dumps(
                {
                    "complete": True,
                    "panel_count": 6,
                    "frozen_run_count": 6,
                    "svg_raster_image_element_count": 0,
                    "svg_text_element_count": 1,
                    "outputs": outputs,
                }
            )
        )
    result = validate_figure(tmp_path)
    assert result["qa_history"]["accepted_render"] == "figures_v2"
    (tmp_path / "figures_v2/panel_source.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="output changed"):
        validate_figure(tmp_path)


def test_failure_archive_description_does_not_claim_accuracy_outputs() -> None:
    description = archive_description(False)
    assert "failed before accuracy evaluation" in description
    assert "no figure or accuracy summary" in description
    partial = archive_description(True, False)
    assert "all frozen cell and stage fates" in partial
    assert "unavailable cells contain no" in partial


def test_resource_context_retains_failed_container_peak(tmp_path) -> None:
    profile = tmp_path / "profile"
    profile.mkdir()
    (profile / "stages.tsv").write_text(
        "stage\texit_code\twall_seconds\nfirst\t1\t2.0\n"
    )
    run_dir = tmp_path / "seed1/default"
    run_dir.mkdir(parents=True)
    (run_dir / "tidehunter.gnu_time.txt").write_text(
        "User time (seconds): 1.0\n"
        "System time (seconds): 0.2\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
        "Maximum resident set size (kbytes): 7776184\n"
        "Exit status: 1\n"
    )
    context = resource_context(tmp_path)
    assert context["failed_profile_stages"][0]["stage"] == "first"
    internal = context["internal_gnu_time"][
        "seed1/default/tidehunter.gnu_time.txt"
    ]
    assert internal["maximum_rss_kb"] == 7776184
