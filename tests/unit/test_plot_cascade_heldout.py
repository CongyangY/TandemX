from __future__ import annotations

import csv
import json
from pathlib import Path

from benchmarks.scripts.plot_cascade_heldout import plot


def test_plot_cascade_heldout_is_six_panel_and_source_backed(tmp_path: Path) -> None:
    evidence = Path("paper/evidence/cascade_native_screen_heldout_v1")
    provenance = plot(evidence, tmp_path / "figure")
    assert provenance["gate_status"] == "failed"
    assert set(provenance["failed_gate_names"]) == {
        "all_comparator_failed_runs",
        "tidehunter_runtime_geometric_mean_ratio",
    }
    rows = list(
        csv.DictReader(
            (tmp_path / "figure" / "panel_source.tsv").open(), delimiter="\t"
        )
    )
    assert {row["panel"] for row in rows} == set("ABCDEF")
    assert len([row for row in rows if row["panel"] == "D"]) == 48
    assert len([row for row in rows if row["panel"] == "F"]) == 12
    svg = (tmp_path / "figure" / "cascade_heldout.svg").read_text()
    assert "<image" not in svg
    assert svg.count("<g id=\"axes_") == 6
    assert svg.count("<text") > 20
    saved = json.loads(
        (tmp_path / "figure" / "figure_provenance.json").read_text()
    )
    assert saved["inputs"]["gates"]["sha256"]


def test_plot_cascade_validation_is_six_panel_and_records_all_gates(
    tmp_path: Path,
) -> None:
    evidence = Path("paper/evidence/cascade_gap_free_validation_v1")
    provenance = plot(evidence, tmp_path / "figure")
    assert provenance["gate_status"] == "passed"
    assert provenance["failed_gate_names"] == []
    assert provenance["evaluation_split"] == "validation"
    assert provenance["evaluation_seeds"] == ["2201"]
    rows = list(
        csv.DictReader(
            (tmp_path / "figure" / "panel_source.tsv").open(), delimiter="\t"
        )
    )
    assert {row["panel"] for row in rows} == set("ABCDEF")
    assert len([row for row in rows if row["panel"] == "D"]) == 16
    assert len([row for row in rows if row["panel"] == "F"]) == 14
    svg = (tmp_path / "figure" / "cascade_validation.svg").read_text()
    assert "<image" not in svg
    assert svg.count("<g id=\"axes_") == 6
    assert svg.count("<text") > 20
    assert "All 14 frozen gates passed" in (
        tmp_path / "figure" / "figure_legend.md"
    ).read_text()
