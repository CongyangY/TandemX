from __future__ import annotations

import csv
import json
from pathlib import Path

from benchmarks.scripts.plot_quantify_depth_gated_validation import plot


def test_plot_depth_gated_validation_is_six_panel_editable_and_paired(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    outdir = tmp_path / "figure"
    provenance = plot(
        root / "paper/evidence/quantify_depth_gated_validation_v1",
        outdir,
    )
    svg = (outdir / "quantify_depth_gated_validation.svg").read_text()
    assert svg.count("<image") == 0
    assert svg.count("<text") >= 80
    assert svg.count("<g id=\"axes_") == 6
    assert (outdir / "quantify_depth_gated_validation.pdf").stat().st_size > 10_000
    assert (outdir / "quantify_depth_gated_validation.png").stat().st_size > 10_000
    rows = list(csv.DictReader((outdir / "panel_source.tsv").open(), delimiter="\t"))
    resource_rows = [row for row in rows if row["panel"] == "F"]
    resource_keys = {(row["method"], row["stratum"]) for row in resource_rows}
    assert len(resource_rows) == 108
    assert len(resource_keys) == 54
    assert json.loads((outdir / "figure_provenance.json").read_text()) == provenance
    assert len(provenance["outputs"]) == 5
