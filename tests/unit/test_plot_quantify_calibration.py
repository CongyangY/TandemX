from __future__ import annotations

import json
from pathlib import Path

from benchmarks.scripts.plot_quantify_calibration import plot


def test_plot_quantify_calibration_is_six_panel_and_editable(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    outdir = tmp_path / "figure"
    provenance = plot(
        root / "paper/evidence/quantify_calibration_development_v1",
        outdir,
    )
    svg = (outdir / "quantify_calibration.svg").read_text()
    assert svg.count("<image") == 0
    assert svg.count("<text") >= 50
    assert svg.count("<g id=\"axes_") == 6
    assert (outdir / "quantify_calibration.pdf").stat().st_size > 10_000
    assert (outdir / "quantify_calibration.png").stat().st_size > 10_000
    saved = json.loads((outdir / "figure_provenance.json").read_text())
    assert saved == provenance
    assert len(saved["outputs"]) == 5
