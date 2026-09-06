"""A small real executable run verifies benchmark output plumbing."""

import json
import shutil
from pathlib import Path

import pytest
import yaml

from benchmarks.challenge.run import run_suite
from benchmarks.challenge.schema import read_table
from tandemx.discover.rust_backend import rust_backend_available


def test_challenge_real_tandemx_cli_and_empty_truth(tmp_path: Path) -> None:
    if not shutil.which("tandemx") or not rust_backend_available():
        pytest.skip("The Rust release and installed TandemX CLI are needed")
    config = {"seeds": {"development": [1101], "heldout": [3101]}, "tools": {"tandemx": "tandemx"},
              "repetitions": 1, "timeout_seconds": 30, "period_range": [30, 200],
              "minimum_span_bp": 100, "minimum_iou": 0.5,
              "scenarios": [{"name": "small", "read_count": 8, "read_length": 800, "copies": 4,
                             "period": 37, "families": 2},
                            {"name": "negative", "read_count": 4, "read_length": 800, "copies": 4,
                             "period": 37, "positive_fraction": 0}]}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    outdir = tmp_path / "run"
    assert run_suite(path, outdir, "development") == 0
    assert json.loads((outdir / "validation.json").read_text())["successful_runs"] == 2
    summary = read_table(outdir / "summary.tsv")
    positive = next(row for row in summary if row["scenario"] == "small")
    negative = next(row for row in summary if row["scenario"] == "negative")
    assert float(positive["array_recall"]) > 0
    assert negative["array_recall"] == "NA"
    assert negative["sequence_family_recall"] == "NA"
    assert positive["deterministic"] == "not_tested_single_run"
    with pytest.raises(FileExistsError):
        run_suite(path, outdir, "development")
