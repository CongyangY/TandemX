import json
from pathlib import Path
import sys

import pytest

from benchmarks.scripts.profile_stage_resources import load_manifest, profile


def test_stage_profiler_records_success_and_failure(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"stages": [
        {"name": "ok", "command": [sys.executable, "-c", "x=bytearray(2000000); print(len(x))"]},
        {"name": "bad", "command": [sys.executable, "-c", "raise SystemExit(7)"]},
        {"name": "unreached", "command": [sys.executable, "-c", "print('no')"]},
    ]}))
    out = tmp_path / "out"
    assert profile(manifest, out, 0.05) == 1
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["complete"] is False
    assert receipt["completed_stage_count"] == 2
    stages = (out / "stages.tsv").read_text().splitlines()
    assert len(stages) == 3
    assert "\t7\t" in stages[2]
    assert (out / "ok.stdout.log").read_text().strip() == "2000000"
    assert (out / "samples.tsv").read_text().startswith("stage\telapsed_seconds\trss_mib")


def test_manifest_rejects_shell_strings_and_duplicate_names(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"stages": [{"name": "x", "command": "echo unsafe"}]}))
    with pytest.raises(ValueError, match="string array"):
        load_manifest(path)
    path.write_text(json.dumps({"stages": [
        {"name": "x", "command": ["true"]},
        {"name": "x", "command": ["true"]},
    ]}))
    with pytest.raises(ValueError, match="duplicate"):
        load_manifest(path)
