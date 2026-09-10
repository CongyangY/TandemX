import csv
import json
import random
import subprocess
import sys
from pathlib import Path

from tandemx.io.validators import validate_project


def test_valid_negative_discovery_and_pipeline(tmp_path: Path) -> None:
    rng = random.Random(33)
    reads = tmp_path / "reads.fa"
    reads.write_text("".join(f">r{i}\n{''.join(rng.choices('ACGT', k=1500))}\n" for i in range(6)))
    out = tmp_path / "out"
    result = subprocess.run([sys.executable, "-m", "tandemx.cli", "run", "--reads", str(reads),
                             "--outdir", str(out), "--steps", "discover,quantify,validate",
                             "--min-period", "30", "--max-period", "500", "--threads", "1"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    summary = json.loads((out / "discover" / "discovery_summary.json").read_text())
    assert summary["status"] == "no_families"
    assert summary["processed_reads"] == 6
    assert summary["candidate_count"] == summary["family_count"] == 0
    assert (out / "discover" / "monomers.fa").read_text() == ""
    assert validate_project(out / "discover")
    with (out / "pipeline_summary.tsv").open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    assert next(r for r in records if r["step"] == "quantify")["notes"] == "skipped_no_discovered_families"
    assert not (out / "quantify" / "copy_number.tsv").exists()
    (out / "quantify").mkdir()
    stale = out / "quantify" / "copy_number.tsv"
    stale.write_text("stale positive output\n")
    before = (stale.read_bytes(), stale.stat().st_mtime_ns)
    refused = subprocess.run([sys.executable, "-m", "tandemx.cli", "run", "--reads", str(reads),
                              "--outdir", str(out), "--steps", "discover,quantify,validate", "--no-resume",
                              "--min-period", "30", "--max-period", "500", "--threads", "1"],
                             capture_output=True, text=True)
    assert refused.returncode != 0
    assert "--no-resume refuses a populated output directory" in refused.stderr
    assert (stale.read_bytes(), stale.stat().st_mtime_ns) == before
    rerun = subprocess.run([sys.executable, "-m", "tandemx.cli", "run", "--reads", str(reads),
                            "--outdir", str(out), "--steps", "discover,quantify,validate", "--force",
                            "--min-period", "30", "--max-period", "500", "--threads", "1"],
                           capture_output=True, text=True)
    assert rerun.returncode == 0, rerun.stderr
    assert not stale.exists()


def test_zero_families_after_support_filter_and_optional_collapse(tmp_path: Path) -> None:
    rng = random.Random(55)
    unit = "".join(rng.choices("ACGT", k=61))
    reads = tmp_path / "reads.fa"
    reads.write_text(f">only\n{unit * 10}\n")
    out = tmp_path / "discover"
    result = subprocess.run([sys.executable, "-m", "tandemx.cli", "discover", "--reads", str(reads),
                             "--outdir", str(out), "--min-support-reads", "2", "--min-period", "30",
                             "--max-period", "100", "--collapse-redundant-families", "--no-progress"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads((out / "discovery_summary.json").read_text())["candidate_count"] == 1
    assert validate_project(out)
