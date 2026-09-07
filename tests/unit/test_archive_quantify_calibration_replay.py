from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.archive_quantify_calibration_replay import (
    METHODS,
    archive,
    compare_runs,
)


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


def _build_run(root: Path, runtime_scale: float, metric_value: float = 2.0) -> None:
    snapshot = root / "source_snapshot"
    source = snapshot / "tandemx" / "quantify" / "mvp.py"
    source.parent.mkdir(parents=True)
    source.write_text("source\n")
    hashes = {"tandemx/quantify/mvp.py": digest_file(source)}
    _write_json(
        root / "environment.json",
        {
            "git_head": "a" * 40,
            "source_snapshot": str(snapshot),
            "file_hashes": hashes,
            "source_digest": hashlib.sha256(
                json.dumps(hashes, sort_keys=True).encode()
            ).hexdigest(),
        },
    )
    _write_json(
        root / "execution.json",
        {
            "exit_code": 0,
            "timed_out": False,
            "runtime_seconds": 20 * runtime_scale,
            "peak_rss_mib": 100,
        },
    )
    executions = []
    metrics = []
    summaries = []
    for index, method in enumerate(METHODS):
        executions.append(
            {
                "seed": 1,
                "condition_id": "c1",
                "method": method,
                "status": "ok",
                "runtime_seconds": (index + 1) * runtime_scale,
                "peak_rss_mib": 10 + index,
            }
        )
        metrics.append(
            {
                "seed": 1,
                "condition_id": "c1",
                "method": method,
                "family_id": "f1",
                "estimated_copy_number": metric_value,
            }
        )
        summaries.append(
            {
                "method": method,
                "coverage": 5,
                "error_model": "none",
                "mean_absolute_relative_error": 0.1,
                "median_runtime_seconds": (index + 1) * runtime_scale,
                "median_peak_rss_mib": 10 + index,
            }
        )
        output = root / "runs" / "s1" / "c1" / method / "output"
        output.mkdir(parents=True)
        (output / "copy_number.tsv").write_text("family_id\testimate\nf1\t2\n")
    write_table(root / "executions.tsv", executions, list(executions[0]))
    write_table(root / "metrics.tsv", metrics, list(metrics[0]))
    write_table(root / "summary.tsv", summaries, list(summaries[0]))
    _write_json(
        root / "validation.json",
        {
            "complete": True,
            "executions": len(executions),
            "successful_executions": len(executions),
            "family_conditions": len(metrics),
        },
    )
    (root / "frozen_config.yaml").write_text("benchmark_id: test\n")
    (root / "stdout.log").write_text("")
    (root / "stderr.log").write_text("")
    controls = root / "controls" / "s1"
    controls.mkdir(parents=True)
    (controls / "single_copy_kmers.tsv").write_text(
        "kmer\texpected_copy_number\nAAAAA\t1\n"
    )
    _write_json(controls / "receipt.json", {"selected_controls": 1})


def test_replay_archive_requires_exact_scientific_parity(tmp_path: Path) -> None:
    baseline, replay = tmp_path / "baseline", tmp_path / "replay"
    _build_run(baseline, 1.0)
    _build_run(replay, 0.5)
    checked = compare_runs(baseline, replay)
    assert checked["comparison"]["status"] == "exact_scientific_output_parity"
    assert checked["comparison"]["copy_number_files_compared"] == len(METHODS)
    assert checked["comparison"]["driver_speedup"] == 2
    assert all(row["median_runtime_speedup"] == 2 for row in checked["resources"])

    outdir = tmp_path / "archive"
    receipt = archive(baseline, replay, outdir)
    assert receipt["metrics_byte_identical"] is True
    manifest = json.loads((outdir / "archive_manifest.json").read_text())
    assert len(manifest) == 13
    for row in manifest:
        assert digest_file(outdir / row["file"]) == row["sha256"]


@pytest.mark.parametrize("fault", ["metrics", "copy_number", "control", "source"])
def test_replay_archive_rejects_changed_evidence(tmp_path: Path, fault: str) -> None:
    baseline, replay = tmp_path / "baseline", tmp_path / "replay"
    _build_run(baseline, 1.0)
    _build_run(replay, 0.5)
    if fault == "metrics":
        text = (replay / "metrics.tsv").read_text().replace("\t2.0\n", "\t3.0\n", 1)
        (replay / "metrics.tsv").write_text(text)
    elif fault == "copy_number":
        path = replay / "runs" / "s1" / "c1" / METHODS[0] / "output" / "copy_number.tsv"
        path.write_text(path.read_text().replace("\t2\n", "\t3\n"))
    elif fault == "control":
        (replay / "controls" / "s1" / "single_copy_kmers.tsv").write_text(
            "kmer\texpected_copy_number\nCCCCC\t1\n"
        )
    else:
        (replay / "source_snapshot" / "tandemx" / "quantify" / "mvp.py").write_text(
            "changed\n"
        )
    with pytest.raises(ValueError):
        compare_runs(baseline, replay)
