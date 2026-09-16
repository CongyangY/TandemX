"""Verify that an interrupted pipeline step does not keep writing in the background."""

from __future__ import annotations

import os
import json
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from tandemx.pipeline import PipelineConfig, run_command_with_live_logs, run_pipeline


def test_live_logs_normal_completion(tmp_path: Path) -> None:
    stdout = tmp_path / "stdout.log"
    stderr = tmp_path / "stderr.log"
    code = run_command_with_live_logs(
        [sys.executable, "-c", "import sys; print('done'); print('warning', file=sys.stderr)"],
        stdout,
        stderr,
    )
    assert code == 0
    assert stdout.read_text() == "done\n"
    assert stderr.read_text() == "warning\n"


def test_pipeline_records_interrupted_step(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">r1\nACGTACGT\n", encoding="utf-8")
    config = PipelineConfig(
        reads=(reads,), assembly=None, genome_size=1000, haploid_depth=None,
        outdir=tmp_path / "run", max_reads=None, max_read_bases=None,
        kmer_backend="python", steps=("discover",), min_period=2, max_period=10,
        top_periods=2, threads=1, resume=False, force=False, profile=False,
    )

    def interrupt(*_args: object, **_kwargs: object) -> int:
        raise KeyboardInterrupt

    monkeypatch.setattr("tandemx.pipeline.run_command_with_live_logs", interrupt)
    records, status = run_pipeline(config)
    assert status == 130
    assert len(records) == 1
    assert records[0].exit_status == 130
    assert not records[0].output_validated
    summary = json.loads((config.outdir / "pipeline_summary.json").read_text())
    assert summary[0]["notes"] == "interrupted_child_process_group_terminated"
    assert "exit_status=130" in (config.outdir / "pipeline.log").read_text()


@pytest.mark.skipif(os.name != "posix", reason="process groups require POSIX")
def test_sigint_terminates_active_step_process(tmp_path: Path) -> None:
    child_pid = tmp_path / "child.pid"
    stdout = tmp_path / "stdout.log"
    stderr = tmp_path / "stderr.log"
    wrapper = (
        "import sys\n"
        "from pathlib import Path\n"
        "from tandemx.pipeline import run_command_with_live_logs\n"
        "run_command_with_live_logs([sys.executable, '-c', "
        "'import os,sys,time,pathlib;pathlib.Path(sys.argv[1]).write_text(str(os.getpid()));time.sleep(30)', "
        "sys.argv[1]], Path(sys.argv[2]), Path(sys.argv[3]))\n"
    )
    parent = subprocess.Popen(
        [sys.executable, "-c", wrapper, str(child_pid), str(stdout), str(stderr)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        deadline = time.monotonic() + 8
        while not child_pid.is_file() and time.monotonic() < deadline:
            if parent.poll() is not None:
                pytest.fail("wrapper exited before child started")
            time.sleep(0.05)
        assert child_pid.is_file()
        pid = int(child_pid.read_text())
        os.kill(parent.pid, signal.SIGINT)
        assert parent.wait(timeout=8) != 0
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.05)
        else:
            pytest.fail("interrupted pipeline child remained alive")
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
