from pathlib import Path

import pytest

from benchmarks.scripts.run_tidecluster_docker_reference import (
    container_path,
    parse_gnu_time,
)


def test_container_path_requires_mounted_root(tmp_path: Path) -> None:
    mounted = tmp_path / "root"
    sample = mounted / "results/sample.fa"
    sample.parent.mkdir(parents=True)
    sample.write_text(">x\nACGT\n")
    assert container_path(sample, mounted, "/input") == "/input/results/sample.fa"
    with pytest.raises(ValueError, match="outside"):
        container_path(tmp_path / "other.fa", mounted, "/input")


def test_parse_successful_gnu_time(tmp_path: Path) -> None:
    path = tmp_path / "time.txt"
    path.write_text(
        "User time (seconds): 1.25\n"
        "System time (seconds): 0.50\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50\n"
        "Maximum resident set size (kbytes): 2048\n"
        "Exit status: 0\n"
    )
    result = parse_gnu_time(path)
    assert result["wall_seconds"] == 62.5
    assert result["maximum_rss_kb"] == 2048


def test_parse_gnu_time_can_retain_failed_stage_resources(tmp_path: Path) -> None:
    path = tmp_path / "failed.time"
    path.write_text(
        "User time (seconds): 1.2\n"
        "System time (seconds): 0.3\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.50\n"
        "Maximum resident set size (kbytes): 7776184\n"
        "Exit status: 1\n"
    )
    with pytest.raises(ValueError, match="successful stage"):
        parse_gnu_time(path)
    result = parse_gnu_time(path, require_success=False)
    assert result["exit_status"] == 1
    assert result["maximum_rss_kb"] == 7776184
