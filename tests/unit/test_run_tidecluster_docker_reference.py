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
