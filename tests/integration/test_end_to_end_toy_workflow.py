from __future__ import annotations

import subprocess
import os
import sys
from pathlib import Path


def test_end_to_end_toy_workflow_script(tmp_path: Path) -> None:
    output_dir = tmp_path / "toy_results"
    result = subprocess.run(
        ["bash", "examples/toy/run_toy_workflow.sh", str(output_dir)],
        check=False,
        text=True,
        capture_output=True,
        env={**os.environ, "TANDEMX_CMD": f"{sys.executable} -m tandemx.cli"},
    )

    assert result.returncode == 0, result.stderr
    assert "Toy workflow complete" in result.stdout
    expected = [
        output_dir / "simulated" / "reads.fa",
        output_dir / "simulated" / "assembly.fa",
        output_dir / "discover" / "families.tsv",
        output_dir / "discover" / "monomers.fa",
        output_dir / "quantify" / "copy_number.tsv",
        output_dir / "locate" / "repeat_density.bedgraph",
        output_dir / "locate" / "arrays.bed",
        output_dir / "locate" / "assembly_vs_read_cn.tsv",
        output_dir / "probe" / "probes.fa",
        output_dir / "probe" / "probes.rank.tsv",
        output_dir / "probe" / "in_silico_fish.tsv",
        output_dir / "visualize" / "catalogue_summary.svg",
        output_dir / "visualize" / "catalogue_summary.pdf",
        output_dir / "visualize" / "assembly_vs_read.svg",
        output_dir / "visualize" / "assembly_vs_read.pdf",
        output_dir / "visualize" / "in_silico_fish.svg",
        output_dir / "visualize" / "in_silico_fish.pdf",
    ]
    for path in expected:
        assert path.is_file(), path

    total_size = sum(path.stat().st_size for path in output_dir.rglob("*") if path.is_file())
    assert total_size < 1_000_000
