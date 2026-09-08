import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.tidecluster.factorial_run import build_stage_specs, validate_datasets


def _config(tmp_path: Path) -> dict:
    genome = tmp_path / "genome"
    genome.mkdir()
    (genome / "genome.fa").write_text(">chr\nACGT\n")
    (genome / "catalogue.fa").write_text(">f1\nACGT\n")
    (genome / "truth_copy_number.tsv").write_text("header\n")
    files = {name: digest_file(genome / name) for name in ("genome.fa", "catalogue.fa", "truth_copy_number.tsv")}
    (genome / "manifest.json").write_text(
        json.dumps({"generator": "streamed_factorial_genome_v1", "seed": 7, "genome_bp": 10_000_000, "files": files})
    )
    return {
        "image": "image:version",
        "cpus": 4,
        "minimum_length": 100,
        "settings": {"default_primary": {"tidehunter_arguments": "-p 40 -P 3000 -c 5 -e 0.25"}},
        "datasets": [
            {
                "seed": 7,
                "genome_dir": str(genome),
                "genome_sha256": files["genome.fa"],
                "catalogue_sha256": files["catalogue.fa"],
                "truth_sha256": files["truth_copy_number.tsv"],
            }
        ],
    }


def test_factorial_stage_manifest_mounts_only_observable_assembly(tmp_path: Path) -> None:
    config = _config(tmp_path)
    datasets = validate_datasets(config)
    stages, runs = build_stage_specs(config, datasets, tmp_path / "out", "docker")
    assert len(stages) == 2 and len(runs) == 1
    commands = "\n".join(" ".join(stage["command"]) for stage in stages)
    assert "genome.fa:/input/genome.fa:ro" in commands
    assert "truth_copy_number.tsv" not in commands
    assert "catalogue.fa" not in commands
    assert "-p 40 -P 3000 -c 5 -e 0.25" in commands


def test_factorial_stage_manifest_rejects_changed_input(tmp_path: Path) -> None:
    config = _config(tmp_path)
    Path(config["datasets"][0]["genome_dir"], "genome.fa").write_text("changed")
    with pytest.raises(ValueError, match="differs"):
        validate_datasets(config)
