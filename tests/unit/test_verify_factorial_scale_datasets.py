from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from benchmarks.scripts.generate_factorial_scale import run as generate
from benchmarks.scripts.verify_factorial_scale_datasets import run as verify


def test_verify_factorial_scale_datasets_checks_every_manifest_payload(tmp_path: Path) -> None:
    histogram = tmp_path / "histogram.tsv"
    histogram.write_text("length_bp\tread_count\n50\t1\n")
    generation_config = {
        "seeds": {"development": [11], "validation": [12], "heldout": [13]},
        "genome_bp": 1000,
        "background_gc": 0.4,
        "periods": [11],
        "copies": [3],
        "gc_fractions": [0.5],
        "unit_substitution_rates": [0],
        "coverages": [1],
        "read_error_models": [
            {
                "label": "clean",
                "substitution_rate": 0,
                "insertion_rate": 0,
                "deletion_rate": 0,
            }
        ],
    }
    config_path = tmp_path / "generation.json"
    config_path.write_text(json.dumps(generation_config))
    dataset = tmp_path / "dataset"
    generate(config_path, histogram, dataset, 12, 100_000, "validation")
    receipt = json.loads((dataset / "generation_receipt.json").read_text())
    validation_config = tmp_path / "validation.yaml"
    validation_config.write_text(
        yaml.safe_dump(
            {
                "seeds": [12],
                "split": "validation",
                "condition_count_per_seed": 1,
                "dataset_generation": {
                    "config_sha256": receipt["config_sha256"],
                    "length_histogram_sha256": receipt["histogram_sha256"],
                    "source_sha256": receipt["source_sha256"],
                },
            }
        )
    )
    output = tmp_path / "verification.json"
    result = verify([dataset], validation_config, output)
    assert result["complete"] is True
    assert result["manifest_count"] == 2
    assert result["payload_file_count"] == 7
    assert json.loads(output.read_text()) == result

    reads = dataset / "reads/condition_001/reads.fa"
    reads.write_text(reads.read_text() + "A\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify([dataset], validation_config, None)
