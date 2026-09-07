from __future__ import annotations

import json

import pytest
import yaml

from benchmarks.challenge.schema import digest_file, read_table
from benchmarks.scripts.evaluate_quantify_depth_gated_validation import (
    run as validate,
    worker,
)
from benchmarks.scripts.generate_factorial_scale import run as generate


def test_depth_gated_validation_runs_public_commands_on_toy_data(tmp_path) -> None:
    histogram = tmp_path / "lengths.tsv"
    histogram.write_text("length_bp\tread_count\n200\t1\n")
    generator_config = {
        "seeds": {
            "development": [],
            "validation": [6461],
            "heldout": [7461],
        },
        "genome_bp": 10_000,
        "background_gc": 0.45,
        "periods": [31],
        "copies": [3, 7],
        "gc_fractions": [0.5],
        "unit_substitution_rates": [0, 0.05],
        "coverages": [3],
        "read_error_models": [
            {
                "label": "clean",
                "substitution_rate": 0,
                "insertion_rate": 0,
                "deletion_rate": 0,
            }
        ],
    }
    generator_path = tmp_path / "generator.json"
    generator_path.write_text(json.dumps(generator_config))
    dataset = tmp_path / "dataset"
    generate(generator_path, histogram, dataset, 6461, 100_000, "validation")

    validation_config = {
        "benchmark_id": "toy_depth_gate",
        "split": "validation",
        "seeds": [6461],
        "forbidden_development_seeds": [6361],
        "condition_count_per_seed": 1,
        "dataset_generation": {
            "config_sha256": digest_file(generator_path),
            "length_histogram_sha256": digest_file(histogram),
            "source_sha256": json.loads(
                (dataset / "generation_receipt.json").read_text()
            )["source_sha256"],
        },
        "k": 11,
        "control_kmer_count": 5,
        "control_candidate_multiplier": 10,
        "control_stride_bp": 97,
        "timeout_seconds": 30,
        "methods": ["baseline_total_bases", "empirical_controls"],
        "candidate": {
            "rule": "empirical_controls_when_control_mean_depth_at_least_threshold_else_total_bases",
            "control_mean_depth_threshold": 0,
        },
        "acceptance_gates": {
            "failed_executions_max": 0,
            "aggregate_mare_reduction_min": -100,
            "candidate_mare_max": 100,
            "minimum_seed_mare_reduction_strictly_greater_than": -100,
            "minimum_coverage_mare_reduction_min": -100,
            "paired_improved_fraction_min": 0,
            "paired_nonworse_fraction_min": 0,
            "selected_condition_count_per_branch_min": 0,
        },
    }
    validation_path = tmp_path / "validation.yaml"
    validation_path.write_text(yaml.safe_dump(validation_config, sort_keys=False))
    outdir = tmp_path / "result"
    validate(validation_path, [dataset], outdir)

    validation = json.loads((outdir / "validation.json").read_text())
    gates = json.loads((outdir / "gate_results.json").read_text())
    assert validation["complete"] and validation["executions"] == 2
    assert validation["raw_family_conditions"] == 8
    assert validation["candidate_family_conditions"] == 4
    assert gates["status"] == "passed"
    assert {row["method"] for row in read_table(outdir / "metrics.tsv")} == {
        "baseline_total_bases",
        "empirical_controls",
        "depth_gated_controls",
    }
    for method in ("baseline_total_bases", "empirical_controls"):
        output = outdir / "runs" / "s6461" / "condition_001" / method / "output"
        assert (output / "run_config.yaml").is_file()
        assert (output / "run.log").is_file()

    receipt_path = dataset / "generation_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["config_sha256"] = "0" * 64
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="wrong split semantics"):
        worker(validation_path, [dataset], tmp_path / "wrong-generation")
