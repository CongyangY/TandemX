from __future__ import annotations

import random
from pathlib import Path

import yaml

from benchmarks.scripts.evaluate_quantify_calibration import (
    METHODS,
    select_single_copy_controls,
    summarize,
)
from tandemx.utils.kmers import canonical_kmer


def test_control_selection_is_unique_and_excludes_array(tmp_path: Path) -> None:
    rng = random.Random(17)
    background = "".join(rng.choices("ACGT", k=500))
    monomer = "ACGTTGCACTGATCGAACCTG"
    genome = background[:200] + monomer * 4 + background[200:]
    genome_path = tmp_path / "genome.fa"
    truth_path = tmp_path / "truth.tsv"
    catalogue_path = tmp_path / "catalogue.fa"
    genome_path.write_text(f">chr\n{genome}\n")
    truth_path.write_text(
        "start\tend\n200\t284\n"
    )
    catalogue_path.write_text(f">f1\n{monomer}\n")
    controls, receipt = select_single_copy_controls(
        genome_path,
        truth_path,
        catalogue_path,
        k=15,
        desired=10,
        stride=7,
        candidate_multiplier=5,
    )
    canonical_genome = [
        canonical_kmer(genome[index : index + 15])
        for index in range(len(genome) - 14)
    ]
    assert len(controls) == len(set(controls)) == 10
    assert all(canonical_genome.count(word) == 1 for word in controls)
    assert receipt["selected_controls"] == 10


def test_calibration_summary_keeps_accuracy_and_resource_grains() -> None:
    metrics = [
        {
            "seed": seed,
            "condition_id": "c1",
            "method": "baseline_total_bases",
            "coverage": 5,
            "error_model": "iid",
            "signed_relative_error": error,
            "absolute_relative_error": abs(error),
            "estimator_minus_sampling_oracle": error / 2,
            "interval_contains_truth": seed == 1,
            "interval_relative_width": 0.4,
        }
        for seed, error in ((1, -0.2), (2, 0.1))
    ]
    executions = [
        {
            "seed": seed,
            "condition_id": "c1",
            "method": "baseline_total_bases",
            "runtime_seconds": runtime,
            "peak_rss_mib": rss,
        }
        for seed, runtime, rss in ((1, 2.0, 10.0), (2, 4.0, 14.0))
    ]
    row = summarize(metrics, executions)[0]
    assert row["independent_genomes"] == 2
    assert row["family_conditions"] == 2
    assert row["mean_signed_relative_error"] == -0.05
    assert row["diagnostic_spread_truth_coverage"] == 0.5
    assert row["median_runtime_seconds"] == 3.0
    assert row["median_peak_rss_mib"] == 12.0


def test_frozen_calibration_config_matches_implementation() -> None:
    path = Path("benchmarks/configs/quantify_calibration_development_v1.yaml")
    config = yaml.safe_load(path.read_text())
    assert config["split"] == "development"
    assert tuple(config["methods"]) == METHODS
    assert len(config["seeds"]) == len(set(config["seeds"])) == 3
