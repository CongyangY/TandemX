from __future__ import annotations

import json
from pathlib import Path

from benchmarks.scripts.evaluate_macadamia_ont_occupancy import (
    calibrate_hifi_mapping,
    evaluate_ont,
    prepare_templates,
    summarize_paf,
)


def test_direction_only_occupancy_flow(tmp_path: Path) -> None:
    catalogue = tmp_path / "monomers.fa"
    prior = tmp_path / "prior.tsv"
    template = tmp_path / "templates.fa"
    catalogue.write_text(
        ">family_id=F1;monomer_id=M1;length_bp=20;confidence=high\nACGTTGCACTGATCGAACCT\n"
        ">family_id=F2;monomer_id=M2;length_bp=20;confidence=high\nGGTACCATGCTAGCTTACGA\n"
    )
    prior.write_text(
        "family_id\teligibility\tread_estimated_bp\tnew_assembly_bp\n"
        "F1\teligible\t2000\t100\n"
        "F2\teligible\t1000\t900\n"
    )
    receipt = tmp_path / "template.json"
    prepare_templates(catalogue, prior, template, receipt, template_length=1000)

    paf = tmp_path / "reads.paf"
    paf.write_text(
        "r1\t1000\t0\t1000\t+\tF1\t1000\t0\t1000\t900\t1000\t60\ttp:A:P\n"
        "r2\t1000\t0\t500\t+\tF1\t1000\t0\t500\t450\t500\t60\ttp:A:P\n"
        "r2\t1000\t500\t1000\t+\tF2\t1000\t0\t500\t450\t500\t60\ttp:A:P\n"
    )
    occupancy = tmp_path / "hifi.tsv"
    summarize_paf(
        paf,
        template,
        occupancy,
        tmp_path / "hifi_summary.json",
        total_library_bases=10_000,
        genome_size_bp=10_000,
        minimum_block_bp=500,
        minimum_identity=0.75,
    )
    calibration = tmp_path / "calibration.json"
    result = calibrate_hifi_mapping(occupancy, prior, calibration)
    assert result["positive_eligible_families"] == 1
    assert result["global_median_raw_HiFi_mapping_over_frozen_HiFi_kmer"] == 0.5

    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "macadamia_genome_size_sensitivity_bp": [8000, 10000],
                "binary_threshold": 0.6,
                "preorthogonal_hifi_newer_assembly_deficit_candidates": {
                    "macadamia": ["F1"]
                },
            }
        )
    )
    output = tmp_path / "ont.tsv"
    summary = evaluate_ont(
        occupancy,
        calibration,
        prior,
        config,
        output,
        tmp_path / "ont_summary.json",
    )
    assert summary["candidate_result_counts"] == {
        "ONT_direction_supports_residual_collapse": 1
    }
    assert "negative_context_no_residual_deficit" in output.read_text()
