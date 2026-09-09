from __future__ import annotations

import json
from pathlib import Path

from benchmarks.scripts.evaluate_orthogonal_illumina_abundance import (
    evaluate,
    prepare_targets,
)


def test_prepare_and_evaluate_frozen_targets(tmp_path: Path) -> None:
    monomer = "ACGTTGCACTGATCGAACCTG"
    catalogue = tmp_path / "monomers.fa"
    controls = tmp_path / "controls.tsv"
    target_fasta = tmp_path / "targets.fa"
    target_map = tmp_path / "targets.tsv"
    prepare_receipt = tmp_path / "prepare.json"
    catalogue.write_text(
        f">family_id=F1;monomer_id=M1;length_bp={len(monomer)};confidence=high\n{monomer}\n"
    )
    controls.write_text("kmer\texpected_copy_number\nAACCGGTTAACCGGT\t1\n")
    prepared = prepare_targets(
        catalogue, controls, target_fasta, target_map, prepare_receipt, k=15
    )
    assert prepared["single_copy_controls"] == 1
    assert prepared["diagnostic_kmers"] > 0

    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "preorthogonal_hifi_newer_assembly_deficit_candidates": {
                    "ey15": ["F1"],
                    "macadamia": [],
                },
                "binary_threshold": 0.6,
                "quantification_bias_fold_threshold": 1.5,
            }
        )
    )
    prior = tmp_path / "prior.tsv"
    prior.write_text(
        "family_id\tread_estimated_bp\tnew_assembly_bp\teligibility\n"
        "F1\t2100\t100\teligible\n"
    )
    kmc_dump = tmp_path / "counts.txt"
    target_rows = [line.split("\t") for line in target_map.read_text().splitlines()[1:]]
    with kmc_dump.open("w") as handle:
        for word, target_type, *_ in target_rows:
            handle.write(f"{word}\t{10 if target_type == 'single_copy_control' else 1000}\n")
    output = tmp_path / "result.tsv"
    summary_path = tmp_path / "summary.json"
    summary = evaluate(
        "ey15",
        config,
        catalogue,
        controls,
        kmc_dump,
        prior,
        output,
        summary_path,
        k=15,
    )
    row = output.read_text().splitlines()[1].split("\t")
    header = output.read_text().splitlines()[0].split("\t")
    result = dict(zip(header, row))
    assert float(result["orthogonal_estimated_bp"]) > 100
    assert result["k_result"] == "orthogonal_supports_residual_collapse_at_this_k"
    assert summary["control_mean_depth"] == 10.0
