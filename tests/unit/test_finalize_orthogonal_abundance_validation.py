from __future__ import annotations

import csv
from pathlib import Path

import pytest

from benchmarks.scripts.finalize_orthogonal_abundance_validation import finalize


HEADER = [
    "species_key",
    "family_id",
    "orthogonal_estimated_bp",
    "frozen_hifi_estimated_bp",
    "newer_assembly_bp",
    "newer_assembly_orthogonal_ratio",
    "hifi_orthogonal_ratio",
    "estimated_under_representation_bp",
    "preorthogonal_candidate",
    "eligibility",
    "k_result",
]


def write_metric(path: Path, species: str, state: str, orthogonal_bp: int) -> None:
    values = [
        species,
        "F1",
        orthogonal_bp,
        2000,
        100,
        100 / orthogonal_bp,
        2000 / orthogonal_bp,
        orthogonal_bp - 100,
        "true",
        "eligible",
        state,
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerow(values)


def test_ey15_cross_k_disagreement_is_unresolved(tmp_path: Path) -> None:
    k21, k31 = tmp_path / "k21.tsv", tmp_path / "k31.tsv"
    write_metric(k21, "ey15", "orthogonal_supports_residual_collapse_at_this_k", 1000)
    write_metric(k31, "ey15", "orthogonal_supports_quantification_bias_at_this_k", 120)
    summary = finalize(
        "ey15", k21, k31, tmp_path / "out.tsv", tmp_path / "summary.json"
    )
    assert summary["final_interpretation_counts"] == {"unresolved": 1}


@pytest.mark.parametrize(
    ("ont_state", "expected"),
    [
        ("ONT_direction_supports_residual_collapse", "orthogonal_supports_residual_collapse"),
        ("ONT_direction_does_not_support_residual_collapse", "unresolved"),
    ],
)
def test_macadamia_residual_requires_ont_direction(
    tmp_path: Path, ont_state: str, expected: str
) -> None:
    k21, k31 = tmp_path / "k21.tsv", tmp_path / "k31.tsv"
    state = "orthogonal_supports_residual_collapse_at_this_k"
    write_metric(k21, "macadamia", state, 1000)
    write_metric(k31, "macadamia", state, 900)
    ont = tmp_path / "ont.tsv"
    ont.write_text(
        "family_id\tpreorthogonal_candidate\tONT_direction_result\n"
        f"F1\ttrue\t{ont_state}\n"
    )
    summary = finalize(
        "macadamia",
        k21,
        k31,
        tmp_path / "out.tsv",
        tmp_path / "summary.json",
        ont_direction=ont,
    )
    assert summary["final_interpretation_counts"] == {expected: 1}
