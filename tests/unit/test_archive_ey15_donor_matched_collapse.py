import csv
import json
from pathlib import Path

import pytest

from benchmarks.scripts.archive_ey15_donor_matched_collapse import (
    validate_evaluation,
    validate_figure,
    validate_resources,
)
from benchmarks.challenge.schema import digest_file


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_validate_evaluation_requires_all_fates_and_independent_pass(tmp_path: Path) -> None:
    directory = tmp_path / "evaluation"
    directory.mkdir()
    rows = [
        {"family_id": "a", "eligibility": "eligible"},
        {"family_id": "b", "eligibility": "not_source_eligible"},
        {"family_id": "c", "eligibility": "technical_failure"},
    ]
    write_tsv(directory / "family_metrics.tsv", rows)
    (directory / "summary.json").write_text(
        json.dumps(
            {
                "interpretation": (
                    "donor_matched_retrospective_reference_proxy_not_absolute_biological_truth"
                ),
                "all_family_rows": 3,
                "eligible_family_rows": 1,
                "not_source_eligible_rows": 1,
                "technical_failure_rows": 1,
            }
        )
    )
    (directory / "independent_verification.json").write_text(
        json.dumps({"verification_passed": True, "failures": []})
    )
    summary = validate_evaluation(tmp_path, "evaluation")
    assert summary["technical_failure_rows"] == 1
    (directory / "independent_verification.json").write_text(
        json.dumps({"verification_passed": False, "failures": ["changed"]})
    )
    with pytest.raises(ValueError, match="verification failed"):
        validate_evaluation(tmp_path, "evaluation")


def test_validate_resources_rejects_incomplete_receipt(tmp_path: Path) -> None:
    for label in ("assembly_stage_resources", "read_stage_resources"):
        directory = tmp_path / label
        directory.mkdir()
        (directory / "receipt.json").write_text(
            json.dumps(
                {"complete": True, "completed_stage_count": 1, "requested_stage_count": 1}
            )
        )
        write_tsv(directory / "stages.tsv", [{"stage": "one", "exit_code": 0}])
    assert len(validate_resources(tmp_path)) == 2
    (tmp_path / "read_stage_resources/receipt.json").write_text(
        json.dumps(
            {"complete": False, "completed_stage_count": 1, "requested_stage_count": 2}
        )
    )
    with pytest.raises(ValueError, match="profile incomplete"):
        validate_resources(tmp_path)


def test_validate_figure_requires_editable_svg_and_matching_hashes(
    tmp_path: Path,
) -> None:
    provenances = {}
    for version in ("figures_v1", "figures_v2"):
        directory = tmp_path / version
        directory.mkdir()
        svg = directory / "ey15_donor_validation.svg"
        png = directory / "ey15_donor_validation.png"
        panel_source = directory / "panel_source.tsv"
        svg.write_text(f"<svg><text>{version}</text></svg>")
        png.write_text(version)
        panel_source.write_text("panel\tvalue\nA\t1\n")
        provenance = {
            "complete": True,
            "panel_count": 6,
            "primary_family_count": 19,
            "svg_raster_image_element_count": 0,
            "svg_text_element_count": 1,
            "outputs": {
                path.name: {
                    "bytes": path.stat().st_size,
                    "sha256": digest_file(path),
                }
                for path in (svg, png, panel_source)
            },
        }
        (directory / "figure_provenance.json").write_text(json.dumps(provenance))
        provenances[version] = provenance
    assert validate_figure(tmp_path)["qa_history"]["panel_source_identical"] is True
    directory = tmp_path / "figures_v2"
    provenance = provenances["figures_v2"]
    provenance["svg_raster_image_element_count"] = 1
    (directory / "figure_provenance.json").write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match="raster image"):
        validate_figure(tmp_path)
