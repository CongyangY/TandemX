import csv
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_macadamia_jansenii_donor_matched_collapse import (
    validate_alignment,
    validate_evaluation,
    validate_resources,
)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_validate_evaluation_preserves_all_fates(tmp_path: Path) -> None:
    directory = tmp_path / "evaluation"
    directory.mkdir()
    write_tsv(
        directory / "family_metrics.tsv",
        [
            {"family_id": "a", "eligibility": "eligible"},
            {"family_id": "b", "eligibility": "not_source_eligible"},
            {"family_id": "c", "eligibility": "technical_failure"},
        ],
    )
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
    assert validate_evaluation(tmp_path, "evaluation")["eligible_family_rows"] == 1
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


def test_validate_alignment_requires_independent_pass(tmp_path: Path) -> None:
    run = tmp_path / "old_new_alignment_audit"
    context = tmp_path / "old_new_alignment_context"
    run.mkdir()
    context.mkdir()
    (run / "receipt.json").write_text(json.dumps({"complete": True}))
    (context / "summary.json").write_text(json.dumps({"eligible_family_count": 1}))
    write_tsv(context / "family_alignment_context.tsv", [{"family_id": "a"}])
    verification = context / "independent_verification.json"
    verification.write_text(json.dumps({"verification_passed": True, "failures": []}))
    assert validate_alignment(tmp_path)["eligible_family_count"] == 1
    verification.write_text(
        json.dumps({"verification_passed": False, "failures": ["changed"]})
    )
    with pytest.raises(ValueError, match="verification failed"):
        validate_alignment(tmp_path)


def test_committed_archive_is_complete_and_hash_checked() -> None:
    root = (
        Path(__file__).parents[2]
        / "paper/evidence/macadamia_jansenii_donor_matched_collapse_v1"
    )
    manifest = json.loads((root / "archive_manifest.json").read_text())
    assert manifest["complete"] is True
    assert len(manifest["files"]) == 40
    assert all(
        digest_file(root / row["file"]) == row["sha256"]
        for row in manifest["files"]
    )
    headline = json.loads((root / "headline_summary.json").read_text())
    assert headline["complete"] is True
    assert headline["primary_total_bases_depth"]["confusion"] == {
        "TP": 2,
        "FN": 0,
        "FP": 0,
        "TN": 41,
    }
    assert headline["sensitivity_reported_28gb_depth"]["confusion"] == {
        "TP": 2,
        "FN": 0,
        "FP": 0,
        "TN": 41,
    }
