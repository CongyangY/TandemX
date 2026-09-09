#!/usr/bin/env python3
"""Validate and archive compact Macadamia donor-matched collapse evidence."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any

from benchmarks.challenge.schema import digest_file, write_table


RESULT_FILES = (
    "assembly_stage_resources/receipt.json",
    "assembly_stage_resources/stages.tsv",
    "read_stage_resources/receipt.json",
    "read_stage_resources/stages.tsv",
    "run/discover/discovery_summary.json",
    "run/discover/family_audit_summary.json",
    "run/discover/run_config.yaml",
    "run/discover/monomers.fa",
    "run/discover/families.tsv",
    "run/discover/monomer_membership.tsv",
    "run/quantify_primary/run_config.yaml",
    "run/quantify_primary/copy_number.tsv",
    "run/quantify_reported_depth_sensitivity/run_config.yaml",
    "run/quantify_reported_depth_sensitivity/copy_number.tsv",
    "run/locate_old/run_config.yaml",
    "run/locate_old/arrays.bed",
    "run/locate_new/run_config.yaml",
    "run/locate_new/arrays.bed",
    "evaluation_primary_total_bases/summary.json",
    "evaluation_primary_total_bases/family_metrics.tsv",
    "evaluation_primary_total_bases/independent_verification.json",
    "evaluation_reported_depth_sensitivity/summary.json",
    "evaluation_reported_depth_sensitivity/family_metrics.tsv",
    "evaluation_reported_depth_sensitivity/independent_verification.json",
    "old_new_alignment_audit/environment.json",
    "old_new_alignment_audit/old_new.paf",
    "old_new_alignment_audit/receipt.json",
    "old_new_alignment_audit/resources/stages.tsv",
    "old_new_alignment_context/summary.json",
    "old_new_alignment_context/family_alignment_context.tsv",
    "old_new_alignment_context/independent_verification.json",
)


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_evaluation(source: Path, name: str) -> dict[str, Any]:
    directory = source / name
    summary = load_json(directory / "summary.json")
    verification = load_json(directory / "independent_verification.json")
    rows = read_tsv(directory / "family_metrics.tsv")
    if verification.get("verification_passed") is not True:
        raise ValueError(f"independent evaluation verification failed: {directory}")
    if verification.get("failures") != []:
        raise ValueError(f"independent evaluation retained failures: {directory}")
    if len(rows) != summary.get("all_family_rows"):
        raise ValueError(f"family row count differs from summary: {directory}")
    if summary.get("interpretation") != (
        "donor_matched_retrospective_reference_proxy_not_absolute_biological_truth"
    ):
        raise ValueError(f"evaluation lost its evidence boundary: {directory}")
    expected = {
        "eligible": "eligible_family_rows",
        "not_source_eligible": "not_source_eligible_rows",
        "technical_failure": "technical_failure_rows",
    }
    for fate, field in expected.items():
        if sum(row["eligibility"] == fate for row in rows) != summary.get(field):
            raise ValueError(f"{fate} row count differs from summary: {directory}")
    return summary


def validate_resources(source: Path) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    for label in ("assembly_stage_resources", "read_stage_resources"):
        directory = source / label
        receipt = load_json(directory / "receipt.json")
        rows = read_tsv(directory / "stages.tsv")
        if receipt.get("complete") is not True:
            raise ValueError(f"resource profile incomplete: {directory}")
        if receipt.get("completed_stage_count") != receipt.get("requested_stage_count"):
            raise ValueError(f"resource stage count incomplete: {directory}")
        if len(rows) != receipt.get("completed_stage_count"):
            raise ValueError(f"resource table count differs from receipt: {directory}")
        if any(int(row["exit_code"]) != 0 for row in rows):
            raise ValueError(f"resource table retains failed completed stage: {directory}")
        result[label] = rows
    return result


def validate_alignment(source: Path) -> dict[str, Any]:
    run = load_json(source / "old_new_alignment_audit/receipt.json")
    summary = load_json(source / "old_new_alignment_context/summary.json")
    verification = load_json(
        source / "old_new_alignment_context/independent_verification.json"
    )
    rows = read_tsv(source / "old_new_alignment_context/family_alignment_context.tsv")
    if run.get("complete") is not True:
        raise ValueError("old/new assembly alignment did not complete")
    if verification.get("verification_passed") is not True:
        raise ValueError("independent old/new alignment verification failed")
    if verification.get("failures") != []:
        raise ValueError("independent old/new alignment verification retains failures")
    if len(rows) != summary.get("eligible_family_count"):
        raise ValueError("old/new alignment family row count differs from summary")
    return summary


def archive(
    source: Path,
    enrollment: Path,
    configs: list[Path],
    outdir: Path,
) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite archive: {outdir}")
    source = source.resolve()
    missing = [name for name in RESULT_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"required Macadamia evidence is incomplete: {missing}")
    enrollment_payload = load_json(enrollment)
    if enrollment_payload.get("complete") is not True:
        raise ValueError("assembly enrollment is incomplete")
    primary = validate_evaluation(source, "evaluation_primary_total_bases")
    sensitivity = validate_evaluation(
        source, "evaluation_reported_depth_sensitivity"
    )
    resources = validate_resources(source)
    alignment = validate_alignment(source)
    if primary["all_family_rows"] != sensitivity["all_family_rows"]:
        raise ValueError("primary and sensitivity family universes differ")
    if primary["eligible_family_rows"] != sensitivity["eligible_family_rows"]:
        raise ValueError("primary and sensitivity eligibility denominators differ")

    outdir.mkdir(parents=True)
    manifest_rows: list[dict[str, Any]] = []
    inputs = {"assembly_enrollment_receipt.json": enrollment}
    inputs.update({path.name: path for path in configs})
    if len(inputs) != len(configs) + 1:
        raise ValueError("config basenames are not unique")
    for destination_name, source_path in inputs.items():
        target = outdir / destination_name
        shutil.copyfile(source_path, target)
        manifest_rows.append(_manifest_row(target, outdir, str(source_path.resolve())))
    for name in RESULT_FILES:
        source_path = source / name
        target = outdir / "results" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        if digest_file(target) != digest_file(source_path):
            raise OSError(f"archive copy hash differs: {source_path}")
        manifest_rows.append(_manifest_row(target, outdir, str(source_path)))

    resource_rows = [
        {"resource_group": group, **row}
        for group, rows in resources.items()
        for row in rows
    ]
    write_table(outdir / "resource_summary.tsv", resource_rows, list(resource_rows[0]))
    headline = {
        "schema_version": 1,
        "complete": True,
        "interpretation": (
            "second_species_same_sample_reference_proxy;"
            "not_same_extraction_or_absolute_biological_truth"
        ),
        "primary_total_bases_depth": primary,
        "sensitivity_reported_28gb_depth": sensitivity,
        "old_new_alignment_context": alignment,
        "resource_stage_count": len(resource_rows),
        "warning": (
            "paper_states_same_sample_but_archival_identifiers_and_collection_dates_differ;"
            "new_assembly_shares_HiFi_evidence_with_read_estimator;"
            "alignment_context_is_not_independent_copy_truth"
        ),
    }
    (outdir / "headline_summary.json").write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (outdir / "README.md").write_text(
        "# Macadamia jansenii donor-matched collapse validation v1\n\n"
        "This compact archive retains frozen configurations, exact assembly "
        "enrollment, discovery, two complete-read quantifications, old/new "
        "localizations, all family fates, independent recomputations, resource "
        "profiles and the post-primary old/new alignment audit.\n\n"
        "The update paper describes the HiFi material as the same sample used "
        "for the earlier CLR assembly, but archival identifiers and collection "
        "dates differ. The comparison is therefore a same-sample reference "
        "proxy, not proof of the same DNA extraction or absolute truth. The "
        "newer assembly also shares HiFi evidence with the read estimator.\n",
        encoding="utf-8",
    )
    for name in ("resource_summary.tsv", "headline_summary.json", "README.md"):
        manifest_rows.append(
            _manifest_row(
                outdir / name,
                outdir,
                "generated_by_archive_macadamia_jansenii_donor_matched_collapse",
            )
        )
    manifest = {"schema_version": 1, "complete": True, "files": manifest_rows}
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _manifest_row(path: Path, root: Path, source: str) -> dict[str, Any]:
    return {
        "file": path.relative_to(root).as_posix(),
        "source": source,
        "bytes": path.stat().st_size,
        "sha256": digest_file(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--enrollment", required=True, type=Path)
    parser.add_argument("--config", required=True, action="append", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.enrollment, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
