#!/usr/bin/env python3
"""Validate and archive compact Ey15 donor-matched collapse evidence."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
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
    "run/quantify_depth107_sensitivity/run_config.yaml",
    "run/quantify_depth107_sensitivity/copy_number.tsv",
    "run/locate_old/run_config.yaml",
    "run/locate_old/arrays.bed",
    "run/locate_new/run_config.yaml",
    "run/locate_new/arrays.bed",
    "run/locate_sensitivity/run_config.yaml",
    "run/locate_sensitivity/arrays.bed",
    "evaluation_primary_total_bases/summary.json",
    "evaluation_primary_total_bases/family_metrics.tsv",
    "evaluation_primary_total_bases/independent_verification.json",
    "evaluation_depth107_sensitivity/summary.json",
    "evaluation_depth107_sensitivity/family_metrics.tsv",
    "evaluation_depth107_sensitivity/independent_verification.json",
    "author_annotation_audit/summary.json",
    "author_annotation_audit/annotation_class_summary.tsv",
    "author_annotation_audit/family_annotation_summary.tsv",
    "author_annotation_audit/independent_raster_validation.json",
    "reference_annotation_context_posthoc_v2/summary.json",
    "reference_annotation_context_posthoc_v2/reference_annotation_rows.tsv",
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
    eligible = [row for row in rows if row["eligibility"] == "eligible"]
    source_ineligible = [row for row in rows if row["eligibility"] == "not_source_eligible"]
    technical = [row for row in rows if row["eligibility"] == "technical_failure"]
    if len(eligible) != summary.get("eligible_family_rows"):
        raise ValueError(f"eligible row count differs from summary: {directory}")
    if len(source_ineligible) != summary.get("not_source_eligible_rows"):
        raise ValueError(f"source-ineligible row count differs from summary: {directory}")
    if len(technical) != summary.get("technical_failure_rows"):
        raise ValueError(f"technical-failure row count differs from summary: {directory}")
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


def validate_context(source: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    author = load_json(source / "author_annotation_audit/summary.json")
    raster = load_json(source / "author_annotation_audit/independent_raster_validation.json")
    context = load_json(source / "reference_annotation_context_posthoc_v2/summary.json")
    context_rows = read_tsv(
        source / "reference_annotation_context_posthoc_v2/reference_annotation_rows.tsv"
    )
    if not all(raster.get("checks", {}).values()):
        raise ValueError("independent author-annotation raster validation failed")
    if context.get("analysis_status") != (
        "posthoc_descriptive_audit_after_reference_state_inspection"
    ):
        raise ValueError("annotation association is not labelled post hoc")
    if context.get("primary_benchmark_changed") is not False:
        raise ValueError("annotation context claims the primary benchmark changed")
    if len(context_rows) != context.get("eligible_families"):
        raise ValueError("annotation context row count differs from summary")
    return author, context


def archive(
    source: Path,
    enrollment: Path,
    preregistration: Path,
    read_manifest: Path,
    assembly_manifest: Path,
    outdir: Path,
) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite archive: {outdir}")
    source = source.resolve()
    missing = [name for name in RESULT_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"required Ey15 evidence is incomplete: {missing}")
    enrollment_payload = load_json(enrollment)
    if enrollment_payload.get("complete") is not True:
        raise ValueError("assembly enrollment is incomplete")
    primary = validate_evaluation(source, "evaluation_primary_total_bases")
    depth107 = validate_evaluation(source, "evaluation_depth107_sensitivity")
    resources = validate_resources(source)
    author, context = validate_context(source)
    if primary["all_family_rows"] != depth107["all_family_rows"]:
        raise ValueError("primary and depth-107 family universes differ")
    if primary["eligible_family_rows"] != depth107["eligible_family_rows"]:
        raise ValueError("primary and depth-107 eligibility denominators differ")

    outdir.mkdir(parents=True)
    inputs = {
        "preregistration.json": preregistration,
        "read_stage_manifest.json": read_manifest,
        "assembly_stage_manifest.json": assembly_manifest,
        "assembly_enrollment_receipt.json": enrollment,
    }
    manifest_rows: list[dict[str, Any]] = []
    for destination_name, source_path in inputs.items():
        target = outdir / destination_name
        shutil.copyfile(source_path, target)
        manifest_rows.append(
            {
                "file": destination_name,
                "source": str(source_path.resolve()),
                "bytes": target.stat().st_size,
                "sha256": digest_file(target),
            }
        )
    for name in RESULT_FILES:
        source_path = source / name
        target = outdir / "results" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        if target.stat().st_size != source_path.stat().st_size:
            raise OSError(f"archive copy size differs: {source_path}")
        if digest_file(target) != digest_file(source_path):
            raise OSError(f"archive copy hash differs: {source_path}")
        manifest_rows.append(
            {
                "file": target.relative_to(outdir).as_posix(),
                "source": str(source_path),
                "bytes": target.stat().st_size,
                "sha256": digest_file(target),
            }
        )

    resource_rows = []
    for group, rows in resources.items():
        for row in rows:
            resource_rows.append({"resource_group": group, **row})
    write_table(
        outdir / "resource_summary.tsv",
        resource_rows,
        list(resource_rows[0]),
    )
    headline = {
        "schema_version": 1,
        "complete": True,
        "interpretation": (
            "single_donor_matched_retrospective_high_quality_reference_proxy;"
            "not_absolute_or_fully_independent_biological_truth"
        ),
        "primary_total_bases_depth": primary,
        "sensitivity_explicit_depth_107": depth107,
        "author_annotation_audit": author,
        "posthoc_reference_annotation_context": context,
        "resource_stage_count": len(resource_rows),
        "warning": (
            "new_assembly_shares_HiFi_evidence_with_read_estimator;"
            "single_material;low_author_centromere_annotation_recall;"
            "posthoc_annotation_association_is_not_primary_validation"
        ),
    }
    (outdir / "headline_summary.json").write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (outdir / "README.md").write_text(
        "# Ey15-2 donor-matched collapse validation v1\n\n"
        "This compact archive contains the frozen preregistration, exact "
        "assembly-enrollment receipt, discovery catalogue, old/new/sensitivity "
        "localizations, two full-read quantifications, all family fates, "
        "independent standard-library recomputations, resource profiles and "
        "author-annotation audits. `headline_summary.json` is the concise "
        "result; `results/` retains the auditable inputs to each claim.\n\n"
        "The newer assembly is a donor-matched high-quality reference proxy, "
        "not absolute truth, and shares HiFi evidence with the read estimator. "
        "The selected author annotations cover only a fraction of the reported "
        "centromeric sequence. The annotation association was defined after "
        "reference-state inspection and is explicitly post hoc.\n",
        encoding="utf-8",
    )
    for name in ("resource_summary.tsv", "headline_summary.json", "README.md"):
        path = outdir / name
        manifest_rows.append(
            {
                "file": name,
                "source": "generated_by_archive_ey15_donor_matched_collapse",
                "bytes": path.stat().st_size,
                "sha256": digest_file(path),
            }
        )
    manifest = {"schema_version": 1, "complete": True, "files": manifest_rows}
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--enrollment", required=True, type=Path)
    parser.add_argument("--preregistration", required=True, type=Path)
    parser.add_argument("--read-manifest", required=True, type=Path)
    parser.add_argument("--assembly-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(
        args.source,
        args.enrollment,
        args.preregistration,
        args.read_manifest,
        args.assembly_manifest,
        args.outdir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

