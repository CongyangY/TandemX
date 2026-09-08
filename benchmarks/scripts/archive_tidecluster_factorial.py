#!/usr/bin/env python3
"""Archive compact success or failure evidence from TideCluster factorial runs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import statistics
from typing import Any

from benchmarks.scripts.run_tidecluster_docker_reference import parse_gnu_time


HEADLINE_METRICS = (
    "array_recall",
    "array_precision",
    "base_union_recall",
    "base_union_precision",
    "matched_boundary_mae_bp",
    "matched_period_mae_bp",
    "cyclic_monomer_recall",
    "homologous_consensus_fraction",
    "tidehunter_wall_seconds",
    "tidehunter_maximum_rss_kb",
    "clustering_wall_seconds",
    "clustering_maximum_rss_kb",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def aggregate_settings(rows: list[dict[str, str]]) -> dict[str, Any]:
    settings = sorted({row["setting"] for row in rows})
    result: dict[str, Any] = {}
    for setting in settings:
        selected = [row for row in rows if row["setting"] == setting]
        metrics: dict[str, Any] = {}
        for name in HEADLINE_METRICS:
            values = [float(row[name]) for row in selected if row.get(name) not in (None, "")]
            metrics[name] = {
                "mean": statistics.fmean(values) if values else None,
                "minimum": min(values) if values else None,
                "maximum": max(values) if values else None,
            }
        result[setting] = {
            "run_count": len(selected),
            "successful_accuracy_run_count": sum(
                row.get("status") == "ok" for row in selected
            ),
            "unavailable_accuracy_run_count": sum(
                row.get("status") != "ok" for row in selected
            ),
            "seeds": sorted(int(row["seed"]) for row in selected),
            "metrics": metrics,
        }
    return result


def candidate_files(source: Path) -> list[Path]:
    paths: set[Path] = set()
    for name in (
        "environment.json",
        "stage_manifest.json",
        "summary.tsv",
        "cell_fates.tsv",
        "stage_fates.tsv",
        "run_receipt.json",
        "independent_verification.json",
    ):
        path = source / name
        if path.is_file():
            paths.add(path)
    for directory in (
        source / "source_snapshot",
        source / "parent_failure_snapshot",
        source / "prior_success_snapshot",
        source / "profile",
    ):
        if directory.is_dir():
            paths.update(path for path in directory.rglob("*") if path.is_file())
    for figure_name in ("figures_v1", "figures_v2"):
        figure_dir = source / figure_name
        if figure_dir.is_dir():
            paths.update(path for path in figure_dir.rglob("*") if path.is_file())
    retained_run_patterns = (
        "*.gnu_time.txt",
        "tc_cmd_args.json",
        "tc_chunks.bed",
        "tc_tidehunter.gff3",
        "tc_clustering.gff3_1.gff3",
        "tc_clustering.gff3",
        "tc_consensus/consensus_sequences_all.fasta",
        "evaluation/*.tsv",
        "evaluation/*.json",
    )
    for seed_dir in source.glob("seed*"):
        if not seed_dir.is_dir():
            continue
        for pattern in retained_run_patterns:
            paths.update(path for path in seed_dir.glob(f"*/{pattern}") if path.is_file())
    return sorted(paths, key=lambda path: path.relative_to(source).as_posix())


def resource_context(source: Path) -> dict[str, Any]:
    stage_path = source / "profile/stages.tsv"
    stage_rows = read_tsv(stage_path) if stage_path.is_file() else []
    internal: dict[str, Any] = {}
    for path in sorted(source.glob("seed*/*/*.gnu_time.txt")):
        relative = path.relative_to(source).as_posix()
        try:
            internal[relative] = parse_gnu_time(path, require_success=False)
        except (OSError, ValueError) as error:
            internal[relative] = {"parse_error": str(error)}
    return {
        "profile_stage_count": len(stage_rows),
        "failed_profile_stages": [
            row for row in stage_rows if int(row.get("exit_code", -1)) != 0
        ],
        "internal_gnu_time": internal,
    }


def archive_description(
    experiment_complete: bool, accuracy_complete: bool = True
) -> str:
    if experiment_complete and accuracy_complete:
        return (
            "This compact archive retains the frozen configuration, native GFF "
            "and consensus outputs, normalization/evaluation tables, external-"
            "stage resources, logs, source snapshot and editable six-panel figure. "
            "The separate verifier recomputed interval, base-union, boundary, "
            "period and cyclic-family endpoints before archival.\n\n"
        )
    if experiment_complete:
        return (
            "This compact archive retains all frozen cell and stage fates, native "
            "outputs for successful cells, failed-stage logs and resources, "
            "normalization/evaluation tables, source snapshots and an editable "
            "six-panel figure. The independent verifier recomputed successful "
            "accuracy cells and confirmed that unavailable cells contain no "
            "invented zero measurements.\n\n"
        )
    return (
        "This compact archive retains the frozen configuration, exact attempted "
        "stage, external-process logs, source snapshot and host/container resource "
        "records. Execution failed before accuracy evaluation, so no figure or "
        "accuracy summary is present and unavailable measurements remain missing "
        "rather than zero.\n\n"
    )


def validate_figure(source: Path) -> dict[str, Any]:
    failed_directory = source / "figures_v1"
    accepted_directory = source / "figures_v2"
    directory = accepted_directory if accepted_directory.is_dir() else failed_directory
    provenance = json.loads(
        (directory / "figure_provenance.json").read_text(encoding="utf-8")
    )
    if provenance.get("complete") is not True or provenance.get("panel_count") != 6:
        raise ValueError("TideCluster validation figure is incomplete")
    if provenance.get("frozen_run_count") != 6:
        raise ValueError("TideCluster validation figure changed the frozen run count")
    if provenance.get("svg_raster_image_element_count") != 0:
        raise ValueError("TideCluster validation SVG contains raster image elements")
    if provenance.get("svg_text_element_count", 0) < 1:
        raise ValueError("TideCluster validation SVG text is not editable")
    for name, record in provenance.get("outputs", {}).items():
        path = directory / name
        if not path.is_file():
            raise FileNotFoundError(f"TideCluster figure output is missing: {path}")
        if path.stat().st_size != record.get("bytes") or digest(path) != record.get(
            "sha256"
        ):
            raise ValueError(f"TideCluster figure output changed: {path}")
    if accepted_directory.is_dir():
        failed_provenance = json.loads(
            (failed_directory / "figure_provenance.json").read_text(encoding="utf-8")
        )
        for name, record in failed_provenance.get("outputs", {}).items():
            path = failed_directory / name
            if not path.is_file() or path.stat().st_size != record.get("bytes"):
                raise ValueError(f"failed TideCluster figure output changed: {path}")
            if digest(path) != record.get("sha256"):
                raise ValueError(f"failed TideCluster figure output changed: {path}")
        if digest(failed_directory / "panel_source.tsv") != digest(
            accepted_directory / "panel_source.tsv"
        ):
            raise ValueError("failed and accepted TideCluster panel sources differ")
        if digest(failed_directory / "tidecluster_factorial_validation.png") == digest(
            accepted_directory / "tidecluster_factorial_validation.png"
        ):
            raise ValueError("failed and accepted TideCluster renders are identical")
        provenance = dict(provenance)
        provenance["qa_history"] = {
            "failed_render": "figures_v1",
            "accepted_render": "figures_v2",
            "failed_reason": "top_and_bottom_legends_overlap_panel_content",
            "panel_source_identical": True,
        }
    return provenance


def archive(source: Path, config_path: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite archive: {outdir}")
    source = source.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    run_receipt_path = source / "run_receipt.json"
    if not run_receipt_path.is_file():
        raise FileNotFoundError(f"TideCluster run receipt is missing: {run_receipt_path}")
    run_receipt = json.loads(run_receipt_path.read_text(encoding="utf-8"))
    experiment_complete = run_receipt.get("complete") is True
    summary_rows: list[dict[str, str]] = []
    independent: dict[str, Any] | None = None
    figure: dict[str, Any] | None = None
    if experiment_complete:
        summary_rows = read_tsv(source / "summary.tsv")
        expected_runs = len(config["datasets"]) * len(config["settings"])
        if len(summary_rows) != expected_runs:
            raise ValueError("TideCluster completed summary changed its frozen run count")
        verification_path = source / "independent_verification.json"
        independent = json.loads(verification_path.read_text(encoding="utf-8"))
        if independent.get("verification_passed") is not True:
            raise ValueError("TideCluster independent accuracy verification failed")
        if independent.get("failures") != []:
            raise ValueError("TideCluster independent verification retained mismatches")
        figure = validate_figure(source)

    outdir.mkdir(parents=True)
    copied: list[dict[str, Any]] = []
    preregistration = outdir / "preregistration.json"
    shutil.copyfile(config_path, preregistration)
    copied.append(
        {
            "file": preregistration.name,
            "source": str(config_path.resolve()),
            "bytes": preregistration.stat().st_size,
            "sha256": digest(preregistration),
        }
    )
    for path in candidate_files(source):
        relative = path.relative_to(source)
        target = outdir / "results" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        if target.stat().st_size != path.stat().st_size or digest(target) != digest(path):
            raise OSError(f"TideCluster archive copy differs: {path}")
        copied.append(
            {
                "file": target.relative_to(outdir).as_posix(),
                "source": str(path),
                "bytes": target.stat().st_size,
                "sha256": digest(target),
            }
        )

    headline = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": experiment_complete,
        "experiment_fate": run_receipt.get(
            "fate",
            "completed_and_independently_verified"
            if experiment_complete
            else "incomplete_unknown_fate",
        ),
        "accuracy_complete": run_receipt.get(
            "accuracy_complete", experiment_complete
        ),
        "run_count": len(summary_rows),
        "setting_summary": aggregate_settings(summary_rows),
        "independent_verification": independent,
        "validation_figure": figure,
        "resource_context": resource_context(source),
        "warning": (
            "same_process_simulated_validation_genomes;technical_repetitions_not_"
            "biological_replicates;failed_processes_are_missing_not_zero_measurements"
        ),
    }
    (outdir / "headline_summary.json").write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    title = run_receipt.get(
        "experiment_id", config.get("experiment_id", "TideCluster factorial validation")
    )
    (outdir / "README.md").write_text(
        f"# {title}\n\n"
        + archive_description(
            experiment_complete,
            run_receipt.get("accuracy_complete", experiment_complete),
        )
        + "The three 10-Mb genomes are same-process simulations and the seeds were "
        "already consumed for TandemX quantification validation. TideCluster "
        "settings were frozen before its outputs were inspected, but the runs "
        "are technical comparisons rather than independent biological replicates.\n",
        encoding="utf-8",
    )
    for name in ("headline_summary.json", "README.md"):
        path = outdir / name
        copied.append(
            {
                "file": name,
                "source": "generated_by_archive_tidecluster_factorial",
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )
    manifest = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": experiment_complete,
        "file_count": len(copied),
        "files": copied,
    }
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
