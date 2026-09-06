"""Archive the frozen depth-gated development and held-out evidence chain."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile

from benchmarks.challenge.schema import digest_file, read_table
from benchmarks.scripts.archive_abundance_evidence import archive as archive_abundance
from benchmarks.scripts.archive_classifier_validation import (
    confusion,
    validate_paired_classifier_rows,
)


BASELINE = "single_k21"
SELECTED = "depth_gated_blend_v3"
DEVELOPMENT_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "comparison_metrics.tsv",
    "comparison_summary.tsv",
    "per_cohort_metrics.tsv",
    "per_seed_metrics.tsv",
    "selection.tsv",
)
MULTIK_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "raw_comparison_metrics.tsv",
    "comparison_metrics.tsv",
    "comparison_summary.tsv",
    "calibration.tsv",
)
HELDOUT_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "comparison_metrics.tsv",
    "comparison_summary.tsv",
    "selected_metrics.tsv",
    "per_seed_metrics.tsv",
)
DEVELOPMENT_HASHES = {
    "development_validation_sha256": "validation.json",
    "development_selection_sha256": "selection.tsv",
    "development_metrics_sha256": "comparison_metrics.tsv",
    "development_environment_sha256": "environment.json",
    "development_config_sha256": "run_config.json",
}


def _json(root: Path, name: str) -> dict:
    return json.loads((root / name).read_text())


def _require_files(root: Path, names: tuple[str, ...]) -> None:
    missing = [name for name in names if not (root / name).is_file()]
    if missing:
        raise ValueError(f"Missing depth-gated evidence in {root}: {', '.join(missing)}")


def _require_digest(path: Path, expected: object, label: str) -> None:
    if not isinstance(expected, str) or digest_file(path) != expected:
        raise ValueError(f"Depth-gated evidence hash differs: {label}")


def _validate_snapshot(root: Path, environment: dict) -> list[str]:
    snapshot = Path(environment.get("source_snapshot", ""))
    source_hashes = environment.get("file_hashes")
    benchmark_hashes = environment.get("benchmark_source_sha256")
    if (
        not snapshot.is_dir()
        or not isinstance(source_hashes, dict)
        or not isinstance(benchmark_hashes, dict)
        or not benchmark_hashes
    ):
        raise ValueError(f"Missing depth-gated source snapshot: {root}")
    source_digest = hashlib.sha256(
        json.dumps(source_hashes, sort_keys=True).encode()
    ).hexdigest()
    if source_digest != environment.get("source_digest"):
        raise ValueError("Depth-gated source digest differs from its manifest")
    for relative, expected in benchmark_hashes.items():
        _require_digest(snapshot / relative, expected, relative)
    return sorted(benchmark_hashes)


def _same_confusion(observed: dict, expected: object, label: str) -> None:
    if not isinstance(expected, dict):
        raise ValueError(f"Missing confusion receipt for {label}")
    for name in ("available", "unavailable", "TP", "FN", "FP", "TN"):
        if observed[name] != expected.get(name):
            raise ValueError(f"Confusion receipt differs for {label}: {name}")
    for name in ("sensitivity", "false_positive_rate", "precision"):
        if not math.isclose(
            float(observed[name]), float(expected.get(name)), rel_tol=0, abs_tol=1e-15
        ):
            raise ValueError(f"Metric receipt differs for {label}: {name}")


def validate_depth_gated_counts(
    selected_rows: list[dict], raw_rows: list[dict], validation: dict, cutoff: float
) -> None:
    """Recompute depth strata and standard-depth blend fallback counts."""
    baseline_rows = [row for row in selected_rows if row["method"] == BASELINE]
    low_depth = sum(float(row["estimated_haploid_depth"]) < cutoff for row in baseline_rows)
    standard_depth = len(baseline_rows) - low_depth
    raw_pairs: dict[tuple[str, ...], dict[str, dict]] = {}
    key_fields = (
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id",
    )
    for row in raw_rows:
        key = tuple(row[name] for name in key_fields)
        raw_pairs.setdefault(key, {})[row["method"]] = row
    fallback = 0
    for pair in raw_pairs.values():
        single = pair.get(BASELINE)
        multi = pair.get("multik_loglinear")
        if single is None or multi is None:
            raise ValueError("Raw depth-gated rows are incomplete")
        if float(single["estimated_haploid_depth"]) < cutoff:
            continue
        unavailable = multi["outcome"] == "NA" or multi["read_estimated_copies"] == "NA"
        nonpositive = (
            not unavailable
            and (
                float(single["read_estimated_copies"]) <= 0
                or float(multi["read_estimated_copies"]) <= 0
            )
        )
        fallback += unavailable or nonpositive
    if (
        low_depth != validation.get("low_depth_rows")
        or standard_depth != validation.get("standard_depth_rows")
        or fallback != validation.get("standard_depth_fallback_rows")
    ):
        raise ValueError("Depth-gated stratum or fallback receipt differs")


def validate_chain(
    development: Path, baseline: Path, multik: Path, heldout: Path
) -> dict:
    roots = {
        "development": development.resolve(),
        "baseline": baseline.resolve(),
        "multik": multik.resolve(),
        "heldout": heldout.resolve(),
    }
    for name, files in (
        ("development", DEVELOPMENT_FILES),
        ("multik", MULTIK_FILES),
        ("heldout", HELDOUT_FILES),
    ):
        _require_files(roots[name], files)
    _require_files(roots["baseline"], ("environment.json", "run_config.json", "validation.json"))

    configs = {name: _json(root, "run_config.json") for name, root in roots.items()}
    envs = {name: _json(root, "environment.json") for name, root in roots.items()}
    vals = {name: _json(root, "validation.json") for name, root in roots.items()}
    source_files = {
        name: _validate_snapshot(roots[name], envs[name])
        for name in ("development", "multik", "heldout")
    }
    model = configs["heldout"].get("classifier_model", {})
    development_seeds = vals["development"].get("consumed_development_seeds")
    heldout_seeds = vals["development"].get("reserved_future_heldout_seeds")
    if (
        vals["development"].get("complete") is not True
        or vals["development"].get("development_only") is not True
        or vals["development"].get("post_failed_heldout_refinement") is not True
        or vals["development"].get("acceptance", {}).get("passed") is not True
        or vals["multik"].get("complete") is not True
        or vals["multik"].get("evaluation_mode")
            != "raw_single_multik_for_predeclared_blend_grid"
        or set(vals["multik"].get("method_confusion", {}))
            != {BASELINE, "multik_loglinear"}
        or vals["heldout"].get("complete") is not True
        or vals["heldout"].get("acceptance", {}).get("passed") is not True
        or vals["heldout"].get("selected_candidate") != SELECTED
        or not isinstance(development_seeds, list)
        or not isinstance(heldout_seeds, list)
        or set(development_seeds) & set(heldout_seeds)
        or model.get("development_seeds") != development_seeds
        or configs["baseline"] != configs["multik"]
        or configs["baseline"] != configs["heldout"]
        or envs["heldout"].get("heldout_seeds") != heldout_seeds
        or envs["heldout"].get("classifier_model") != SELECTED
        or envs["heldout"].get("no_heldout_fit_or_selection") is not True
    ):
        raise ValueError("Depth-gated development/held-out split or gate status differs")

    for key, filename in DEVELOPMENT_HASHES.items():
        _require_digest(roots["development"] / filename, model.get(key), filename)
    for filename, key in (
        ("validation.json", "previous_validation_sha256"),
        ("environment.json", "previous_environment_sha256"),
        ("run_config.json", "previous_config_sha256"),
    ):
        _require_digest(roots["baseline"] / filename, envs["multik"].get(key), filename)
        _require_digest(roots["multik"] / filename, envs["heldout"].get(key), filename)
    _require_digest(
        roots["multik"] / "comparison_metrics.tsv",
        vals["multik"].get("comparison_metrics_sha256"),
        "multik/comparison_metrics.tsv",
    )
    _require_digest(
        roots["heldout"] / "comparison_metrics.tsv",
        vals["heldout"].get("comparison_metrics_sha256"),
        "heldout/comparison_metrics.tsv",
    )
    _require_digest(
        roots["heldout"] / "selected_metrics.tsv",
        vals["heldout"].get("selected_metrics_sha256"),
        "heldout/selected_metrics.tsv",
    )

    selected_rows = read_table(roots["heldout"] / "comparison_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
        "estimated_haploid_depth", "read_estimated_copies",
    })
    expected_pairs = int(vals["heldout"].get("paired_family_conditions", 0))
    grouped = validate_paired_classifier_rows(
        selected_rows, BASELINE, SELECTED, expected_pairs
    )
    _same_confusion(
        confusion(grouped[BASELINE]), vals["heldout"].get("baseline_confusion"),
        "heldout baseline",
    )
    _same_confusion(
        confusion(grouped[SELECTED]), vals["heldout"].get("selected_confusion"),
        "heldout selected",
    )
    raw_rows = read_table(roots["multik"] / "comparison_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
        "estimated_haploid_depth", "read_estimated_copies",
    })
    validate_paired_classifier_rows(
        raw_rows, BASELINE, "multik_loglinear", expected_pairs
    )
    validate_depth_gated_counts(
        selected_rows, raw_rows, vals["heldout"], float(model["low_depth_cutoff"])
    )
    if {int(row["seed"]) for row in selected_rows} != set(heldout_seeds):
        raise ValueError("Held-out depth-gated rows contain unexpected seeds")
    return {"roots": roots, "source_files": source_files}


def _copy(origin: Path, target: Path, relative: Path) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(origin, target)
    expected = digest_file(origin)
    if target.stat().st_size != origin.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Depth-gated archive copy differs: {origin}")
    return {
        "file": relative.as_posix(),
        "source": str(origin.resolve()),
        "sha256": expected,
        "bytes": target.stat().st_size,
    }


def archive(
    development: Path, baseline: Path, multik: Path, heldout: Path, outdir: Path
) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    checked = validate_chain(development, baseline, multik, heldout)
    roots: dict[str, Path] = checked["roots"]
    with tempfile.TemporaryDirectory(prefix="tandemx-depth-gated-") as temporary:
        temporary_baseline = Path(temporary) / "baseline"
        archive_abundance(roots["baseline"], temporary_baseline)
        outdir.mkdir(parents=True)
        shutil.copytree(temporary_baseline, outdir / "baseline")

    manifest = []
    baseline_manifest = _json(outdir / "baseline", "archive_manifest.json")
    for row in baseline_manifest:
        manifest.append({**row, "file": f"baseline/{row['file']}"})
    baseline_manifest_path = outdir / "baseline" / "archive_manifest.json"
    manifest.append({
        "file": "baseline/archive_manifest.json",
        "source": str((roots["baseline"] / "runs").resolve()),
        "sha256": digest_file(baseline_manifest_path),
        "bytes": baseline_manifest_path.stat().st_size,
    })
    groups = {
        "development": DEVELOPMENT_FILES,
        "multik": MULTIK_FILES,
        "heldout": HELDOUT_FILES,
    }
    for group, files in groups.items():
        for filename in files:
            relative = Path(group) / filename
            manifest.append(_copy(roots[group] / filename, outdir / relative, relative))
    for group, files in checked["source_files"].items():
        environment = _json(roots[group], "environment.json")
        snapshot = Path(environment["source_snapshot"])
        for filename in files:
            relative = Path("source_snapshot") / group / filename
            manifest.append(_copy(snapshot / filename, outdir / relative, relative))
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--multik", required=True, type=Path)
    parser.add_argument("--heldout", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.development, args.baseline, args.multik, args.heldout, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
