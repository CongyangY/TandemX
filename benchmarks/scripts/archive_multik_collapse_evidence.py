"""Archive a verified held-out multi-k assembly-comparison evaluation."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, read_table


METHODS = ("single_k21", "multik_loglinear", "multik_fallback_depth_rule")
ROOT_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "comparison_metrics.tsv",
    "comparison_summary.tsv",
    "calibration.tsv",
)
SOURCE_FILES = (
    "benchmarks/abundance/evaluate_multik_collapse.py",
    "tandemx/compare/mvp.py",
    "tandemx/quantify/multik.py",
    "tandemx/quantify/mvp.py",
    "tandemx/utils/kmers.py",
    "rust-core/src/lib.rs",
    "rust-core/src/weighted_words.rs",
)
CALIBRATION_FILES = {
    "validation.json": "calibration_validation_sha256",
    "calibration.tsv": "calibration_table_sha256",
    "environment.json": "calibration_environment_sha256",
}
BASELINE_FILES = {
    "validation.json": "previous_validation_sha256",
    "environment.json": "previous_environment_sha256",
    "run_config.json": "previous_config_sha256",
}


def _confusion(rows: list[dict]) -> dict[str, dict[str, int]]:
    result = {}
    for method in METHODS:
        counts = Counter(row["outcome"] for row in rows if row["method"] == method)
        result[method] = {
            "available": sum(counts[name] for name in ("TP", "FN", "FP", "TN")),
            "unavailable": counts["NA"],
            "TP": counts["TP"],
            "FN": counts["FN"],
            "FP": counts["FP"],
            "TN": counts["TN"],
        }
    return result


def _copy(origin: Path, target: Path, relative: Path) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(origin, target)
    expected = digest_file(origin)
    if target.stat().st_size != origin.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Archive copy differs: {origin}")
    return {
        "file": relative.as_posix(),
        "source": str(origin.resolve()),
        "sha256": expected,
        "bytes": target.stat().st_size,
    }


def _validate_sources(source: Path, environment: dict) -> Path:
    snapshot = Path(environment["source_snapshot"])
    file_hashes = environment.get("file_hashes", {})
    expected_digest = hashlib.sha256(
        json.dumps(file_hashes, sort_keys=True).encode()
    ).hexdigest()
    if expected_digest != environment.get("source_digest"):
        raise ValueError("Source digest does not match its per-file manifest")
    benchmark_hashes = environment.get("benchmark_source_sha256", {})
    for relative, expected in {**file_hashes, **benchmark_hashes}.items():
        path = snapshot / relative
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"Frozen source differs: {relative}")
    if digest_file(source / "run_config.json") != environment.get("previous_config_sha256"):
        raise ValueError("Evaluation config does not match the baseline config hash")
    return snapshot


def _validate_rows(source: Path, config: dict, validation: dict) -> list[dict]:
    rows = read_table(source / "comparison_metrics.tsv", {
        "seed", "coverage", "substitution_rate", "assembly_fraction", "family_id",
        "method", "outcome", "decision_threshold",
    })
    if (
        validation.get("complete") is not True
        or validation.get("split") != "heldout"
        or validation.get("heldout_used") is not True
        or validation.get("leave_one_genome_out_folds") != 0
        or len(rows) != validation.get("comparison_family_rows")
        or len(rows) != validation.get("expected_comparison_family_rows")
    ):
        raise ValueError("Held-out multi-k validation is incomplete")
    heldout = {str(seed) for seed in config.get("seeds", {}).get("heldout", [])}
    development = {str(seed) for seed in config.get("seeds", {}).get("development", [])}
    if not heldout or not development or heldout & development or {row["seed"] for row in rows} != heldout:
        raise ValueError("Development and held-out seed semantics differ")
    by_method = {method: [row for row in rows if row["method"] == method] for method in METHODS}
    if set(row["method"] for row in rows) != set(METHODS) or len({len(group) for group in by_method.values()}) != 1:
        raise ValueError("Method rows are missing or unbalanced")
    keys = {
        method: {(row["seed"], row["coverage"], row["substitution_rate"],
                  row["assembly_fraction"], row["family_id"]) for row in group}
        for method, group in by_method.items()
    }
    if len({frozenset(value) for value in keys.values()}) != 1:
        raise ValueError("Methods do not evaluate identical conditions")
    if _confusion(rows) != validation.get("method_confusion"):
        raise ValueError("Validation confusion counts differ from the metric rows")
    return rows


def _validate_frozen_model(source: Path, config: dict, environment: dict) -> tuple[Path, Path]:
    model = config.get("collapse_model", {})
    if (
        model.get("method") != "multik_loglinear_else_single_k21"
        or model.get("k_values") != [15, 21, 27, 31]
        or model.get("low_depth_cutoff") != 2.0
        or model.get("low_depth_threshold") != 0.5
        or model.get("standard_depth_threshold") != 0.6
        or model.get("unavailable_rule") != "fallback_single_k21"
        or model.get("calibration_seeds") != config.get("seeds", {}).get("development")
    ):
        raise ValueError("Frozen collapse model semantics differ")
    calibration_rows = read_table(source / "calibration.tsv", {
        "fold", "training_seeds", "evaluation_seeds", "selected_low_depth_threshold",
        "baseline_TP", "baseline_FN", "baseline_FP", "baseline_TN",
        "calibrated_TP", "calibrated_FN", "calibrated_FP", "calibrated_TN",
    })
    if len(calibration_rows) != 1 or calibration_rows[0]["fold"] != "predeclared_heldout_validation":
        raise ValueError("Held-out output fitted or omitted its predeclared model receipt")
    if calibration_rows[0]["selected_low_depth_threshold"] != str(model["low_depth_threshold"]):
        raise ValueError("Applied threshold differs from the frozen model")
    calibration = Path(model["calibration_result"])
    baseline = Path(environment["previous_result"])
    for name, key in CALIBRATION_FILES.items():
        if not (calibration / name).is_file() or digest_file(calibration / name) != model.get(key):
            raise ValueError(f"Calibration evidence differs: {name}")
    for name, key in BASELINE_FILES.items():
        if not (baseline / name).is_file() or digest_file(baseline / name) != environment.get(key):
            raise ValueError(f"Baseline receipt differs: {name}")
    return calibration, baseline


def archive(source: Path, outdir: Path) -> list[dict]:
    source = source.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    for name in ROOT_FILES:
        if not (source / name).is_file():
            raise ValueError(f"Missing multi-k evidence file: {name}")
    environment = json.loads((source / "environment.json").read_text())
    config = json.loads((source / "run_config.json").read_text())
    validation = json.loads((source / "validation.json").read_text())
    if environment.get("scope") != "predeclared_heldout_known_catalogue_point_estimate_validation":
        raise ValueError("Require the predeclared held-out evaluation scope")
    snapshot = _validate_sources(source, environment)
    _validate_rows(source, config, validation)
    calibration, baseline = _validate_frozen_model(source, config, environment)

    outdir.mkdir(parents=True)
    manifest = []
    for name in ROOT_FILES:
        relative = Path(name)
        manifest.append(_copy(source / name, outdir / relative, relative))
    for name in SOURCE_FILES:
        expected = (environment.get("benchmark_source_sha256", {}).get(name)
                    or environment.get("file_hashes", {}).get(name))
        if expected is None or digest_file(snapshot / name) != expected:
            raise ValueError(f"Selected source is absent or differs: {name}")
        relative = Path("source_snapshot") / name
        manifest.append(_copy(snapshot / name, outdir / relative, relative))
    for name in CALIBRATION_FILES:
        relative = Path("calibration_source") / name
        manifest.append(_copy(calibration / name, outdir / relative, relative))
    for name in BASELINE_FILES:
        relative = Path("baseline_receipt") / name
        manifest.append(_copy(baseline / name, outdir / relative, relative))
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
