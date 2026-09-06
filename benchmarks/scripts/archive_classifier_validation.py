"""Archive the frozen classifier development and held-out evidence chain."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, read_table


DEVELOPMENT_V1_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "candidate_metrics.tsv",
    "candidate_summary.tsv",
    "cross_validated_metrics.tsv",
    "selected_metrics.tsv",
    "selection.tsv",
)
DEVELOPMENT_V2_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "candidate_robust_summary.tsv",
    "selected_metrics.tsv",
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
BASELINE_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "localization_metrics.tsv",
    "resource_metrics.tsv",
)
ROBUST_HASHES = {
    "development_validation_sha256": "validation.json",
    "development_selection_sha256": "selection.tsv",
    "development_candidate_summary_sha256": "candidate_robust_summary.tsv",
    "development_selected_metrics_sha256": "selected_metrics.tsv",
    "development_environment_sha256": "environment.json",
    "development_config_sha256": "run_config.json",
}


def confusion(rows: list[dict]) -> dict[str, int | float]:
    counts = Counter(row["outcome"] for row in rows)
    if set(counts) - {"TP", "FN", "FP", "TN", "NA"}:
        raise ValueError("Unknown classifier outcome")
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    return {
        "available": tp + fn + fp + tn,
        "unavailable": counts["NA"],
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else math.nan,
        "false_positive_rate": fp / (fp + tn) if fp + tn else math.nan,
        "precision": tp / (tp + fp) if tp + fp else math.nan,
    }


def validate_paired_classifier_rows(
    rows: list[dict], baseline: str, selected: str, expected_pairs: int
) -> dict[str, list[dict]]:
    methods = {row["method"] for row in rows}
    if methods != {baseline, selected}:
        raise ValueError("Classifier archive requires exactly the frozen method pair")
    key_fields = (
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id",
    )
    grouped = {method: [row for row in rows if row["method"] == method]
               for method in methods}
    keys = {
        method: {tuple(row[name] for name in key_fields) for row in group}
        for method, group in grouped.items()
    }
    if (
        any(len(group) != expected_pairs for group in grouped.values())
        or any(len(keys[method]) != expected_pairs for method in methods)
        or keys[baseline] != keys[selected]
    ):
        raise ValueError("Classifier method rows are incomplete, duplicated or unpaired")
    return grouped


def _json(root: Path, name: str) -> dict:
    return json.loads((root / name).read_text())


def _require_files(root: Path, names: tuple[str, ...]) -> None:
    missing = [name for name in names if not (root / name).is_file()]
    if missing:
        raise ValueError(f"Missing classifier evidence in {root}: {', '.join(missing)}")


def _require_digest(path: Path, expected: object, label: str) -> None:
    if not isinstance(expected, str) or digest_file(path) != expected:
        raise ValueError(f"Classifier evidence hash differs: {label}")


def _validate_snapshot(root: Path, environment: dict) -> list[str]:
    snapshot = Path(environment.get("source_snapshot", ""))
    source_hashes = environment.get("file_hashes")
    benchmark_hashes = environment.get("benchmark_source_sha256")
    if (
        not snapshot.is_dir() or not isinstance(source_hashes, dict)
        or not isinstance(benchmark_hashes, dict) or not benchmark_hashes
    ):
        raise ValueError(f"Missing frozen classifier source snapshot: {root}")
    digest = hashlib.sha256(json.dumps(source_hashes, sort_keys=True).encode()).hexdigest()
    if digest != environment.get("source_digest"):
        raise ValueError("Frozen source digest differs from its manifest")
    for relative, expected in benchmark_hashes.items():
        _require_digest(snapshot / relative, expected, relative)
    return sorted(benchmark_hashes)


def _validate_manifest(root: Path) -> None:
    manifest = _json(root, "archive_manifest.json")
    if not isinstance(manifest, list) or not manifest:
        raise ValueError("Baseline archive manifest is empty")
    for row in manifest:
        relative = row.get("file")
        if not isinstance(relative, str):
            raise ValueError("Baseline archive manifest has an invalid path")
        path = root / relative
        if not path.is_file() or path.stat().st_size != row.get("bytes"):
            raise ValueError(f"Baseline archive file is absent or truncated: {relative}")
        _require_digest(path, row.get("sha256"), relative)


def _copy(origin: Path, target: Path, relative: Path) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(origin, target)
    expected = digest_file(origin)
    if target.stat().st_size != origin.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Classifier archive copy differs: {origin}")
    return {
        "file": relative.as_posix(),
        "source": str(origin.resolve()),
        "sha256": expected,
        "bytes": target.stat().st_size,
    }


def _same_metrics(observed: dict, expected: object, label: str) -> None:
    if not isinstance(expected, dict):
        raise ValueError(f"Missing confusion receipt for {label}")
    for name in ("available", "unavailable", "TP", "FN", "FP", "TN"):
        if observed[name] != expected.get(name):
            raise ValueError(f"Confusion receipt differs for {label}: {name}")
    for name in ("sensitivity", "false_positive_rate", "precision"):
        if not math.isclose(float(observed[name]), float(expected.get(name)), abs_tol=1e-15):
            raise ValueError(f"Metric receipt differs for {label}: {name}")


def validate_chain(
    development_v1: Path,
    development_v2: Path,
    baseline: Path,
    multik: Path,
    heldout: Path,
) -> dict:
    roots = {
        "development_v1": development_v1.resolve(),
        "development_v2": development_v2.resolve(),
        "baseline": baseline.resolve(),
        "multik": multik.resolve(),
        "heldout": heldout.resolve(),
    }
    for name, files in (
        ("development_v1", DEVELOPMENT_V1_FILES),
        ("development_v2", DEVELOPMENT_V2_FILES),
        ("baseline", BASELINE_FILES),
        ("multik", MULTIK_FILES),
        ("heldout", HELDOUT_FILES),
    ):
        _require_files(roots[name], files)
    _validate_manifest(roots["baseline"])

    configs = {name: _json(root, "run_config.json") for name, root in roots.items()}
    envs = {name: _json(root, "environment.json") for name, root in roots.items()}
    vals = {name: _json(root, "validation.json") for name, root in roots.items()}
    source_files = {
        name: _validate_snapshot(root, envs[name])
        for name, root in roots.items() if name != "baseline"
    }
    development = vals["development_v2"].get("development_seeds")
    reserved = vals["development_v2"].get("reserved_heldout_seeds")
    heldout_seeds = envs["heldout"].get("heldout_seeds")
    if (
        vals["development_v1"].get("complete") is not True
        or vals["development_v1"].get("acceptance", {}).get("passed") is not False
        or vals["development_v2"].get("complete") is not True
        or vals["development_v2"].get("acceptance", {}).get("passed") is not True
        or vals["heldout"].get("complete") is not True
        or vals["heldout"].get("acceptance", {}).get("passed") is not False
        or not isinstance(development, list) or not isinstance(reserved, list)
        or reserved != heldout_seeds or set(development) & set(reserved)
        or configs["baseline"] != configs["multik"]
        or configs["baseline"] != configs["heldout"]
    ):
        raise ValueError("Classifier development/held-out split or gate status differs")

    candidate_source_hashes = envs["development_v2"].get("candidate_source_hashes", {})
    v1_map = {
        "validation_sha256": "validation.json",
        "candidate_metrics_sha256": "candidate_metrics.tsv",
        "candidate_summary_sha256": "candidate_summary.tsv",
        "selection_sha256": "selection.tsv",
        "environment_sha256": "environment.json",
        "run_config_sha256": "run_config.json",
    }
    for key, filename in v1_map.items():
        _require_digest(roots["development_v1"] / filename,
                        candidate_source_hashes.get(key), f"development_v1/{filename}")

    model = configs["heldout"].get("classifier_model", {})
    selected = vals["heldout"].get("selected_candidate")
    expected_selected = f"blend_a{float(model.get('blend_alpha')):g}_t{float(model.get('decision_threshold')):g}"
    if (
        model.get("development_seeds") != development
        or selected != vals["development_v2"].get("selected_candidate")
        or selected != expected_selected
        or envs["heldout"].get("classifier_model") != selected
    ):
        raise ValueError("Frozen classifier identity or development seeds differ")
    for key, filename in ROBUST_HASHES.items():
        _require_digest(roots["development_v2"] / filename, model.get(key),
                        f"development_v2/{filename}")

    for filename, key in (
        ("validation.json", "previous_validation_sha256"),
        ("environment.json", "previous_environment_sha256"),
        ("run_config.json", "previous_config_sha256"),
    ):
        _require_digest(roots["baseline"] / filename, envs["multik"].get(key),
                        f"baseline/{filename}")
        _require_digest(roots["multik"] / filename, envs["heldout"].get(key),
                        f"multik/{filename}")

    _require_digest(roots["multik"] / "comparison_metrics.tsv",
                    vals["multik"].get("comparison_metrics_sha256"),
                    "multik/comparison_metrics.tsv")
    _require_digest(roots["heldout"] / "comparison_metrics.tsv",
                    vals["heldout"].get("comparison_metrics_sha256"),
                    "heldout/comparison_metrics.tsv")
    _require_digest(roots["heldout"] / "selected_metrics.tsv",
                    vals["heldout"].get("selected_metrics_sha256"),
                    "heldout/selected_metrics.tsv")

    rows = read_table(roots["heldout"] / "comparison_metrics.tsv", {
        "seed", "unit_substitution_rate", "array_fragments", "coverage",
        "substitution_rate", "assembly_fraction", "family_id", "method", "outcome",
    })
    expected_pairs = int(vals["heldout"].get("paired_family_conditions", 0))
    grouped = validate_paired_classifier_rows(rows, "single_k21", selected, expected_pairs)
    _same_metrics(confusion(grouped["single_k21"]), vals["heldout"].get("baseline_confusion"),
                  "heldout baseline")
    _same_metrics(confusion(grouped[selected]), vals["heldout"].get("selected_confusion"),
                  "heldout selected")
    if {int(row["seed"]) for row in rows} != set(heldout_seeds):
        raise ValueError("Held-out classifier rows contain unexpected seeds")
    return {"roots": roots, "source_files": source_files, "selected": selected}


def archive(
    development_v1: Path,
    development_v2: Path,
    baseline: Path,
    multik: Path,
    heldout: Path,
    outdir: Path,
) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    checked = validate_chain(development_v1, development_v2, baseline, multik, heldout)
    roots: dict[str, Path] = checked["roots"]
    file_groups = {
        "development_v1": DEVELOPMENT_V1_FILES,
        "development_v2": DEVELOPMENT_V2_FILES,
        "baseline": BASELINE_FILES,
        "multik": MULTIK_FILES,
        "heldout": HELDOUT_FILES,
    }
    outdir.mkdir(parents=True)
    manifest = []
    for group, files in file_groups.items():
        for filename in files:
            relative = Path(group) / filename
            manifest.append(_copy(roots[group] / filename, outdir / relative, relative))
    for group, files in checked["source_files"].items():
        environment = _json(roots[group], "environment.json")
        snapshot = Path(environment["source_snapshot"])
        for filename in files:
            relative = Path("source_snapshot") / group / filename
            manifest.append(_copy(snapshot / filename, outdir / relative, relative))
    manifest_path = outdir / "archive_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-v1", required=True, type=Path)
    parser.add_argument("--development-v2", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--multik", required=True, type=Path)
    parser.add_argument("--heldout", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.development_v1, args.development_v2, args.baseline,
            args.multik, args.heldout, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
