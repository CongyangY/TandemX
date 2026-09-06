"""Archive compact, verified evidence from a conditional abundance experiment."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

from benchmarks.challenge.schema import digest_file, read_table, write_table


ROOT_FILES = (
    "environment.json",
    "run_config.json",
    "validation.json",
    "run.log",
    "copy_number_metrics.tsv",
    "copy_number_summary.tsv",
    "localization_metrics.tsv",
    "comparison_metrics.tsv",
    "comparison_summary.tsv",
)
DETECTOR_FILES = (
    "benchmarks/challenge/run.py",
    "benchmarks/challenge/schema.py",
    "tandemx/cli.py",
    "tandemx/compare/mvp.py",
    "tandemx/locate/mvp.py",
    "tandemx/quantify/mvp.py",
    "tandemx/utils/kmers.py",
    "rust-core/src/lib.rs",
    "rust-core/src/weighted_words.rs",
)


def _positive_unique(values: object, name: str) -> list:
    if not isinstance(values, list) or not values or len(set(values)) != len(values):
        raise ValueError(f"{name} must be a nonempty unique list")
    return values


def _validate_source(environment: dict, source: Path) -> Path:
    snapshot = Path(environment.get("source_snapshot", ""))
    file_hashes = environment.get("file_hashes")
    benchmark_hashes = environment.get("benchmark_source_sha256")
    if not snapshot.is_dir() or not isinstance(file_hashes, dict) or not file_hashes:
        raise ValueError("Missing frozen detector source manifest")
    if not isinstance(benchmark_hashes, dict) or not benchmark_hashes:
        raise ValueError("Missing frozen benchmark source manifest")
    expected_digest = hashlib.sha256(json.dumps(file_hashes, sort_keys=True).encode()).hexdigest()
    if expected_digest != environment.get("source_digest"):
        raise ValueError("Source digest does not match the per-file manifest")
    for relative, expected in {**file_hashes, **benchmark_hashes}.items():
        path = snapshot / relative
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"Frozen source differs: {relative}")
    if not re.fullmatch(r"[0-9a-f]{40}", str(environment.get("git_head", ""))):
        raise ValueError("Invalid Git source revision")
    if digest_file(source / "run_config.json") != environment.get("config_sha256"):
        raise ValueError("Configuration hash differs from environment receipt")
    return snapshot


def _validate_matrix(source: Path, environment: dict, config: dict, validation: dict) -> None:
    split = environment.get("split")
    seeds_by_split = config.get("seeds")
    if not isinstance(seeds_by_split, dict) or split not in seeds_by_split:
        raise ValueError("Environment split is absent from the configuration")
    seeds = _positive_unique(seeds_by_split[split], "split seeds")
    if any(seed in seeds for name, group in seeds_by_split.items() if name != split for seed in group):
        raise ValueError("Seed groups are not disjoint")
    periods = _positive_unique(config.get("periods"), "periods")
    copies = _positive_unique(config.get("copies"), "copies")
    if len(periods) != len(copies):
        raise ValueError("Period and copy grids differ")
    coverages = _positive_unique(config.get("coverages"), "coverages")
    errors = _positive_unique(config.get("substitution_rates"), "substitution rates")
    fractions = _positive_unique(config.get("assembly_fractions"), "assembly fractions")
    family_count = len(periods)
    expected = {
        "copy_number_family_rows": len(seeds) * len(coverages) * len(errors) * family_count,
        "localization_family_rows": len(seeds) * len(fractions) * family_count,
        "comparison_family_rows": len(seeds) * len(coverages) * len(errors) * len(fractions) * family_count,
    }
    expected_executions = (
        len(seeds) * len(fractions)
        + len(seeds) * len(coverages) * len(errors)
        + len(seeds) * len(coverages) * len(errors) * len(fractions)
    )
    if (
        validation.get("complete") is not True
        or validation.get("successful") != expected_executions
        or validation.get("executions") != expected_executions
        or any(validation.get(key) != value for key, value in expected.items())
    ):
        raise ValueError("Execution validation does not match the configured matrix")

    tables = {
        "copy_number_metrics.tsv": (expected["copy_number_family_rows"], {"seed", "family_id"}),
        "localization_metrics.tsv": (expected["localization_family_rows"], {"seed", "family_id"}),
        "comparison_metrics.tsv": (expected["comparison_family_rows"], {"seed", "family_id", "outcome"}),
        "copy_number_summary.tsv": (len(coverages) * len(errors), {"coverage", "substitution_rate"}),
        "comparison_summary.tsv": (
            len(coverages) * len(errors) * len(fractions),
            {"coverage", "substitution_rate", "assembly_fraction", "TP", "FN", "FP", "TN"},
        ),
    }
    for name, (row_count, required) in tables.items():
        rows = read_table(source / name, required)
        if len(rows) != row_count:
            raise ValueError(f"Unexpected row count in {name}")
        if "seed" in required and {int(row["seed"]) for row in rows} != set(seeds):
            raise ValueError(f"Unexpected seed set in {name}")


def _receipt_rows(source: Path, expected: int, config: dict, split: str) -> list[dict]:
    seed_count = len(config["seeds"][split])
    expected_labels = Counter({
        "locate": seed_count * len(config["assembly_fractions"]),
        "quantify": seed_count * len(config["coverages"]) * len(config["substitution_rates"]),
        "compare": seed_count * len(config["coverages"]) * len(config["substitution_rates"])
        * len(config["assembly_fractions"]),
    })
    rows = []
    labels = Counter()
    receipt_paths = sorted((source / "runs").glob("**/receipt.json"))
    if len(receipt_paths) != expected:
        raise ValueError("Receipt count differs from the validated execution count")
    for path in receipt_paths:
        receipt = json.loads(path.read_text())
        label = receipt.get("label")
        command = receipt.get("command")
        numeric = [receipt.get(name) for name in (
            "runtime_seconds", "peak_rss_mib", "cpu_user_seconds", "cpu_system_seconds"
        )]
        if (
            label not in expected_labels
            or receipt.get("exit_code") != 0
            or receipt.get("timed_out") is not False
            or not isinstance(command, list)
            or not command
            or any(not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 for value in numeric)
        ):
            raise ValueError(f"Invalid successful execution receipt: {path}")
        labels[label] += 1
        rows.append({
            "receipt": path.relative_to(source).as_posix(),
            "receipt_sha256": digest_file(path),
            "label": label,
            "runtime_seconds": receipt["runtime_seconds"],
            "peak_rss_mib": receipt["peak_rss_mib"],
            "cpu_user_seconds": receipt["cpu_user_seconds"],
            "cpu_system_seconds": receipt["cpu_system_seconds"],
            "command_json": json.dumps(command, separators=(",", ":")),
        })
    if labels != expected_labels:
        raise ValueError("Execution labels differ from the configured matrix")
    return rows


def _copy(origin: Path, target: Path, relative: Path) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(origin, target)
    expected = digest_file(origin)
    if target.stat().st_size != origin.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Archive copy differs: {origin}")
    return {"file": relative.as_posix(), "source": str(origin.resolve()),
            "sha256": expected, "bytes": target.stat().st_size}


def archive(source: Path, outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    for name in ROOT_FILES:
        if not (source / name).is_file():
            raise ValueError(f"Missing abundance evidence file: {name}")
    environment = json.loads((source / "environment.json").read_text())
    config = json.loads((source / "run_config.json").read_text())
    validation = json.loads((source / "validation.json").read_text())
    snapshot = _validate_source(environment, source)
    _validate_matrix(source, environment, config, validation)
    receipts = _receipt_rows(source, validation["executions"], config, environment["split"])

    outdir.mkdir(parents=True)
    manifest = []
    for name in ROOT_FILES:
        relative = Path(name)
        manifest.append(_copy(source / relative, outdir / relative, relative))
    selected = sorted(set(environment["benchmark_source_sha256"]) | set(DETECTOR_FILES))
    for name in selected:
        if name not in environment["file_hashes"] and name not in environment["benchmark_source_sha256"]:
            raise ValueError(f"Selected source is absent from the frozen manifest: {name}")
        relative = Path("source_snapshot") / name
        manifest.append(_copy(snapshot / name, outdir / relative, relative))

    resource_path = outdir / "resource_metrics.tsv"
    write_table(resource_path, receipts, list(receipts[0]))
    manifest.append({"file": resource_path.name, "source": f"{(source / 'runs').resolve()}/**/receipt.json",
                     "sha256": digest_file(resource_path), "bytes": resource_path.stat().st_size})
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
