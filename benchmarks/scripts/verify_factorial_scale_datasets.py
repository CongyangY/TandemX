"""Verify streamed factorial-scale datasets without modifying their contents."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import yaml

from benchmarks.challenge.schema import digest_file


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"Required JSON file is missing: {path}")
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _require_equal(observed: object, expected: object, label: str) -> None:
    if observed != expected:
        raise ValueError(f"{label} differs: observed={observed!r}, expected={expected!r}")


def _verify_payloads(folder: Path, files: dict[str, str]) -> tuple[int, int]:
    if not files:
        raise ValueError(f"Manifest has no payload files: {folder}")
    byte_count = 0
    for relative, expected in files.items():
        path = folder / relative
        if not path.is_file():
            raise ValueError(f"Manifest payload is missing: {path}")
        _require_equal(digest_file(path), expected, f"SHA-256 for {path}")
        byte_count += path.stat().st_size
    return len(files), byte_count


def verify_dataset(
    dataset: Path, expected: dict[str, Any] | None = None
) -> dict[str, object]:
    """Verify one complete dataset and return a compact audit summary."""
    dataset = dataset.resolve()
    receipt_path = dataset / "generation_receipt.json"
    receipt = _load_json(receipt_path)
    if receipt.get("complete") is not True:
        raise ValueError(f"Generation receipt is incomplete: {receipt_path}")
    seed = int(receipt["seed"])
    condition_rows = receipt.get("conditions_completed")
    if not isinstance(condition_rows, list) or not condition_rows:
        raise ValueError(f"No completed conditions recorded: {receipt_path}")
    condition_ids = [str(row["condition_id"]) for row in condition_rows]
    if len(set(condition_ids)) != len(condition_ids):
        raise ValueError(f"Duplicate condition IDs: {receipt_path}")

    _require_equal(
        digest_file(dataset / "run_config.json"),
        receipt["config_sha256"],
        f"configuration SHA-256 for seed {seed}",
    )
    _require_equal(
        digest_file(dataset / "length_histogram.tsv"),
        receipt["histogram_sha256"],
        f"length-histogram SHA-256 for seed {seed}",
    )
    if expected is not None:
        _require_equal(seed in expected["seeds"], True, f"allowed seed {seed}")
        _require_equal(receipt.get("split"), expected["split"], f"split for seed {seed}")
        _require_equal(
            receipt.get("validation_used"), True, f"validation-used flag for seed {seed}"
        )
        _require_equal(
            receipt.get("heldout_used"), False, f"held-out-used flag for seed {seed}"
        )
        generation = expected["dataset_generation"]
        _require_equal(
            receipt["config_sha256"], generation["config_sha256"], "frozen config hash"
        )
        _require_equal(
            receipt["histogram_sha256"],
            generation["length_histogram_sha256"],
            "frozen length-histogram hash",
        )
        _require_equal(
            receipt["source_sha256"], generation["source_sha256"], "frozen source hashes"
        )
        _require_equal(
            len(condition_rows),
            int(expected["condition_count_per_seed"]),
            f"condition count for seed {seed}",
        )

    genome_dir = dataset / "genome"
    genome_manifest_path = genome_dir / "manifest.json"
    genome_manifest = _load_json(genome_manifest_path)
    _require_equal(int(genome_manifest["seed"]), seed, f"genome seed for {dataset}")
    _require_equal(
        len(genome_manifest.get("array_specs", [])),
        int(receipt["family_count"]),
        f"genome family count for seed {seed}",
    )
    _require_equal(
        int(genome_manifest["genome_bp"]),
        int(receipt["genome_bp"]),
        f"genome size for seed {seed}",
    )
    genome_files, genome_bytes = _verify_payloads(genome_dir, genome_manifest["files"])
    genome_manifest_hash = digest_file(genome_manifest_path)

    payload_files = genome_files
    payload_bytes = genome_bytes
    read_count = 0
    observed_bases = 0
    manifest_count = 1
    for condition in condition_rows:
        condition_id = str(condition["condition_id"])
        read_dir = dataset / "reads" / condition_id
        manifest_path = read_dir / "manifest.json"
        _require_equal(
            digest_file(manifest_path),
            condition["manifest_sha256"],
            f"read-manifest SHA-256 for seed {seed} {condition_id}",
        )
        manifest = _load_json(manifest_path)
        if manifest.get("complete") is not True:
            raise ValueError(f"Read manifest is incomplete: {manifest_path}")
        _require_equal(
            manifest["genome_manifest_sha256"],
            genome_manifest_hash,
            f"genome-manifest link for seed {seed} {condition_id}",
        )
        _require_equal(
            manifest["length_histogram_sha256"],
            receipt["histogram_sha256"],
            f"histogram link for seed {seed} {condition_id}",
        )
        _require_equal(
            int(manifest["read_count"]),
            int(condition["read_count"]),
            f"read count for seed {seed} {condition_id}",
        )
        _require_equal(
            int(manifest["total_bases"]),
            int(condition["total_bases"]),
            f"observed bases for seed {seed} {condition_id}",
        )
        count, size = _verify_payloads(read_dir, manifest["files"])
        payload_files += count
        payload_bytes += size
        read_count += int(manifest["read_count"])
        observed_bases += int(manifest["total_bases"])
        manifest_count += 1

    return {
        "dataset": str(dataset),
        "seed": seed,
        "split": receipt["split"],
        "condition_count": len(condition_rows),
        "manifest_count": manifest_count,
        "payload_file_count": payload_files,
        "payload_bytes": payload_bytes,
        "read_count": read_count,
        "observed_bases": observed_bases,
        "generation_receipt_sha256": digest_file(receipt_path),
        "status": "passed",
    }


def run(
    datasets: list[Path], validation_config: Path | None, output: Path | None
) -> dict[str, object]:
    expected = yaml.safe_load(validation_config.read_text()) if validation_config else None
    rows = [verify_dataset(dataset, expected) for dataset in datasets]
    seeds = [int(row["seed"]) for row in rows]
    if len(set(seeds)) != len(seeds):
        raise ValueError("Dataset seeds are duplicated")
    if expected is not None:
        _require_equal(seeds, [int(value) for value in expected["seeds"]], "dataset seed order")
    result = {
        "complete": True,
        "status": "passed",
        "dataset_count": len(rows),
        "manifest_count": sum(int(row["manifest_count"]) for row in rows),
        "payload_file_count": sum(int(row["payload_file_count"]) for row in rows),
        "payload_bytes": sum(int(row["payload_bytes"]) for row in rows),
        "read_count": sum(int(row["read_count"]) for row in rows),
        "observed_bases": sum(int(row["observed_bases"]) for row in rows),
        "datasets": rows,
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("datasets", nargs="+", type=Path)
    parser.add_argument("--validation-config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.datasets, args.validation_config, args.output), indent=2))


if __name__ == "__main__":
    main()
