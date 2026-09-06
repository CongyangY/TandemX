"""Archive compact, verified evidence from a clustering interface replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file


LABELS = ("word_bridge", "sequence_native")
SOURCE_FILES = (
    "benchmarks/__init__.py",
    "benchmarks/scripts/replay_clustering_isolated.py",
    "benchmarks/scripts/replay_clustering.py",
    "tandemx/discover/clustering.py",
    "tandemx/discover/rust_backend.py",
    "rust-core/src/representative_index.rs",
)


def archive(source: Path, outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    environment = json.loads((source / "environment.json").read_text())
    validation = json.loads((source / "validation.json").read_text())
    if environment.get("comparison") != "native_sequence_interface_ablation":
        raise ValueError("Require a native sequence-interface ablation")
    if not validation.get("complete") or not validation.get("exact_output_parity"):
        raise ValueError("Clustering replay is incomplete or lacks exact parity")
    measurements = {row.get("label"): row for row in validation.get("measurements", [])}
    if set(measurements) != set(LABELS):
        raise ValueError("Require exactly the word_bridge and sequence_native measurements")

    snapshot = Path(environment["source_snapshot"])
    expected_source_digest = hashlib.sha256(
        json.dumps(environment["file_hashes"], sort_keys=True).encode()
    ).hexdigest()
    if expected_source_digest != environment.get("source_digest"):
        raise ValueError("Source digest does not match the per-file manifest")
    for relative, expected in environment["file_hashes"].items():
        path = snapshot / relative
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"Frozen source differs: {relative}")
    for relative, expected in environment.get("helper_hashes", {}).items():
        path = snapshot / relative
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"Frozen helper differs: {relative}")

    relative_paths = [Path("environment.json"), Path("validation.json")]
    output_hashes = set()
    for label in LABELS:
        row = measurements[label]
        if row.get("exit_code") != 0 or row.get("timed_out") or row.get("index_backend") != label:
            raise ValueError(f"Invalid child execution: {label}")
        payload = source / f"{label}.json"
        receipt = source / f"{label}.receipt.json"
        if not payload.is_file() or not receipt.is_file() or digest_file(payload) != row.get("output_sha256"):
            raise ValueError(f"Payload or hash differs: {label}")
        receipt_data = json.loads(receipt.read_text())
        if receipt_data.get("output_sha256") != row.get("output_sha256"):
            raise ValueError(f"Receipt hash differs: {label}")
        output_hashes.add(row["output_sha256"])
        relative_paths.extend(
            Path(name) for name in (f"{label}.receipt.json", f"{label}.stdout.log", f"{label}.stderr.log")
        )
    if len(output_hashes) != 1:
        raise ValueError("Complete clustering payloads are not byte-identical")

    missing_sources = [relative for relative in SOURCE_FILES if not (snapshot / relative).is_file()]
    if missing_sources:
        raise ValueError(f"Required source files are missing: {missing_sources}")
    outdir.mkdir(parents=True)
    manifest = []
    for relative in relative_paths:
        origin = source / relative
        target = outdir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(origin, target)
        manifest.append(_manifest_row(relative, origin, target))
    for name in SOURCE_FILES:
        relative = Path("source_snapshot") / name
        origin = snapshot / name
        target = outdir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(origin, target)
        manifest.append(_manifest_row(relative, origin, target))
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _manifest_row(relative: Path, origin: Path, target: Path) -> dict:
    expected = digest_file(origin)
    if target.stat().st_size != origin.stat().st_size or digest_file(target) != expected:
        raise OSError(f"Archive copy differs: {origin}")
    return {"file": relative.as_posix(), "source": str(origin.resolve()),
            "sha256": expected, "bytes": target.stat().st_size}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
