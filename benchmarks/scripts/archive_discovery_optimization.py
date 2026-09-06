"""Archive compact, hash-checked full-discovery optimization replays."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, write_table


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def _copy(source: Path, destination: Path, relative: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    expected = digest_file(source)
    if destination.stat().st_size != source.stat().st_size or digest_file(destination) != expected:
        raise OSError(f"Discovery optimization archive copy differs: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source.resolve()),
        "sha256": expected,
        "bytes": destination.stat().st_size,
    }


def validate_pair(label: str, baseline: Path, replay: Path) -> dict:
    baseline = baseline.resolve()
    replay = replay.resolve()
    baseline_environment_path = baseline / "environment.json"
    baseline_execution_path = baseline / "tandemx" / "execution.json"
    baseline_summary_path = baseline / "tandemx" / "discover" / "discovery_summary.json"
    replay_environment_path = replay / "environment.json"
    replay_validation_path = replay / "validation.json"
    replay_summary_path = replay / "discover" / "discovery_summary.json"
    required = (
        baseline_environment_path,
        baseline_execution_path,
        baseline_summary_path,
        replay_environment_path,
        replay_validation_path,
        replay_summary_path,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Missing discovery optimization evidence: {', '.join(missing)}")

    baseline_environment = _json(baseline_environment_path)
    baseline_execution = _json(baseline_execution_path)
    baseline_summary = _json(baseline_summary_path)
    replay_environment = _json(replay_environment_path)
    replay_validation = _json(replay_validation_path)
    replay_summary = _json(replay_summary_path)
    replay_execution = replay_validation.get("execution", {})
    products = replay_validation.get("products", {})
    if (
        replay_validation.get("complete") is not True
        or baseline_execution.get("exit_code") != 0
        or baseline_execution.get("timed_out") is not False
        or replay_execution.get("exit_code") != 0
        or replay_execution.get("timed_out") is not False
        or not isinstance(products, dict)
        or len(products) < 6
        or not all(row.get("byte_identical") is True for row in products.values())
        or Path(replay_environment.get("previous_run", "")).resolve() != baseline
        or replay_environment.get("baseline_environment_sha256")
            != digest_file(baseline_environment_path)
        or replay_environment.get("baseline_execution_sha256")
            != digest_file(baseline_execution_path)
        or replay_environment.get("input_sha256")
            != baseline_environment.get("input", {}).get("fasta_sha256")
    ):
        raise ValueError(f"Incomplete or mismatched discovery replay: {label}")
    for name, receipt in products.items():
        baseline_product = baseline / "tandemx" / "discover" / name
        replay_product = replay / "discover" / name
        if (
            not baseline_product.is_file()
            or not replay_product.is_file()
            or digest_file(baseline_product) != receipt.get("expected_sha256")
            or digest_file(replay_product) != receipt.get("observed_sha256")
        ):
            raise ValueError(f"Discovery replay product hash differs: {label}/{name}")
    summary_fields = (
        "processed_reads", "processed_bases", "candidate_count", "family_count"
    )
    if any(baseline_summary.get(name) != replay_summary.get(name) for name in summary_fields):
        raise ValueError(f"Discovery replay summary differs: {label}")
    baseline_runtime = float(baseline_execution["runtime_seconds"])
    replay_runtime = float(replay_execution["runtime_seconds"])
    baseline_rss = float(baseline_execution["peak_rss_mib"])
    replay_rss = float(replay_execution["peak_rss_mib"])
    if min(baseline_runtime, replay_runtime, baseline_rss, replay_rss) <= 0:
        raise ValueError(f"Discovery replay resources must be positive: {label}")
    return {
        "label": label,
        "baseline": baseline,
        "replay": replay,
        "row": {
            "dataset": label,
            **{name: baseline_summary[name] for name in summary_fields},
            "baseline_runtime_seconds": baseline_runtime,
            "optimized_runtime_seconds": replay_runtime,
            "runtime_change_percent": 100 * (replay_runtime / baseline_runtime - 1),
            "speedup": baseline_runtime / replay_runtime,
            "baseline_peak_rss_mib": baseline_rss,
            "optimized_peak_rss_mib": replay_rss,
            "peak_rss_change_percent": 100 * (replay_rss / baseline_rss - 1),
            "compared_products": len(products),
            "all_products_byte_identical": True,
            "uncompared_new_products": ";".join(
                replay_validation.get("uncompared_new_products", [])
            ),
            "warning": "single_historical_baseline_and_single_replay;engineering_not_publication_timing",
        },
    }


def archive(pairs: list[tuple[str, Path, Path]], outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    checked = [validate_pair(label, baseline, replay) for label, baseline, replay in pairs]
    outdir.mkdir(parents=True)
    rows = [item["row"] for item in checked]
    write_table(outdir / "comparison.tsv", rows, list(rows[0]))
    manifest = []
    for item in checked:
        label = item["label"]
        sources = {
            "baseline_environment.json": item["baseline"] / "environment.json",
            "baseline_execution.json": item["baseline"] / "tandemx" / "execution.json",
            "baseline_discovery_summary.json": (
                item["baseline"] / "tandemx" / "discover" / "discovery_summary.json"
            ),
            "replay_environment.json": item["replay"] / "environment.json",
            "replay_validation.json": item["replay"] / "validation.json",
            "replay_discovery_summary.json": (
                item["replay"] / "discover" / "discovery_summary.json"
            ),
            "replay_stdout.log": item["replay"] / "stdout.log",
            "replay_stderr.log": item["replay"] / "stderr.log",
        }
        for name, source in sources.items():
            relative = Path(label) / name
            manifest.append(_copy(source, outdir / relative, relative))
    comparison = outdir / "comparison.tsv"
    manifest.append({
        "file": "comparison.tsv",
        "source": "derived_from_hash_checked_receipts",
        "sha256": digest_file(comparison),
        "bytes": comparison.stat().st_size,
    })
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pair", action="append", nargs=3, metavar=("LABEL", "BASELINE", "REPLAY"),
        required=True,
    )
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    archive(
        [(label, Path(baseline), Path(replay)) for label, baseline, replay in args.pair],
        args.outdir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
