"""Archive compact evidence from a completed three-tool real-read comparison."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file, iter_table


ROOT_FILES = ("environment.json", "summary.tsv")
TOOL_FILES = {
    "tandemx": (
        "execution.json", "stderr.log", "stdout.log", "discover/run_config.yaml",
        "discover/run.log", "discover/discovery_summary.json",
        "discover/family_audit_summary.json",
    ),
    "trf": ("execution.json", "stderr.log"),
    "tidehunter": ("execution.json", "stderr.log", "stdout.log"),
}


def archive(source: Path, outdir: Path) -> list[dict]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    environment_path = source / "environment.json"
    summary_path = source / "summary.tsv"
    if not environment_path.is_file() or not summary_path.is_file():
        raise ValueError("Environment and summary are required")
    environment = json.loads(environment_path.read_text())
    rows = list(
        iter_table(
            summary_path,
            {"tool", "exit_code", "timed_out", "normalization", "runtime_seconds", "peak_rss_mib"},
        )
    )
    if (
        {row["tool"] for row in rows} != set(TOOL_FILES)
        or len(rows) != len(TOOL_FILES)
        or any(
            row["exit_code"] != "0"
            or row["timed_out"].lower() != "false"
            or row["normalization"] != "ok"
            for row in rows
        )
    ):
        raise ValueError("Require exactly three successful, normalized tool rows")
    if environment.get("accuracy") != "not_assessed_without_curated_independent_truth":
        raise ValueError("Real-read evidence boundary is missing or unexpected")

    relative_paths = [Path(name) for name in ROOT_FILES]
    for tool, names in TOOL_FILES.items():
        relative_paths.extend(Path(tool) / name for name in names)
    missing = [path for path in relative_paths if not (source / path).is_file()]
    if missing:
        raise ValueError(f"Required compact evidence is missing: {missing}")

    outdir.mkdir(parents=True)
    manifest = []
    for relative in relative_paths:
        source_path = source / relative
        target = outdir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        source_digest = digest_file(source_path)
        if target.stat().st_size != source_path.stat().st_size or digest_file(target) != source_digest:
            raise OSError(f"Archive copy differs: {source_path}")
        manifest.append(
            {
                "file": relative.as_posix(),
                "source": str(source_path.resolve()),
                "sha256": source_digest,
                "bytes": target.stat().st_size,
            }
        )
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
