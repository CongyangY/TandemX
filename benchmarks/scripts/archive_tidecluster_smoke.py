"""Validate and archive a compact TideCluster assembly-comparator smoke run."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import shutil

from benchmarks.challenge.schema import digest_file


VERSIONS = {
    "tidecluster": "1.21.2",
    "tidehunter": "1.4.3",
    "kitehor": "0.13.2",
    "mmseqs2_commit": "747c64cc8db3b4803a0f1194a3f75b3ba9f81bcb",
    "blastn": "2.16.0+",
    "gnu_time": "1.10",
}
ROOT_FILES = (
    "tc_tidehunter.gff3",
    "tc_clustering.gff3",
    "tc_cmd_args.json",
    "tc_consensus/consensus_sequences_all.fasta",
    "tidehunter.gnu_time.txt",
    "clustering.gnu_time.txt",
)
PROFILE_FILES = (
    "profile/receipt.json",
    "profile/stages.tsv",
    "profile/tidehunter.stdout.log",
    "profile/tidehunter.stderr.log",
    "profile/clustering.stdout.log",
    "profile/clustering.stderr.log",
)
NORMALIZED_FILES = (
    "normalized/normalized_arrays.tsv",
    "normalized/matches.tsv",
    "normalized/family_recovery.tsv",
    "normalized/metrics.json",
    "normalized/receipt.json",
)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _gnu_time(path: Path) -> dict[str, float | int]:
    text = path.read_text(encoding="utf-8")
    fields = {
        "user_seconds": r"User time \(seconds\): ([0-9.]+)",
        "system_seconds": r"System time \(seconds\): ([0-9.]+)",
        "maximum_rss_kb": r"Maximum resident set size \(kbytes\): ([0-9]+)",
        "exit_status": r"Exit status: ([0-9]+)",
    }
    result: dict[str, float | int] = {}
    for name, pattern in fields.items():
        match = re.search(pattern, text)
        if match is None:
            raise ValueError(f"GNU time field {name} is missing: {path}")
        result[name] = int(match.group(1)) if name in {"maximum_rss_kb", "exit_status"} else float(match.group(1))
    elapsed = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([0-9:.]+)", text)
    if elapsed is None:
        raise ValueError(f"GNU time elapsed field is missing: {path}")
    parts = [float(value) for value in elapsed.group(1).split(":")]
    result["wall_seconds"] = sum(value * (60 ** index) for index, value in enumerate(reversed(parts)))
    if result["exit_status"] != 0 or result["maximum_rss_kb"] <= 0:
        raise ValueError(f"GNU time does not describe a successful measured command: {path}")
    return result


def _count_gff(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))


def _validate_sha256_image_id(image_id: str) -> str:
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", image_id):
        raise ValueError("Docker image ID must be a sha256: digest")
    return image_id


def validate(source: Path, container_dir: Path, image_id: str) -> dict[str, object]:
    image_id = _validate_sha256_image_id(image_id)
    required = [source / name for name in (*ROOT_FILES, *PROFILE_FILES, *NORMALIZED_FILES)]
    required.extend((container_dir / "Dockerfile", container_dir / "README.md"))
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Required TideCluster evidence is missing: {missing}")

    profile = json.loads((source / "profile/receipt.json").read_text())
    stages = _rows(source / "profile/stages.tsv")
    if (
        profile.get("complete") is not True
        or profile.get("requested_stage_count") != 2
        or profile.get("completed_stage_count") != 2
        or {row["stage"] for row in stages} != {"tidehunter", "clustering"}
        or any(int(row["exit_code"]) != 0 for row in stages)
    ):
        raise ValueError("The two profiled TideCluster stages are not complete")

    normalized_receipt = json.loads((source / "normalized/receipt.json").read_text())
    metrics = json.loads((source / "normalized/metrics.json").read_text())
    if normalized_receipt.get("complete") is not True:
        raise ValueError("Normalized TideCluster receipt is incomplete")
    for name, expected in normalized_receipt.get("outputs", {}).items():
        path = source / "normalized" / name
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"Normalized TideCluster output changed: {name}")
    for name, metadata in normalized_receipt.get("inputs", {}).items():
        path = Path(name)
        if not path.is_file() or digest_file(path) != metadata.get("sha256") or path.stat().st_size != metadata.get("bytes"):
            raise ValueError(f"TideCluster evaluation input changed: {path}")

    expected_metrics = {
        "truth_array_count": 3,
        "predicted_array_count": 3,
        "matched_array_count": 3,
        "array_recall": 1.0,
        "array_precision": 1.0,
        "matched_period_mae_bp": 0.0,
        "truth_family_count": 3,
        "recovered_family_count": 3,
        "sequence_family_recall": 1.0,
    }
    if any(metrics.get(key) != value for key, value in expected_metrics.items()):
        raise ValueError("TideCluster smoke metrics differ from the declared controlled result")
    if (
        not 0 <= metrics.get("base_union_precision", -1) <= 1
        or metrics.get("base_union_recall") != 1.0
        or metrics.get("warning") != "known_planted_assembly_truth;smoke_test_not_publication_scale"
        or _count_gff(source / "tc_tidehunter.gff3") != 3
        or _count_gff(source / "tc_clustering.gff3") != 3
    ):
        raise ValueError("TideCluster smoke coordinates, union metrics or evidence boundary disagree")

    internal = {
        stage: _gnu_time(source / f"{stage}.gnu_time.txt")
        for stage in ("tidehunter", "clustering")
    }
    return {
        "complete": True,
        "image": {
            "id": image_id,
            "os": "linux",
            "architecture": "amd64",
            "versions": VERSIONS,
        },
        "smoke_test": {
            "scope": "controlled_planted_assembly_truth",
            "publication_scale": False,
            "metrics": metrics,
        },
        "internal_gnu_time": internal,
        "host_profile_scope": profile.get("rss_scope"),
        "warning": "toy_smoke_only;container_internal_rss_and_host_docker_client_rss_are_not_interchangeable",
    }


def archive(source: Path, container_dir: Path, image_id: str, outdir: Path) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    summary = validate(source, container_dir, image_id)
    relative_paths = [Path(name) for name in (*ROOT_FILES, *PROFILE_FILES, *NORMALIZED_FILES)]
    outdir.mkdir(parents=True)
    manifest = []
    for relative in relative_paths:
        source_path = source / relative
        target = outdir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        observed = digest_file(target)
        if observed != digest_file(source_path) or target.stat().st_size != source_path.stat().st_size:
            raise OSError(f"Archive copy differs: {source_path}")
        manifest.append({"file": relative.as_posix(), "source": str(source_path.resolve()), "sha256": observed, "bytes": target.stat().st_size})
    for name in ("Dockerfile", "README.md"):
        source_path = container_dir / name
        target = outdir / "container" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        manifest.append({"file": f"container/{name}", "source": str(source_path.resolve()), "sha256": digest_file(target), "bytes": target.stat().st_size})
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest.append({"file": "summary.json", "source": "generated_by_archive_tidecluster_smoke", "sha256": digest_file(outdir / "summary.json"), "bytes": (outdir / "summary.json").stat().st_size})
    archive_manifest = {"complete": True, "files": manifest}
    (outdir / "archive_manifest.json").write_text(json.dumps(archive_manifest, indent=2) + "\n")
    return archive_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--container-dir", required=True, type=Path)
    parser.add_argument("--image-id", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.container_dir, args.image_id, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
