"""Archive compact evidence for a completed checksum-pinned reference FASTA."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from benchmarks.challenge.schema import digest_file


COMPACT_FILES = ("source_metadata.html", "reference_plan.json", "reference_receipt.json")


def _sha256(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or set(value) - set("0123456789abcdef")
    ):
        raise ValueError(f"Invalid {label} SHA-256")
    return value


def archive(source: Path, outdir: Path) -> dict:
    """Validate a local reference and copy only its compact provenance records."""
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    paths = {name: source / name for name in COMPACT_FILES}
    if any(not path.is_file() for path in paths.values()):
        raise ValueError("Reference metadata, plan and receipt are required")

    plan = json.loads(paths["reference_plan.json"].read_text())
    receipt = json.loads(paths["reference_receipt.json"].read_text())
    transfer = receipt.get("transfer", {})
    qc = receipt.get("qc", {})
    filename = plan.get("filename")
    if not isinstance(filename, str) or not filename or Path(filename).name != filename:
        raise ValueError("Reference plan has an invalid filename")
    reference = source / filename
    if not reference.is_file():
        raise ValueError("The external reference FASTA is missing")

    expected_hash = _sha256(plan.get("expected_sha256"), "plan")
    if (
        receipt.get("complete") is not True
        or qc.get("complete") is not True
        or receipt.get("plan_sha256") != digest_file(paths["reference_plan.json"])
        or plan.get("metadata_sha256") != digest_file(paths["source_metadata.html"])
        or _sha256(transfer.get("sha256"), "transfer") != expected_hash
        or _sha256(qc.get("input_sha256"), "QC") != expected_hash
    ):
        raise ValueError("Reference provenance or completion checks disagree")

    reference_bytes = reference.stat().st_size
    if transfer.get("bytes") != reference_bytes or not 0 < reference_bytes <= plan.get("max_bytes", 0):
        raise ValueError("Reference byte count disagrees with its plan or receipt")
    contigs = qc.get("contigs")
    base_counts = qc.get("base_counts")
    if (
        not isinstance(contigs, list)
        or len(contigs) != qc.get("contig_count")
        or len({row.get("contig") for row in contigs}) != len(contigs)
        or sum(row.get("length_bp", -1) for row in contigs) != qc.get("total_bases")
        or not isinstance(base_counts, dict)
        or sum(base_counts.values()) != qc.get("total_bases")
    ):
        raise ValueError("Reference contig or base-count totals disagree")
    observed_hash = digest_file(reference)
    if observed_hash != expected_hash:
        raise ValueError("External reference FASTA changed after completed QC")

    outdir.mkdir(parents=True)
    archived_files = []
    for name, path in paths.items():
        target = outdir / name
        shutil.copyfile(path, target)
        source_hash = digest_file(path)
        if target.stat().st_size != path.stat().st_size or digest_file(target) != source_hash:
            raise OSError(f"Archive copy differs: {path}")
        archived_files.append({
            "file": name,
            "source": str(path.resolve()),
            "sha256": source_hash,
            "bytes": target.stat().st_size,
        })
    manifest = {
        "complete": True,
        "archived_files": archived_files,
        "external_reference": {
            "filename": filename,
            "source": str(reference.resolve()),
            "sha256": observed_hash,
            "bytes": reference_bytes,
            "copied_to_repository": False,
        },
        "qc": {
            "contig_count": qc["contig_count"],
            "total_bases": qc["total_bases"],
            "n_bases": qc.get("n_bases"),
            "other_ambiguous_bases": qc.get("other_ambiguous_bases"),
        },
        "warning": receipt.get("warning"),
    }
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
