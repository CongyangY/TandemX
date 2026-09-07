"""Completion receipts distinguish valid negative discovery from failed runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

FILES = (
    "candidate_reads.tsv",
    "families.tsv",
    "monomers.fa",
    "family_similarity.tsv",
    "family_hierarchy.tsv",
)
OPTIONAL_FILES = ("collapsed_families.tsv", "collapsed_monomers.fa", "family_collapse.tsv",
                  "candidate_monomers.fa", "monomer_membership.tsv", "family_audit_summary.json")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_discovery_summary(outdir: Path, processed_reads: int, processed_bases: int,
                            candidate_count: int, family_count: int) -> None:
    if processed_reads < 1 or processed_bases < 1 or min(candidate_count, family_count) < 0:
        raise ValueError("Discovery completion requires processed input and nonnegative output counts")
    summary = {"schema_version": 1, "status": "completed" if family_count else "no_families",
               "processed_reads": processed_reads, "processed_bases": processed_bases,
               "candidate_count": candidate_count, "family_count": family_count,
               "warning": "none" if family_count else "no_families_passed_configured_filters",
               "output_sha256": {name: file_sha256(outdir / name) for name in (*FILES, *OPTIONAL_FILES)
                                 if name in FILES or (outdir / name).is_file()}}
    target = outdir / "discovery_summary.json"
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    temporary.replace(target)


def verified_empty_output(path: Path) -> bool:
    """Accept only an explicit completed zero-result receipt and intact outputs."""
    if path.name not in {"candidate_reads.tsv", "families.tsv", "monomers.fa", "collapsed_families.tsv", "collapsed_monomers.fa",
                         "candidate_monomers.fa", "monomer_membership.tsv", "tidehunter_import.tsv"}:
        return False
    if path.name == "tidehunter_import.tsv":
        receipt = path.parent / "import_summary.json"
        try:
            summary = json.loads(receipt.read_text(encoding="utf-8"))
            return (
                summary["schema_version"] == 1
                and summary["complete"] is True
                and summary["candidate_count"] == 0
                and summary["family_count"] == 0
                and file_sha256(path) == summary["output_sha256"][path.name]
            )
        except (OSError, ValueError, KeyError, TypeError):
            return False
    receipt = path.parent / "discovery_summary.json"
    if not receipt.is_file():
        return False
    try:
        summary = json.loads(receipt.read_text(encoding="utf-8"))
        if (summary["schema_version"] != 1 or summary["status"] != "no_families"
                or summary["family_count"] != 0 or summary["processed_reads"] < 1
                or summary["processed_bases"] < 1 or summary["candidate_count"] < 0):
            return False
        if path.name in {"candidate_reads.tsv", "candidate_monomers.fa", "monomer_membership.tsv", "tidehunter_import.tsv"} and summary["candidate_count"] != 0:
            return False
        names = (*FILES, path.name) if path.name in OPTIONAL_FILES else FILES
        return all(file_sha256(path.parent / name) == summary["output_sha256"][name] for name in names)
    except (OSError, ValueError, KeyError, TypeError):
        # A malformed receipt never weakens ordinary schema validation.
        return False


def has_verified_empty_catalog(outdir: Path) -> bool:
    path = outdir / "monomers.fa"
    return path.is_file() and path.stat().st_size == 0 and verified_empty_output(path)
