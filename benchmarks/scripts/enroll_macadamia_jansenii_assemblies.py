#!/usr/bin/env python3
"""Verify and enroll the frozen Macadamia CLR and HiFi assemblies."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.qc_reference_fasta import reference_qc


EXPECTED_ARCHIVE_BYTES = 2_564_233_808
EXPECTED_ARCHIVE_MD5 = "f2be4b47799a1b37ca54dda6e70abc63"
EXPECTED_NEW_BYTES = 737_617_956
EXPECTED_NEW_MD5 = "e44711320edf9382a44b227c312f99f6"
OLD_MEMBER = "PacBio/PacBio_falcon_purgeHaplotigs_cns_p_ctg.fasta"
EXPECTED_QC = {
    "old": {"contig_count": 762, "total_bases": 758_277_953, "n50_bp": 1_585_406},
    "new": {"contig_count": 284, "total_bases": 737_613_980, "n50_bp": 4_485_544},
}


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def n50(contigs: list[dict[str, object]]) -> int:
    lengths = sorted((int(row["length_bp"]) for row in contigs), reverse=True)
    half = sum(lengths) / 2
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    raise ValueError("cannot compute N50 for an empty assembly")


def validate_members(archive: tarfile.TarFile) -> None:
    seen: set[str] = set()
    for member in archive.getmembers():
        path = PurePosixPath(member.name)
        if (
            member.name in seen
            or path.is_absolute()
            or ".." in path.parts
            or member.issym()
            or member.islnk()
            or member.isdev()
        ):
            raise ValueError(f"unsafe or duplicate TAR member: {member.name}")
        seen.add(member.name)


def enroll(
    archive_path: Path,
    new_path: Path,
    outdir: Path,
    *,
    expected_archive_bytes: int = EXPECTED_ARCHIVE_BYTES,
    expected_archive_md5: str = EXPECTED_ARCHIVE_MD5,
    expected_new_bytes: int = EXPECTED_NEW_BYTES,
    expected_new_md5: str = EXPECTED_NEW_MD5,
    old_member: str = OLD_MEMBER,
    expected_qc: dict[str, dict[str, int]] = EXPECTED_QC,
) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"output directory already exists: {outdir}")
    if (
        archive_path.stat().st_size != expected_archive_bytes
        or file_md5(archive_path) != expected_archive_md5
    ):
        raise ValueError("old assembly archive failed the frozen size or MD5 gate")
    if new_path.stat().st_size != expected_new_bytes or file_md5(new_path) != expected_new_md5:
        raise ValueError("new assembly failed the frozen size or MD5 gate")

    temporary = outdir.with_name(outdir.name + ".partial")
    if temporary.exists():
        raise ValueError(f"retained partial enrollment exists: {temporary}")
    temporary.mkdir(parents=True)
    old_target = temporary / "old_assembly.fa"
    new_target = temporary / "new_assembly.fa"
    with tarfile.open(archive_path, "r:gz") as archive:
        validate_members(archive)
        try:
            member = archive.getmember(old_member)
        except KeyError as error:
            raise ValueError(f"frozen old assembly member is absent: {old_member}") from error
        source = archive.extractfile(member)
        if source is None:
            raise ValueError(f"frozen old assembly member is not a file: {old_member}")
        with source, old_target.open("wb") as destination:
            shutil.copyfileobj(source, destination, length=1024 * 1024)
        if old_target.stat().st_size != member.size:
            raise OSError("old assembly extraction was incomplete")
    shutil.copyfile(new_path, new_target)

    qc: dict[str, dict[str, object]] = {}
    for role, path in (("old", old_target), ("new", new_target)):
        result = reference_qc(path)
        result["n50_bp"] = n50(result["contigs"])
        observed = {key: result[key] for key in ("contig_count", "total_bases", "n50_bp")}
        if observed != expected_qc[role]:
            raise ValueError(
                f"{role} assembly does not match the frozen published-result gate: "
                f"observed={observed} expected={expected_qc[role]}"
            )
        qc[role] = result

    receipt = {
        "schema_version": 1,
        "complete": True,
        "scope": "paper_same_sample_Macadamia_CLR_Falcon_vs_HiFi_IPA_reference_proxy",
        "source_dois": {"old_data": "10.5524/100812", "new_data": "10.5524/100906"},
        "source_paper_dois": {
            "old": "10.1093/gigascience/giaa146",
            "new": "10.46471/gigabyte.24",
        },
        "inputs": {
            "old_archive": {
                "path": str(archive_path.resolve()),
                "bytes": expected_archive_bytes,
                "md5": expected_archive_md5,
                "sha256": digest_file(archive_path),
                "archive_member": old_member,
            },
            "new_assembly": {
                "path": str(new_path.resolve()),
                "bytes": expected_new_bytes,
                "md5": expected_new_md5,
                "sha256": digest_file(new_path),
            },
        },
        "files": {
            role: {"path": path.name, "bytes": path.stat().st_size, "sha256": digest_file(path)}
            for role, path in (("old", old_target), ("new", new_target))
        },
        "fasta_qc": qc,
        "donor_match_basis": (
            "the update paper calls the HiFi material the same sample previously used for "
            "the CLR comparison"
        ),
        "provenance_warning": (
            "archive BioSample identifiers, aliases, and collection dates differ; same DNA "
            "extraction is not established"
        ),
        "evidence_boundary": (
            "paper-level same-sample high-quality-reference proxy; the new assembly shares "
            "HiFi evidence with TandemX and is not absolute independent copy-number truth"
        ),
        "script_sha256": digest_file(Path(__file__)),
    }
    (temporary / "enrollment_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.rename(outdir)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-archive", required=True, type=Path)
    parser.add_argument("--new-assembly", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    receipt = enroll(args.old_archive, args.new_assembly, args.outdir)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
