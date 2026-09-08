#!/usr/bin/env python3
"""Verify and enroll the exact Ey15-2 old/new assembly files from Zenodo."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import zipfile

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.qc_reference_fasta import reference_qc


EXPECTED_BUNDLE_BYTES = 1_605_112_570
EXPECTED_BUNDLE_MD5 = "fd3b9bd2f12195aa036a16d2b25d660a"
SOURCE_URL = (
    "https://zenodo.org/records/7326462/files/"
    "SupportingData_A.thaliana_CLR_vs_HiFi.zip?download=1"
)
MEMBERS = {
    "old_assembly.fa": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/Bionano_optical_maps_based_assemblies/"
        "9994.CLR_Canu/CLR_Canu-Arrow.bionano.ragtag.fa"
    ),
    "old_assembly.fa.fai": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/Bionano_optical_maps_based_assemblies/"
        "9994.CLR_Canu/CLR_Canu-Arrow.bionano.ragtag.fa.fai"
    ),
    "new_assembly.fa": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/Bionano_optical_maps_based_assemblies/"
        "9994.HiFi_Hifiasm/HiFi_Hifiasm.bionano.ragtag.fa"
    ),
    "new_assembly.fa.fai": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/Bionano_optical_maps_based_assemblies/"
        "9994.HiFi_Hifiasm/HiFi_Hifiasm.bionano.ragtag.fa.fai"
    ),
    "sensitivity_assembly.fa": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/"
        "Ey15-2_HiFi-Hifiasm_plus_CLR-Canu_assembly/"
        "Hifiasm_Canu-Arrow_2.patch.scaffold.Chr.fa"
    ),
    "sensitivity_repeat_annotation.gff": (
        "SupportingData_A.thaliana_CLR_vs_HiFi/"
        "Ey15-2_HiFi-Hifiasm_plus_CLR-Canu_assembly/"
        "9994.Hifiasm_Canu-Arrow.patch.Repeats_merged.gff"
    ),
}


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_archive_members(archive: zipfile.ZipFile) -> None:
    seen: set[str] = set()
    for info in archive.infolist():
        path = PurePosixPath(info.filename)
        mode = info.external_attr >> 16
        if (
            info.filename in seen
            or path.is_absolute()
            or ".." in path.parts
            or (mode and stat.S_ISLNK(mode))
        ):
            raise ValueError(f"Unsafe or duplicate ZIP member: {info.filename}")
        seen.add(info.filename)


def read_fai(path: Path) -> tuple[int, int]:
    count = total = 0
    identifiers: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 5 or not fields[1].isdigit() or int(fields[1]) <= 0:
                raise ValueError(f"Malformed FAI row {line_number}: {path}")
            if not fields[0] or fields[0] in identifiers:
                raise ValueError(f"Empty or duplicate FAI identifier: {path}")
            identifiers.add(fields[0])
            count += 1
            total += int(fields[1])
    if not count:
        raise ValueError(f"Empty FAI: {path}")
    return count, total


def enroll(
    bundle: Path,
    outdir: Path,
    *,
    expected_bytes: int = EXPECTED_BUNDLE_BYTES,
    expected_md5: str = EXPECTED_BUNDLE_MD5,
    members: dict[str, str] = MEMBERS,
) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    if bundle.stat().st_size != expected_bytes or file_md5(bundle) != expected_md5:
        raise ValueError("Bundle does not match the frozen source size and MD5")
    temporary = outdir.with_name(outdir.name + ".partial")
    if temporary.exists():
        raise ValueError(f"Retained partial enrollment exists: {temporary}")
    temporary.mkdir(parents=True)
    extracted: list[dict[str, object]] = []
    try:
        with zipfile.ZipFile(bundle) as archive:
            validate_archive_members(archive)
            names = set(archive.namelist())
            missing = sorted(set(members.values()) - names)
            if missing:
                raise ValueError(f"Frozen archive members are missing: {missing}")
            for target_name, source_name in members.items():
                target = temporary / target_name
                with archive.open(source_name) as source, target.open("wb") as destination:
                    shutil.copyfileobj(source, destination, length=1024 * 1024)
                info = archive.getinfo(source_name)
                if target.stat().st_size != info.file_size:
                    raise OSError(f"Extracted member size differs: {source_name}")
                extracted.append(
                    {
                        "file": target_name,
                        "archive_member": source_name,
                        "bytes": target.stat().st_size,
                        "sha256": digest_file(target),
                    }
                )
        qc = {}
        for role in ("old", "new", "sensitivity"):
            fasta = temporary / f"{role}_assembly.fa"
            qc[role] = reference_qc(fasta)
            if qc[role]["contig_count"] != 5:
                raise ValueError(f"{role} assembly does not contain five nuclear chromosomes")
        for role in ("old", "new"):
            fai_count, fai_bases = read_fai(temporary / f"{role}_assembly.fa.fai")
            if (fai_count, fai_bases) != (qc[role]["contig_count"], qc[role]["total_bases"]):
                raise ValueError(f"{role} FASTA disagrees with its archived FAI")
        receipt = {
            "schema_version": 1,
            "complete": True,
            "scope": "same_sample_Ey15-2_CLR_Canu_vs_HiFi_Hifiasm_reference_proxy",
            "source_url": SOURCE_URL,
            "source_doi": "10.5281/zenodo.7326462",
            "source_paper_doi": "10.1093/nar/gkac1115",
            "bundle": {
                "path": str(bundle.resolve()),
                "bytes": expected_bytes,
                "md5": expected_md5,
                "sha256": digest_file(bundle),
                "zip_integrity": "passed_by_complete_member_read",
            },
            "files": extracted,
            "fasta_qc": qc,
            "evidence_boundary": (
                "donor_matched_high_quality_reference_proxy;new_assembly_shares_HiFi_evidence;"
                "not_absolute_independent_biological_copy_truth"
            ),
            "script_sha256": digest_file(Path(__file__)),
        }
        (temporary / "enrollment_receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.rename(outdir)
        return receipt
    except Exception:
        # The directory and extracted bytes are evidence of the failed attempt.
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    enroll(args.bundle, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
