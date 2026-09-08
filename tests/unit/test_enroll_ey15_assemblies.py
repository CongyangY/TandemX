import hashlib
from pathlib import Path
import zipfile

import pytest

from benchmarks.scripts.enroll_ey15_assemblies import enroll


FASTA = b">Chr1\nACGT\n>Chr2\nACGT\n>Chr3\nACGT\n>Chr4\nACGT\n>Chr5\nACGT\n"
FAI = b"Chr1\t4\t6\t4\t5\nChr2\t4\t17\t4\t5\nChr3\t4\t28\t4\t5\nChr4\t4\t39\t4\t5\nChr5\t4\t50\t4\t5\n"


def make_bundle(path: Path, unsafe: bool = False) -> tuple[dict[str, str], str]:
    members = {
        "old_assembly.fa": "root/old.fa",
        "old_assembly.fa.fai": "root/old.fa.fai",
        "new_assembly.fa": "root/new.fa",
        "new_assembly.fa.fai": "root/new.fa.fai",
        "sensitivity_assembly.fa": "root/final.fa",
        "sensitivity_repeat_annotation.gff": "root/final.gff",
    }
    with zipfile.ZipFile(path, "w") as archive:
        for name in ("root/old.fa", "root/new.fa", "root/final.fa"):
            archive.writestr(name, FASTA)
        for name in ("root/old.fa.fai", "root/new.fa.fai"):
            archive.writestr(name, FAI)
        archive.writestr("root/final.gff", "Chr1\ttest\trepeat\t1\t4\t.\t+\t.\tName=x\n")
        if unsafe:
            archive.writestr("../escape", "bad")
    return members, hashlib.md5(path.read_bytes()).hexdigest()


def test_enroll_extracts_only_frozen_members_and_validates_fai(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle.zip"
    members, md5 = make_bundle(bundle)
    outdir = tmp_path / "enrolled"
    receipt = enroll(bundle, outdir, expected_bytes=bundle.stat().st_size, expected_md5=md5, members=members)
    assert receipt["complete"]
    assert receipt["fasta_qc"]["old"]["total_bases"] == 20
    assert {path.name for path in outdir.iterdir()} == {
        *members.keys(), "enrollment_receipt.json"
    }
    with pytest.raises(ValueError, match="already exists"):
        enroll(bundle, outdir, expected_bytes=bundle.stat().st_size, expected_md5=md5, members=members)


def test_enroll_rejects_wrong_source_hash(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle.zip"
    members, _md5 = make_bundle(bundle)
    with pytest.raises(ValueError, match="frozen source"):
        enroll(bundle, tmp_path / "out", expected_bytes=bundle.stat().st_size, expected_md5="0" * 32, members=members)


def test_enroll_rejects_unsafe_unselected_member(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle.zip"
    members, md5 = make_bundle(bundle, unsafe=True)
    with pytest.raises(ValueError, match="Unsafe"):
        enroll(bundle, tmp_path / "out", expected_bytes=bundle.stat().st_size, expected_md5=md5, members=members)
