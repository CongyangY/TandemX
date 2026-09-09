import hashlib
import io
from pathlib import Path
import tarfile

import pytest

from benchmarks.scripts.enroll_macadamia_jansenii_assemblies import enroll


def md5(payload: bytes) -> str:
    return hashlib.md5(payload).hexdigest()


def make_archive(path: Path, member: str, payload: bytes) -> None:
    with tarfile.open(path, "w:gz") as archive:
        info = tarfile.TarInfo(member)
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))


def test_enrolls_exact_members_and_checks_sequence_statistics(tmp_path: Path) -> None:
    old = b">old1\nAAAA\n>old2\nCCCCCC\n"
    new = b">new\nACGTN\n"
    archive = tmp_path / "old.tar.gz"
    make_archive(archive, "old.fa", old)
    new_path = tmp_path / "new.fa"
    new_path.write_bytes(new)

    receipt = enroll(
        archive,
        new_path,
        tmp_path / "enrolled",
        expected_archive_bytes=archive.stat().st_size,
        expected_archive_md5=md5(archive.read_bytes()),
        expected_new_bytes=len(new),
        expected_new_md5=md5(new),
        old_member="old.fa",
        expected_qc={
            "old": {"contig_count": 2, "total_bases": 10, "n50_bp": 6},
            "new": {"contig_count": 1, "total_bases": 5, "n50_bp": 5},
        },
    )

    assert receipt["complete"] is True
    assert receipt["fasta_qc"]["old"]["n50_bp"] == 6
    assert (tmp_path / "enrolled" / "old_assembly.fa").read_bytes() == old


def test_rejects_unsafe_archive_member(tmp_path: Path) -> None:
    old = b">old\nAAAA\n"
    new = b">new\nCCCC\n"
    archive = tmp_path / "old.tar.gz"
    make_archive(archive, "../old.fa", old)
    new_path = tmp_path / "new.fa"
    new_path.write_bytes(new)

    with pytest.raises(ValueError, match="unsafe"):
        enroll(
            archive,
            new_path,
            tmp_path / "enrolled",
            expected_archive_bytes=archive.stat().st_size,
            expected_archive_md5=md5(archive.read_bytes()),
            expected_new_bytes=len(new),
            expected_new_md5=md5(new),
            old_member="../old.fa",
            expected_qc={},
        )
