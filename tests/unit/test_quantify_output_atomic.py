"""Failure injection for quantification's public TSV commit boundary."""

from __future__ import annotations

import errno
from pathlib import Path

import pytest

from tandemx.quantify import mvp


def temporary_files(outdir: Path) -> list[Path]:
    return list(outdir.glob(".copy_number.tsv.*.tmp"))


@pytest.mark.parametrize("error_number", [errno.EIO, errno.ENOSPC])
def test_fsync_failure_leaves_no_new_public_result_or_temp(tmp_path, monkeypatch, error_number):
    target = tmp_path / "copy_number.tsv"

    def fail_fsync(_fd):
        raise OSError(error_number, "injected disk failure")

    monkeypatch.setattr(mvp.os, "fsync", fail_fsync)
    with pytest.raises(OSError) as error:
        mvp.write_copy_number(target, [])
    assert error.value.errno == error_number
    assert not target.exists()
    assert temporary_files(tmp_path) == []


def test_replace_failure_preserves_prior_result_and_cleans_temp(tmp_path, monkeypatch):
    target = tmp_path / "copy_number.tsv"
    target.write_bytes(b"prior completed result\n")

    def fail_replace(_source, _target):
        raise OSError(errno.EIO, "injected replace failure")

    monkeypatch.setattr(mvp.os, "replace", fail_replace)
    with pytest.raises(OSError, match="injected replace failure"):
        mvp.write_copy_number(target, [])
    assert target.read_bytes() == b"prior completed result\n"
    assert temporary_files(tmp_path) == []


def test_partial_temporary_write_failure_never_exposes_result(tmp_path, monkeypatch):
    target = tmp_path / "copy_number.tsv"
    create_temporary = mvp.tempfile.NamedTemporaryFile

    class FailingWriter:
        def __init__(self, actual):
            self.actual = actual
            self.name = actual.name

        def __enter__(self):
            self.actual.__enter__()
            return self

        def __exit__(self, *args):
            return self.actual.__exit__(*args)

        def write(self, data):
            self.actual.write(data[:12])
            raise OSError(errno.EIO, "injected partial write")

    monkeypatch.setattr(mvp.tempfile, "NamedTemporaryFile",
                        lambda **kwargs: FailingWriter(create_temporary(**kwargs)))
    with pytest.raises(OSError, match="injected partial write"):
        mvp.write_copy_number(target, [])
    assert not target.exists()
    assert temporary_files(tmp_path) == []


def test_successful_write_commits_header_and_cleans_temp(tmp_path):
    target = tmp_path / "copy_number.tsv"
    mvp.write_copy_number(target, [])
    assert target.read_text().startswith("family_id\tmonomer_length\t")
    assert target.read_text().endswith("confidence\twarning\n")
    assert temporary_files(tmp_path) == []


def test_keyboard_interrupt_cleans_temp(tmp_path, monkeypatch):
    target = tmp_path / "copy_number.tsv"

    def interrupt(_fd):
        raise KeyboardInterrupt

    monkeypatch.setattr(mvp.os, "fsync", interrupt)
    with pytest.raises(KeyboardInterrupt):
        mvp.write_copy_number(target, [])
    assert not target.exists()
    assert temporary_files(tmp_path) == []
