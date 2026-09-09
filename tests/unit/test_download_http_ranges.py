import hashlib
from pathlib import Path

import pytest

from benchmarks.scripts.download_http_ranges import (
    assemble_and_verify,
    planned_ranges,
)


def test_planned_ranges_cover_suffix_without_gaps() -> None:
    assert planned_ranges(3, 14, 4) == [(3, 6), (7, 10), (11, 13)]
    assert planned_ranges(14, 14, 4) == []


@pytest.mark.parametrize(
    ("start", "total", "chunk"),
    [(-1, 10, 2), (11, 10, 2), (0, 0, 2), (0, 10, 0)],
)
def test_planned_ranges_reject_invalid_values(start: int, total: int, chunk: int) -> None:
    with pytest.raises(ValueError, match="invalid range plan"):
        planned_ranges(start, total, chunk)


def test_assemble_and_verify_preserves_exact_bytes(tmp_path: Path) -> None:
    prefix = tmp_path / "prefix"
    prefix.write_bytes(b"abc")
    first = tmp_path / "part1"
    first.write_bytes(b"def")
    second = tmp_path / "part2"
    second.write_bytes(b"ghi")
    output = tmp_path / "complete.bin"
    payload = b"abcdefghi"

    receipt = assemble_and_verify(
        output,
        prefix,
        [first, second],
        len(payload),
        hashlib.md5(payload).hexdigest(),
    )

    assert output.read_bytes() == payload
    assert receipt["sha256"] == hashlib.sha256(payload).hexdigest()
    assert receipt["publisher_md5_verified"] is True


def test_assemble_and_verify_retains_failed_assembly(tmp_path: Path) -> None:
    part = tmp_path / "part"
    part.write_bytes(b"wrong")
    output = tmp_path / "complete.bin"

    with pytest.raises(ValueError, match="failed size or publisher MD5"):
        assemble_and_verify(output, None, [part], 5, "0" * 32)

    assert not output.exists()
    assert output.with_name(output.name + ".assembling").exists()
