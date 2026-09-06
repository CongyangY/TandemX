import gzip
import io

import pytest

from benchmarks.scripts.fetch_ena_subset import BoundedReader, extract_fastq_prefix


def test_bounded_compressed_prefix_extracts_complete_records() -> None:
    source = io.BytesIO(gzip.compress(b"@read/ccs\nACGT\n+\n!!!!\n@second\nTT\n+\n!!\n"))
    reader = BoundedReader(source, 1000)
    destination = io.BytesIO()
    with gzip.GzipFile(fileobj=reader) as fastq:
        count, bases, headers = extract_fastq_prefix(fastq, destination, 1)
    assert (count, bases) == (1, 4)
    assert headers == ["@read/ccs"]
    assert destination.getvalue() == b">read/ccs\nACGT\n"


@pytest.mark.parametrize("record", [b"@r\nACGT\n+\n!!!\n", b"@r\nX\n+\n!\n", b"@r\nA\n-\n!\n", b""])
def test_malformed_or_truncated_fastq_is_rejected(record: bytes) -> None:
    with pytest.raises(ValueError):
        extract_fastq_prefix(io.BytesIO(record), io.BytesIO(), 1)


def test_duplicate_reads_and_exhausted_budget_are_rejected() -> None:
    record = b"@r\nA\n+\n!\n"
    with pytest.raises(ValueError, match="Duplicate"):
        extract_fastq_prefix(io.BytesIO(record * 2), io.BytesIO(), 2)
    reader = BoundedReader(io.BytesIO(b"12345"), 3)
    assert reader.read(10) == b"123"
    with pytest.raises(ValueError, match="budget"):
        reader.read(1)
