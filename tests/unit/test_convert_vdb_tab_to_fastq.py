import gzip
import io
import json
from pathlib import Path

import pytest

from benchmarks.scripts.convert_vdb_tab_to_fastq import convert


def test_convert_streams_biological_rows_to_unique_fastq(tmp_path: Path) -> None:
    source = io.BytesIO(
        b"ACTG\t0,1,2,93\t4\tSRA_READ_TYPE_BIOLOGICAL\n"
        b"NN\t30,31\t2\tSRA_READ_TYPE_BIOLOGICAL\n"
    )
    output = tmp_path / "subset.fastq.gz"
    receipt = tmp_path / "receipt.json"
    result = convert(source, output, receipt, accession="SRR12345678", start_spot=11, end_spot=12)
    assert result["status"] == "complete"
    assert result["record_count"] == 2
    assert result["total_bases"] == 6
    with gzip.open(output, "rt", encoding="ascii") as handle:
        assert handle.read() == "@SRR12345678.11\nACTG\n+\n!\"#~\n@SRR12345678.12\nNN\n+\n?@\n"
    persisted = json.loads(receipt.read_text())
    assert persisted["ids_are_unique"] is True
    assert persisted["selection"] == "fixed_contiguous_spot_range_not_random_sampling"


def test_converter_rejects_nonbiological_row_and_keeps_partial(tmp_path: Path) -> None:
    source = io.BytesIO(
        b"AC\t10,10\t2\tSRA_READ_TYPE_BIOLOGICAL\n"
        b"GT\t10,10\t2\tSRA_READ_TYPE_TECHNICAL\n"
    )
    output = tmp_path / "subset.fastq.gz"
    receipt = tmp_path / "receipt.json"
    with pytest.raises(ValueError, match="disallowed READ_TYPE"):
        convert(source, output, receipt, accession="SRR12345678", start_spot=1, end_spot=2)
    assert not output.exists()
    assert output.with_name(output.name + ".partial").is_file()
    persisted = json.loads(receipt.read_text())
    assert persisted["status"] == "failed_partial_retained"
    assert persisted["completed_records"] == 1
