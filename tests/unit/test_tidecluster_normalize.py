from __future__ import annotations

from pathlib import Path

import pytest

from benchmarks.tidecluster.normalize import (
    normalize_resolved_tidecluster,
    normalize_tidecluster,
)


def test_tidecluster_gff_is_joined_and_normalized(tmp_path: Path) -> None:
    tidehunter = tmp_path / "tidehunter.gff3"
    clustering = tmp_path / "clustering.gff3"
    tidehunter.write_text(
        "##gff-version 3\nchr1\tTideHunter\ttandem_repeat\t11\t30\t100\t.\t.\t"
        "ID=x;consensus_sequence=ACGTA;consensus_length=5;copy_number=4.0\n"
    )
    clustering.write_text(
        "##gff-version 3\nchr1\tTideCluster\ttandem_repeat\t11\t30\t1\t.\t.\t"
        "Name=TRC_1;repeat_type=TR\n"
    )
    records = normalize_tidecluster(tidehunter, clustering)
    assert (records[0].start, records[0].end, records[0].period) == (10, 30, 5)
    assert records[0].family_id == "TRC_1"


def test_tidecluster_join_rejects_missing_tidehunter_provenance(tmp_path: Path) -> None:
    tidehunter = tmp_path / "tidehunter.gff3"
    clustering = tmp_path / "clustering.gff3"
    tidehunter.write_text("##gff-version 3\n")
    clustering.write_text(
        "chr1\tTideCluster\ttandem_repeat\t11\t30\t1\t.\t.\tName=TRC_1\n"
    )
    with pytest.raises(ValueError, match="lacks TideHunter provenance"):
        normalize_tidecluster(tidehunter, clustering)


def test_resolved_join_handles_merged_tidehunter_intervals(tmp_path: Path) -> None:
    tidehunter = tmp_path / "tidehunter.gff3"
    intermediate = tmp_path / "intermediate.gff3"
    clustering = tmp_path / "clustering.gff3"
    tidehunter.write_text(
        "##gff-version 3\n"
        "chr1\tTideHunter\ttandem_repeat\t11\t30\t1\t.\t.\t"
        "ID=rep1;consensus_sequence=ACGTACGTAC;consensus_length=10;copy_number=2\n"
        "chr1\tTideHunter\ttandem_repeat\t25\t40\t1\t.\t.\t"
        "ID=rep2;consensus_sequence=TTTTT;consensus_length=5;copy_number=3\n"
    )
    intermediate.write_text(
        "##gff-version 3\n"
        "chr1\tTideCluster\ttandem_repeat\t11\t40\t1\t.\t.\tName=rep1\n"
    )
    clustering.write_text(
        "##gff-version 3\n"
        "chr1\tTideCluster\ttandem_repeat\t11\t40\t1\t.\t.\tName=TRC_1\n"
    )
    records = normalize_resolved_tidecluster(tidehunter, intermediate, clustering)
    assert len(records) == 1
    assert records[0].start == 10 and records[0].end == 40
    assert records[0].period == 10
    assert records[0].representative_tidehunter_id == "rep1"
    assert records[0].copy_number is None
    assert records[0].copy_number_source == "unavailable_after_interval_merge_or_resolution"
