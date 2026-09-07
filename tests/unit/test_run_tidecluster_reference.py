from __future__ import annotations

from pathlib import Path

from benchmarks.scripts.run_tidecluster_reference import (
    normalize_real_outputs,
    summarize_records,
)
from benchmarks.tidecluster.normalize import TideClusterRecord


def test_summarize_records_merges_overlaps_and_validates_periods() -> None:
    records = [
        TideClusterRecord("tile1", 10, 50, "TRC_1", 20, "ACGT" * 5, 2.0),
        TideClusterRecord("tile1", 40, 70, "TRC_1", 20, "ACGT" * 5, 1.5),
        TideClusterRecord("tile2", 0, 25, "TRC_2", 25, "A" * 25, 1.0),
    ]
    summary = summarize_records(records, {"tile1": 100, "tile2": 50})
    assert summary["predicted_array_count"] == 3
    assert summary["predicted_family_count"] == 2
    assert summary["predicted_positive_sequence_count"] == 2
    assert summary["predicted_union_bp"] == 85
    assert summary["predicted_union_base_fraction"] == 85 / 150
    assert summary["period_bp_median"] == 20


def test_normalize_real_tidecluster_outputs(tmp_path: Path) -> None:
    tidehunter = tmp_path / "tc_tidehunter.gff3"
    intermediate = tmp_path / "tc_clustering.gff3_1.gff3"
    clustering = tmp_path / "tc_clustering.gff3"
    tidehunter.write_text(
        "##gff-version 3\n"
        "tile1\tTideHunter\ttandem_repeat\t11\t30\t1\t.\t.\t"
        "ID=rep1;consensus_sequence=ACGTACGTAC;consensus_length=10;copy_number=2\n"
    )
    intermediate.write_text(
        "##gff-version 3\n"
        "tile1\tTideCluster\ttandem_repeat\t11\t30\t1\t.\t.\tName=rep1\n"
    )
    clustering.write_text(
        "##gff-version 3\n"
        "tile1\tTideCluster\ttandem_repeat\t11\t30\t1\t.\t.\tName=TRC_1\n"
    )
    summary = normalize_real_outputs(
        tidehunter, intermediate, clustering, {"tile1": 100}, tmp_path / "normalized"
    )
    assert summary["predicted_array_count"] == 1
    assert summary["predicted_union_bp"] == 20
    text = (tmp_path / "normalized/normalized_arrays.tsv").read_text()
    assert text.startswith("sequence_id\tstart\tend\t")
    assert "descriptive_real_reference_output_without_accuracy_truth" in text
