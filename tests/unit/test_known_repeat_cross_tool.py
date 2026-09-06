import argparse
import csv
from pathlib import Path

import pytest

from benchmarks.scripts.evaluate_known_repeats_across_tools import evaluate, parse_named_path


def test_cross_tool_known_query_scoring_separates_prediction_units(tmp_path: Path) -> None:
    known = tmp_path / "known.fa"
    known.write_text(">Q1\nACGTTGCAACGT\n")
    catalog = tmp_path / "catalog.fa"
    catalog.write_text(">family1\nTGCAACGTACGT\n>other\nAAAAAAAAAAAA\n")
    table = tmp_path / "arrays.tsv"
    table.write_text("read_id\tstart\tend\tsequence\nread1\t5\t17\tTTTTTTTTTTTT\n")
    outdir = tmp_path / "evaluation"

    receipt = evaluate(
        known,
        ["Q1"],
        {"tandemx": catalog},
        {"trf": table},
        outdir,
        "test material",
        "synthetic unit test only",
        1.0,
    )

    rows = list(csv.DictReader((outdir / "known_query_summary.tsv").open(), delimiter="\t"))
    by_tool = {row["tool"]: row for row in rows}
    assert by_tool["tandemx"]["prediction_unit"] == "family_catalog"
    assert by_tool["tandemx"]["recovered_query_count"] == "1"
    assert by_tool["trf"]["prediction_unit"] == "array_consensus"
    assert by_tool["trf"]["recovered_query_count"] == "0"
    details = list(csv.DictReader((outdir / "known_query_details.tsv").open(), delimiter="\t"))
    tandemx = next(row for row in details if row["tool"] == "tandemx")
    assert tandemx["matched_source_occurrence_count"] == "1"
    assert tandemx["matched_source_examples"] == "family1"
    assert tandemx["supported_distinct_consensus_count"] == "1"
    assert tandemx["best_supported_cyclic_edit_similarity"] == "1.0"
    assert receipt["complete"] and set(receipt["outputs"]) == {
        "known_query_summary.tsv", "known_query_details.tsv"
    }


def test_cross_tool_known_query_scoring_rejects_ambiguous_inputs(tmp_path: Path) -> None:
    known = tmp_path / "known.fa"
    known.write_text(">Q1\nACGT\n")
    catalog = tmp_path / "catalog.fa"
    catalog.write_text(">one\nACGT\n")
    with pytest.raises(ValueError, match="unique"):
        evaluate(known, ["Q1", "Q1"], {"x": catalog}, {}, tmp_path / "out", "m", "b")
    with pytest.raises(ValueError, match="unique"):
        evaluate(known, ["Q1"], {"x": catalog}, {"x": catalog}, tmp_path / "out", "m", "b")
    with pytest.raises(argparse.ArgumentTypeError):
        parse_named_path("missing_separator")
