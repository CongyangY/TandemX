from __future__ import annotations

import csv

from tandemx.discover.hierarchy import hierarchy_edge, write_family_hierarchy
from tandemx.discover.mvp import FamilySimilarity


def similarity(
    family_a: str,
    family_b: str,
    length_a: int,
    length_b: int,
    ratio: float,
    relationship: str = "possible_higher_order_or_partial",
) -> FamilySimilarity:
    return FamilySimilarity(
        family_a=family_a,
        family_b=family_b,
        length_a_bp=length_a,
        length_b_bp=length_b,
        kmer_jaccard=0.6,
        shared_kmer_fraction=0.9,
        local_identity=0.95,
        local_overlap_bp=min(length_a, length_b),
        local_overlap_fraction_shorter=1.0,
        length_ratio=ratio,
        orientation="forward",
        relationship=relationship,
        redundant_candidate=relationship == "likely_redundant",
        notes="test evidence",
    )


def test_hierarchy_edge_directs_shorter_to_longer_and_labels_period_multiple() -> None:
    edge = hierarchy_edge(similarity("TXF342", "TXF171", 342, 171, 2.0), 1)

    assert edge is not None
    assert (edge.shorter_family_id, edge.longer_family_id) == ("TXF171", "TXF342")
    assert edge.nearest_integer_multiple == 2
    assert edge.edge_type == "putative_period_multiple"
    assert edge.status == "candidate"
    assert edge.warning == "heuristic_period_multiple_not_validated_hor"


def test_hierarchy_edge_keeps_noninteger_related_pair_unresolved() -> None:
    edge = hierarchy_edge(similarity("TXF171", "TXF250", 171, 250, 250 / 171), 2)

    assert edge is not None
    assert edge.edge_type == "unresolved_related_or_partial"
    assert edge.status == "unresolved"
    assert edge.warning == "related_sequence_without_supported_integer_period_multiple"


def test_hierarchy_excludes_redundant_and_distinct_pairs() -> None:
    assert hierarchy_edge(
        similarity("TXF1", "TXF2", 171, 170, 171 / 170, "likely_redundant"), 1
    ) is None
    assert hierarchy_edge(
        similarity("TXF1", "TXF2", 171, 320, 320 / 171, "distinct"), 1
    ) is None


def test_hierarchy_writes_all_171_342_684_candidate_edges(tmp_path) -> None:
    path = tmp_path / "family_hierarchy.tsv"
    count = write_family_hierarchy(
        path,
        [
            similarity("TXF171", "TXF342", 171, 342, 2.0),
            similarity("TXF171", "TXF684", 171, 684, 4.0),
            similarity("TXF342", "TXF684", 342, 684, 2.0),
        ],
    )

    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert count == 3
    assert [row["hierarchy_edge_id"] for row in rows] == [
        "TXH000000001",
        "TXH000000002",
        "TXH000000003",
    ]
    assert [row["nearest_integer_multiple"] for row in rows] == ["2", "4", "2"]
    assert {row["edge_type"] for row in rows} == {"putative_period_multiple"}
