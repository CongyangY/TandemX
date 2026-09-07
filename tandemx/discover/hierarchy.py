"""Emit an explicit, uncertainty-labelled family architecture graph.

The graph is intentionally pairwise.  A near-integer length ratio plus local
sequence support is evidence for a possible period multiple, but it is not a
validated higher-order repeat (HOR) structure or an ancestral relationship.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from tandemx.discover.mvp import FamilySimilarity


@dataclass(frozen=True)
class FamilyHierarchyEdge:
    hierarchy_edge_id: str
    shorter_family_id: str
    longer_family_id: str
    shorter_length_bp: int
    longer_length_bp: int
    nearest_integer_multiple: int
    length_ratio: float
    multiple_error: float
    local_identity: float
    local_overlap_fraction_shorter: float
    shared_kmer_fraction: float
    orientation: str
    edge_type: str
    status: str
    warning: str


def hierarchy_header() -> str:
    return (
        "hierarchy_edge_id\tshorter_family_id\tlonger_family_id\t"
        "shorter_length_bp\tlonger_length_bp\tnearest_integer_multiple\t"
        "length_ratio\tmultiple_error\tlocal_identity\t"
        "local_overlap_fraction_shorter\tshared_kmer_fraction\torientation\t"
        "edge_type\tstatus\twarning"
    )


def hierarchy_edge(
    similarity: FamilySimilarity, ordinal: int
) -> FamilyHierarchyEdge | None:
    """Convert a related-family warning into a directional evidence edge."""
    if ordinal < 1:
        raise ValueError("Hierarchy edge ordinals must be positive")
    if similarity.relationship != "possible_higher_order_or_partial":
        return None
    if similarity.length_a_bp <= similarity.length_b_bp:
        shorter_id, longer_id = similarity.family_a, similarity.family_b
        shorter_length, longer_length = similarity.length_a_bp, similarity.length_b_bp
    else:
        shorter_id, longer_id = similarity.family_b, similarity.family_a
        shorter_length, longer_length = similarity.length_b_bp, similarity.length_a_bp
    if shorter_length < 1 or longer_length < shorter_length:
        raise ValueError("Family hierarchy edges require positive ordered lengths")
    nearest = max(1, round(similarity.length_ratio))
    multiple_error = abs(similarity.length_ratio - nearest)
    period_multiple = nearest >= 2 and multiple_error <= 0.05
    return FamilyHierarchyEdge(
        hierarchy_edge_id=f"TXH{ordinal:09d}",
        shorter_family_id=shorter_id,
        longer_family_id=longer_id,
        shorter_length_bp=shorter_length,
        longer_length_bp=longer_length,
        nearest_integer_multiple=nearest,
        length_ratio=similarity.length_ratio,
        multiple_error=multiple_error,
        local_identity=similarity.local_identity,
        local_overlap_fraction_shorter=similarity.local_overlap_fraction_shorter,
        shared_kmer_fraction=similarity.shared_kmer_fraction,
        orientation=similarity.orientation,
        edge_type=(
            "putative_period_multiple"
            if period_multiple
            else "unresolved_related_or_partial"
        ),
        status="candidate" if period_multiple else "unresolved",
        warning=(
            "heuristic_period_multiple_not_validated_hor"
            if period_multiple
            else "related_sequence_without_supported_integer_period_multiple"
        ),
    )


def format_hierarchy_edge(edge: FamilyHierarchyEdge) -> str:
    return "\t".join(
        [
            edge.hierarchy_edge_id,
            edge.shorter_family_id,
            edge.longer_family_id,
            str(edge.shorter_length_bp),
            str(edge.longer_length_bp),
            str(edge.nearest_integer_multiple),
            f"{edge.length_ratio:.4f}",
            f"{edge.multiple_error:.4f}",
            f"{edge.local_identity:.4f}",
            f"{edge.local_overlap_fraction_shorter:.4f}",
            f"{edge.shared_kmer_fraction:.4f}",
            edge.orientation,
            edge.edge_type,
            edge.status,
            edge.warning,
        ]
    )


def write_family_hierarchy(
    path: Path, similarities: list[FamilySimilarity]
) -> int:
    """Write the candidate architecture graph and return its edge count."""
    edges = [
        edge
        for ordinal, similarity in enumerate(
            (
                item
                for item in similarities
                if item.relationship == "possible_higher_order_or_partial"
            ),
            1,
        )
        if (edge := hierarchy_edge(similarity, ordinal)) is not None
    ]
    path.write_text(
        "\n".join([hierarchy_header(), *(format_hierarchy_edge(edge) for edge in edges)])
        + "\n",
        encoding="utf-8",
    )
    return len(edges)
