"""Compare monomer-transition multiplicities from a bounded, anchored array.

This route starts *after* independent monomer decomposition. It never aligns a
read to the assembly path and does not use the label-path DP prototype. START
and END nodes retain flank phase. An edge multiset deliberately discards longer
order, so equal graph spectra with different paths must be reported ambiguous.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Sequence


START = "<LEFT_FLANK>"
END = "<RIGHT_FLANK>"


@dataclass(frozen=True)
class ReadPath:
    molecule_id: str
    labels: tuple[str, ...]
    assignment_confidence: tuple[float, ...]
    error_profile: str
    both_flanks_verified: bool = True
    haplotype: str | None = None


@dataclass(frozen=True)
class GraphAudit:
    status: str
    event_class: str
    support_count: int
    discordant_count: int
    rejected_count: int
    informative_error_profiles: int
    graph_distance: int | None
    implicated_transitions: tuple[tuple[str, str, int, int], ...]
    ambiguity: tuple[str, ...]
    warning: str


def transition_counts(labels: Sequence[str]) -> Counter[tuple[str, str]]:
    """Include unique flank sentinels so copy count and phase affect the graph."""
    if any(not label or label in (START, END) for label in labels):
        raise ValueError("Monomer labels must be nonempty and distinct from flank sentinels")
    return Counter(zip((START, *labels), (*labels, END)))


def _edit_distance(left: str, right: str) -> int:
    previous = list(range(len(right) + 1))
    for index, base in enumerate(left, 1):
        current = [index]
        for column, other in enumerate(right, 1):
            current.append(min(previous[column] + 1, current[column - 1] + 1,
                               previous[column - 1] + (base != other)))
        previous = current
    return previous[-1]


def _close_catalogue_labels(
    assembly: Sequence[str], candidate: Sequence[str],
    monomers: Mapping[str, str] | None,
    similarity_cutoff: float,
) -> bool:
    """Reject near-equivalent catalogue labels that drive the graph difference."""
    if monomers is None:
        return False
    involved = set(assembly) | set(candidate)
    motifs = {label: monomers[label].upper() for label in involved}
    for first in involved:
        for second in involved:
            if first >= second:
                continue
            a, b = motifs[first], motifs[second]
            if _edit_distance(a, b) / max(len(a), len(b)) <= similarity_cutoff:
                return True
    return False


def _result(
    status: str, event_class: str, support: int, discordant: int,
    rejected: int, profiles: int, distance: int | None,
    transitions: tuple[tuple[str, str, int, int], ...] = (),
    ambiguity: tuple[str, ...] = (),
) -> GraphAudit:
    return GraphAudit(
        status, event_class, support, discordant, rejected, profiles, distance,
        transitions, ambiguity,
        "research_only;requires_external_flank_molecule_assignment_error_haplotype_qc;"
        "graph_spectrum_is_not_complete_order_or_assembly_error_truth",
    )


def audit_transition_graph(
    assembly_labels: Sequence[str],
    reads: Sequence[ReadPath],
    *,
    monomer_sequences: Mapping[str, str] | None = None,
    min_molecules: int = 3,
    min_error_profiles: int = 2,
    min_assignment_confidence: float = 0.9,
    similarity_cutoff: float = 0.15,
    synthetic_error_free_mode: bool = False,
) -> GraphAudit:
    """Audit local, full-flank paths with conservative graph-level abstention.

    Confidence is an upstream calibrated per-copy assignment value, not a value
    estimated here. Distinct error profiles are required for a discordance call
    to reduce one-platform correlated-error artifacts. They do not guarantee
    independence. Parameters are development defaults, not calibrated cutoffs.
    """
    if any(not label for label in assembly_labels):
        raise ValueError("Assembly labels must be nonempty strings")
    if min_molecules < 2 or min_error_profiles < 1:
        raise ValueError("Require at least two molecules and one error profile")
    if min_error_profiles < 2 and not synthetic_error_free_mode:
        raise ValueError("One error profile is allowed only for exact synthetic reads")
    if not 0 <= min_assignment_confidence <= 1 or not 0 <= similarity_cutoff <= 1:
        raise ValueError("Confidence and similarity thresholds must be in [0, 1]")
    if monomer_sequences is not None:
        all_labels = set(assembly_labels) | {x for read in reads for x in read.labels}
        if all_labels - monomer_sequences.keys():
            raise ValueError("Monomer catalogue does not cover every observed label")
        if any(not monomer_sequences[x] or set(monomer_sequences[x].upper()) - set("ACGT")
               for x in all_labels):
            raise ValueError("Monomer catalogue requires nonempty ACGT sequences")
    assembly_path = tuple(assembly_labels)
    assembly_graph = transition_counts(assembly_path)
    seen: set[str] = set()
    accepted: list[ReadPath] = []
    rejected = 0
    for read in reads:
        if not read.molecule_id or read.molecule_id in seen:
            raise ValueError("Nonempty unique molecule IDs are required")
        seen.add(read.molecule_id)
        if not read.error_profile:
            raise ValueError("Each read needs an error-profile identifier")
        if len(read.labels) != len(read.assignment_confidence):
            raise ValueError("Each monomer must have an assignment confidence")
        if any(not label for label in read.labels):
            raise ValueError("Read labels must be nonempty strings")
        if any(not 0 <= confidence <= 1 for confidence in read.assignment_confidence):
            raise ValueError("Assignment confidences must be in [0, 1]")
        if synthetic_error_free_mode and read.error_profile != "exact_synthetic":
            raise ValueError("Synthetic mode requires exact_synthetic error profiles")
        if not read.both_flanks_verified or any(
            confidence < min_assignment_confidence for confidence in read.assignment_confidence
        ):
            rejected += 1
            continue
        accepted.append(read)
    if len(accepted) < min_molecules:
        return _result("INSUFFICIENT_READ_SUPPORT", "unresolved", len(accepted), 0,
                       rejected, len({r.error_profile for r in accepted}), None,
                       ambiguity=("too_few_qualified_full_flank_molecules",))
    spectra = [transition_counts(read.labels) for read in accepted]
    graph_discordant_count = sum(spectrum != assembly_graph for spectrum in spectra)
    haplotypes = {r.haplotype for r in accepted if r.haplotype is not None}
    if len(haplotypes) > 1:
        return _result("AMBIGUOUS", "mixed_haplotypes", len(accepted),
                       graph_discordant_count,
                       rejected, len({r.error_profile for r in accepted}), None,
                       ambiguity=("multiple_haplotypes_unphased_for_assembly_comparison",))
    graph_signatures = Counter(
        tuple(sorted(spectrum.items())) for spectrum in spectra
    )
    dominant_signature, dominant_count = graph_signatures.most_common(1)[0]
    minority_count = len(accepted) - dominant_count
    if minority_count >= 2 or dominant_count / len(accepted) <= 2 / 3:
        return _result("AMBIGUOUS", "mixed_transition_graphs", len(accepted),
                       graph_discordant_count, rejected,
                       len({r.error_profile for r in accepted}), None,
                       ambiguity=("multiple_supported_graph_spectra",))
    dominant = dict(dominant_signature)
    dominant_reads = [read for read, spectrum in zip(accepted, spectra)
                      if dict(spectrum) == dominant]
    profiles = len({read.error_profile for read in dominant_reads})
    if profiles < min_error_profiles:
        return _result("AMBIGUOUS", "profile_limited", len(accepted),
                       graph_discordant_count, rejected,
                       profiles, None,
                       ambiguity=("dominant_graph_lacks_independent_error_profiles",))
    edges = tuple(sorted((left, right, assembly_graph[(left, right)],
                          dominant.get((left, right), 0))
                         for left, right in set(assembly_graph) | set(dominant)
                         if assembly_graph[(left, right)] != dominant.get((left, right), 0)))
    distance = sum(abs(before - after) for _, _, before, after in edges)
    path_set = {read.labels for read in dominant_reads}
    if len(path_set) > 1 or (not edges and next(iter(path_set)) != assembly_path):
        return _result("AMBIGUOUS", "order_unidentifiable_from_transitions",
                       len(accepted), graph_discordant_count, rejected,
                       profiles, distance, edges,
                       ("equal_edge_spectra_do_not_identify_full_monomer_order",))
    if not edges:
        return _result("SUPPORTED", "graph_concordant", len(accepted),
                       graph_discordant_count,
                       rejected, profiles, 0)
    if monomer_sequences is None:
        return _result("AMBIGUOUS", "catalogue_unverified", len(accepted),
                       dominant_count, rejected, profiles, distance, edges,
                       ("monomer_sequence_separability_unverified",))
    if _close_catalogue_labels(assembly_path, next(iter(path_set)),
                               monomer_sequences, similarity_cutoff):
        return _result("AMBIGUOUS", "label_identity_unresolved", len(accepted),
                       dominant_count, rejected, profiles, distance, edges,
                       ("near_equivalent_monomer_labels",))
    class_name = ("copy_count_or_graph_multiplicity" if len(next(iter(path_set))) != len(assembly_path)
                  else "transition_rewiring")
    return _result("DISCORDANT", class_name, len(accepted), dominant_count,
                   rejected, profiles, distance, edges)


def to_common_prediction(case_id: str, audit: GraphAudit) -> dict[str, object]:
    """Map a graph result to the controlled-edit prediction schema.

    Binary scores implement the frozen >=0.5 event threshold. They are not
    probabilities, and graph evidence alone cannot predict bp edits or a full
    edited path. Abstentions have null scores and stay in the denominator.
    """
    if not case_id:
        raise ValueError("case_id is required")
    decided = audit.status in ("SUPPORTED", "DISCORDANT")
    return {
        "case_id": case_id,
        "status": "ok" if decided else "abstain",
        "event_score": (1.0 if audit.status == "DISCORDANT" else 0.0)
        if decided else None,
        "event_type": None,
        "predicted_edited_label_path": None,
        "predicted_edited_interval_bp": None,
        "predicted_signed_bp_delta": None,
        "predicted_edited_orientation_path": None,
        "predicted_copy_spans_bp": None,
        "reason": audit.event_class if not decided else None,
        "audit_state": audit.status,
    }
