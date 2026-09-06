"""Abundance-ordered monomer clustering with fixed observed representatives.

These are operational sequence clusters, not ancestral or taxonomic families.
Every assignment has a witnessed circular global edit-similarity lower bound.
Representatives never drift and membership is not propagated transitively.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
from math import fsum
from typing import TYPE_CHECKING, Sequence
import csv
from pathlib import Path

from tandemx.discover.distance import circular_words, cyclic_merge_evidence

if TYPE_CHECKING:
    from tandemx.discover.mvp import CandidateRepeat, RepeatFamily


MEMBERSHIP_FIELDS = ["read_id", "candidate_id", "cluster_id", "family_id", "representative_sha256",
                     "edit_distance_upper_bound", "similarity_lower_bound", "minimum_cluster_identity",
                     "compatible_cluster_count", "alternative_cluster_ids", "status", "warning"]


def resolve_clustering_method(requested: str, discovery_method: str) -> str:
    if requested not in {"auto", "legacy", "sequence"}:
        raise ValueError("--clustering-method must be auto, legacy, or sequence")
    return ("sequence" if discovery_method == "elastic" else "legacy") if requested == "auto" else requested


def write_membership(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MEMBERSHIP_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _index_words(sequence: str) -> set[str]:
    words = circular_words(sequence, min(9, len(sequence)))
    return {min(word, word.translate(str.maketrans("ACGT", "TGCA"))[::-1]) for word in words}


def cluster_monomers(candidates: Sequence[CandidateRepeat], minimum_support: int,
                     minimum_identity: float = 0.95, backend: str = "python"
                     ) -> tuple[list[RepeatFamily], list[dict]]:
    from tandemx.discover.mvp import RepeatFamily, is_low_complexity, orient_monomer

    if not 0 < minimum_identity <= 1 or minimum_support < 1:
        raise ValueError("Need cluster identity in (0,1] and positive support")
    grouped: dict[str, list[CandidateRepeat]] = defaultdict(list)
    unresolved = []
    for candidate in candidates:
        if not candidate.sequence or set(candidate.sequence.upper()) - set("ACGTN"):
            raise ValueError("Candidate monomers must contain nonempty ACGTN sequences")
        if candidate.sequence.upper().count("N") / len(candidate.sequence) > 1 - minimum_identity + 1e-12:
            unresolved.append({"read_id": candidate.read_id, "candidate_id": candidate.candidate_id,
                               "cluster_id": "NA", "family_id": "NA", "representative_sha256": "NA",
                               "edit_distance_upper_bound": len(candidate.sequence), "similarity_lower_bound": 0,
                               "minimum_cluster_identity": minimum_identity, "compatible_cluster_count": 0,
                               "alternative_cluster_ids": "", "status": "unresolved_sequence",
                               "warning": "ambiguous_bases_not_clustered"})
            continue
        grouped[orient_monomer(candidate.sequence)].append(candidate)
    ordered = sorted(grouped, key=lambda seq: (-len({c.read_id for c in grouped[seq]}),
                     -sum(c.repeat_span_bp for c in grouped[seq]),
                     -fsum(c.score for c in grouped[seq]) / len(grouped[seq]), seq))
    representatives: list[str] = []
    members: list[list[CandidateRepeat]] = []
    index: dict[str, set[int]] = defaultdict(set)
    assignment = {}
    for sequence in ordered:
        words = _index_words(sequence)
        if minimum_identity <= 0.9 or len(sequence) < 20:
            possible = range(len(representatives))
        else:
            possible = sorted({j for word in words for j in index.get(word, ())})
        compatible = []
        for j in possible:
            evidence = cyclic_merge_evidence(representatives[j], sequence, minimum_identity, backend)
            if evidence is not None:
                compatible.append((j, *evidence))
        if compatible:
            # Prefer the strongest witnessed lower bound, then the pre-established
            # abundance rank. Ambiguity is retained in the membership table.
            selected = min(compatible, key=lambda item: (-item[2], item[0]))
            j, distance, similarity = selected
        else:
            j = len(representatives)
            representatives.append(sequence)
            members.append([])
            for word in words:
                index[word].add(j)
            distance = sequence.count("N")
            similarity = 1 - distance / len(sequence)
        members[j].extend(grouped[sequence])
        assignment[sequence] = (j, distance, similarity, [item[0] for item in compatible if item[0] != j])
    support = [len({c.read_id for c in group}) for group in members]
    retained = sorted((j for j in range(len(members)) if support[j] >= minimum_support),
                      key=lambda j: (-support[j], -sum(c.repeat_span_bp for c in members[j]), representatives[j]))
    family_ids = {j: f"TXF{i:06d}" for i, j in enumerate(retained, 1)}
    families = []
    for j in retained:
        group, sequence = members[j], representatives[j]
        mean = fsum(c.score for c in group) / len(group)
        low = is_low_complexity(sequence)
        warnings = ["operational_monomer_cluster", "observed_consensus_representative",
                    "uncalibrated_confidence", f"cluster_identity={minimum_identity:g}"]
        if low:
            warnings.append("low_complexity_family")
        if "N" in sequence:
            warnings.append("ambiguous_representative_bases")
        families.append(RepeatFamily(family_ids[j], family_ids[j].replace("TXF", "TXM"), sequence,
                        len(sequence), support[j], sum(c.repeat_span_bp for c in group), mean, low,
                        "high" if support[j] >= max(3, minimum_support) and mean >= 0.9 and not low and "N" not in sequence else "medium",
                        ";".join(warnings)))
    rows = unresolved
    for sequence, group in grouped.items():
        j, distance, similarity, alternatives = assignment[sequence]
        for candidate in group:
            rows.append({"read_id": candidate.read_id, "candidate_id": candidate.candidate_id,
                         "cluster_id": f"TXG{j + 1:06d}", "family_id": family_ids.get(j, "NA"),
                         "representative_sha256": hashlib.sha256(representatives[j].encode()).hexdigest(),
                         "edit_distance_upper_bound": distance, "similarity_lower_bound": similarity,
                         "minimum_cluster_identity": minimum_identity,
                         "compatible_cluster_count": 1 + len(alternatives),
                         "alternative_cluster_ids": ";".join(f"TXG{k + 1:06d}" for k in alternatives),
                         "status": "assigned" if j in family_ids else "below_minimum_support",
                         "warning": "multiple_compatible_clusters" if alternatives else ""})
    return families, sorted(rows, key=lambda row: row["candidate_id"])
