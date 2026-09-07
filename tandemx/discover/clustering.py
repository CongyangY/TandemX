"""Abundance-ordered monomer clustering with fixed observed representatives.

These are operational sequence clusters, not ancestral or taxonomic families.
Every assignment has a witnessed circular global edit-similarity lower bound.
Representatives never drift and membership is not propagated transitively.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from array import array
import hashlib
from math import floor, fsum
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
    return ("sequence" if discovery_method in {"elastic", "cascade"} else "legacy") if requested == "auto" else requested


def write_membership(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MEMBERSHIP_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _index_words(sequence: str) -> Counter[str]:
    words = circular_words(sequence, min(9, len(sequence)))
    result: Counter[str] = Counter()
    for word,count in words.items():
        result[min(word, word.translate(str.maketrans("ACGT", "TGCA"))[::-1])] += count
    return result


def _append_index(index: dict[str, array], words: Counter[str], representative_id: int) -> None:
    """Store each (representative, multiplicity) in one exact unsigned 64-bit word."""
    if not 0 <= representative_id <= 0xFFFFFFFF or any(not 0 < n <= 0xFFFFFFFF for n in words.values()):
        raise ValueError('Representative index and word multiplicities must fit unsigned 32-bit fields')
    prefix = representative_id << 32
    for word, count in words.items():
        postings = index.get(word)
        if postings is None:
            postings = index[word] = array('Q')
        postings.append(prefix | count)


def _indexed_candidate_ids(length: int, words: Counter[str], index: dict[str, array],
                           representative_lengths: list[int], minimum_identity: float) -> list[int]:
    """Exact necessary-condition gate for the existing oriented q-gram test.

    Collapsing each word/RC pair can only increase multiset overlap. Therefore
    insufficient canonical overlap implies that both oriented comparisons would
    fail their unchanged q-gram bound. Length decisions use identical rounding.
    """
    if minimum_identity <= .9 or length < 20:
        return list(range(len(representative_lengths)))
    thresholds: dict[int,int | None] = {}
    shared: Counter[int] = Counter()
    for word,count in words.items():
        for packed in index.get(word, ()):
            j, other_count = packed >> 32, packed & 0xFFFFFFFF
            other_length=representative_lengths[j]
            if other_length not in thresholds:
                size=max(length,other_length)
                limit=floor((1-minimum_identity)*size+1e-9)
                thresholds[other_length]=(max(0,size-min(9,length,other_length)*limit)
                                          if abs(length-other_length)<=limit else None)
            if thresholds[other_length] is not None:
                shared[j]+=min(count,other_count)
    return sorted(j for j,count in shared.items() if count>=thresholds[representative_lengths[j]])


def cluster_monomers(candidates: Sequence[CandidateRepeat], minimum_support: int,
                     minimum_identity: float = 0.95, backend: str = "python"
                     ) -> tuple[list[RepeatFamily], list[dict]]:
    from tandemx.discover.mvp import RepeatFamily, is_low_complexity, orient_monomer

    if not 0 < minimum_identity <= 1 or minimum_support < 1:
        raise ValueError("Need cluster identity in (0,1] and positive support")
    if backend not in {"python", "rust"}:
        raise ValueError("Clustering backend must be python or rust")
    native_index = None
    if backend == "rust":
        from tandemx.discover.rust_backend import RustRepresentativeIndex
        native_index = RustRepresentativeIndex()
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
    representative_lengths: list[int] = []
    members: list[list[CandidateRepeat]] = []
    index: dict[str, array] = {}
    assignment = {}
    for sequence in ordered:
        words = None if native_index is not None else _index_words(sequence)
        possible = (native_index.candidates_sequence(sequence, minimum_identity) if native_index is not None
                    else _indexed_candidate_ids(len(sequence), words, index, representative_lengths, minimum_identity))
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
            representative_lengths.append(len(sequence))
            members.append([])
            if native_index is None:
                _append_index(index, words, j)
            elif native_index.append_sequence(sequence) != j:
                raise RuntimeError("Native representative index lost insertion order")
            distance = sequence.count("N")
            similarity = 1 - distance / len(sequence)
        members[j].extend(grouped[sequence])
        assignment[sequence] = (j, distance, similarity, [item[0] for item in compatible if item[0] != j])
    # The index is no longer needed while building family/member output objects.
    del index, native_index
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
