"""Post hoc annotation of discovered monomers against known repeat libraries."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from tandemx.io.sequences import read_sequence_records
from tandemx.simulate.toy import reverse_complement


ANNOTATION_FIELDS = (
    "family_id",
    "monomer_length",
    "best_known_id",
    "best_known_length",
    "best_orientation",
    "shared_kmer_fraction",
    "jaccard",
    "dice",
    "containment_discovered_in_known",
    "containment_known_in_discovered",
    "local_identity",
    "local_overlap_bp",
    "annotation_status",
    "notes",
)


@dataclass(frozen=True)
class RepeatRecord:
    identifier: str
    sequence: str


@dataclass(frozen=True)
class RepeatAnnotation:
    family_id: str
    monomer_length: int
    best_known_id: str
    best_known_length: int
    best_orientation: str
    shared_kmer_fraction: float
    jaccard: float
    dice: float
    containment_discovered_in_known: float
    containment_known_in_discovered: float
    local_identity: float
    local_overlap_bp: int
    annotation_status: str
    notes: str


def parse_family_id(description: str, fallback: str) -> str:
    for part in description.split(";"):
        if part.startswith("family_id="):
            return part.split("=", 1)[1]
    return fallback


def read_discovered_catalog(path: Path) -> list[RepeatRecord]:
    return [
        RepeatRecord(identifier=parse_family_id(record.description, record.id), sequence=record.sequence)
        for record in read_sequence_records(path)
    ]


def read_known_repeats(path: Path) -> list[RepeatRecord]:
    return [
        RepeatRecord(identifier=record.id, sequence=record.sequence)
        for record in read_sequence_records(path)
    ]


def kmers(sequence: str, k: int) -> set[str]:
    if k <= 0:
        raise ValueError("--kmer-size must be positive")
    if len(sequence) < k:
        return set()
    return {
        sequence[index : index + k]
        for index in range(len(sequence) - k + 1)
        if "N" not in sequence[index : index + k]
    }


def kmer_metrics(discovered_sequence: str, known_sequence: str, k: int) -> tuple[float, float, float, float, float]:
    discovered = kmers(discovered_sequence, k)
    known = kmers(known_sequence, k)
    if not discovered or not known:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    shared = len(discovered & known)
    union = len(discovered | known)
    shared_fraction = shared / min(len(discovered), len(known))
    jaccard = shared / union if union else 0.0
    dice = 2.0 * shared / (len(discovered) + len(known))
    containment_discovered = shared / len(discovered)
    containment_known = shared / len(known)
    return shared_fraction, jaccard, dice, containment_discovered, containment_known


def best_local_identity(sequence_a: str, sequence_b: str) -> tuple[float, int]:
    min_overlap = min(50, len(sequence_a), len(sequence_b))
    best_identity = 0.0
    best_overlap = 0
    if min_overlap == 0:
        return best_identity, best_overlap
    for offset in range(-len(sequence_b) + 1, len(sequence_a)):
        start_a = max(0, offset)
        start_b = max(0, -offset)
        overlap = min(len(sequence_a) - start_a, len(sequence_b) - start_b)
        if overlap < min_overlap:
            continue
        matches = sum(
            1
            for index in range(overlap)
            if sequence_a[start_a + index] == sequence_b[start_b + index]
        )
        identity = matches / overlap
        if identity > best_identity or (identity == best_identity and overlap > best_overlap):
            best_identity = identity
            best_overlap = overlap
    return best_identity, best_overlap


def classify_annotation(
    *,
    monomer_length: int,
    known_length: int,
    jaccard: float,
    dice: float,
    containment_discovered_in_known: float,
    containment_known_in_discovered: float,
    local_identity: float,
    local_overlap_bp: int,
) -> tuple[str, str]:
    if dice == 0 and local_identity < 0.5:
        return "no_known_match", "No meaningful k-mer or local-identity match to the known-repeat library."

    shorter = max(1, min(monomer_length, known_length))
    longer = max(monomer_length, known_length)
    length_ratio = longer / shorter
    near_integer_ratio = abs(length_ratio - round(length_ratio)) <= 0.10 and round(length_ratio) >= 2
    overlap_fraction = local_overlap_bp / shorter

    if dice >= 0.8 and local_identity >= 0.8:
        return "strong_known_match", "Strong post hoc match; known repeat was not used during discovery."
    if (
        near_integer_ratio
        and max(containment_discovered_in_known, containment_known_in_discovered) >= 0.65
        and local_identity >= 0.70
    ):
        return (
            "possible_higher_order_match",
            "One sequence may represent a higher-order or dimer-like unit; do not collapse automatically.",
        )
    if max(containment_discovered_in_known, containment_known_in_discovered) >= 0.50 and overlap_fraction >= 0.50:
        return "possible_partial_match", "Partial or nested known-repeat similarity; inspect manually."
    if dice >= 0.20 or local_identity >= 0.65:
        return "weak_known_match", "Weak post hoc match; insufficient alone for biological assignment."
    return "no_known_match", "No confident known-repeat annotation."


def annotate_discovered_repeats(
    catalog: Sequence[RepeatRecord],
    known_repeats: Sequence[RepeatRecord],
    *,
    kmer_size: int = 11,
) -> list[RepeatAnnotation]:
    if kmer_size <= 0:
        raise ValueError("--kmer-size must be positive")
    if not catalog:
        raise ValueError("Discovered catalog contains no monomers")
    if not known_repeats:
        raise ValueError("Known repeat library contains no records")

    annotations: list[RepeatAnnotation] = []
    for monomer in catalog:
        best: RepeatAnnotation | None = None
        for known in known_repeats:
            effective_k = min(kmer_size, len(monomer.sequence), len(known.sequence))
            for orientation, known_sequence in (
                ("forward", known.sequence),
                ("reverse", reverse_complement(known.sequence)),
            ):
                shared, jaccard, dice, containment_discovered, containment_known = kmer_metrics(
                    monomer.sequence,
                    known_sequence,
                    effective_k,
                )
                local_identity, local_overlap = best_local_identity(monomer.sequence, known_sequence)
                status, notes = classify_annotation(
                    monomer_length=len(monomer.sequence),
                    known_length=len(known.sequence),
                    jaccard=jaccard,
                    dice=dice,
                    containment_discovered_in_known=containment_discovered,
                    containment_known_in_discovered=containment_known,
                    local_identity=local_identity,
                    local_overlap_bp=local_overlap,
                )
                candidate = RepeatAnnotation(
                    family_id=monomer.identifier,
                    monomer_length=len(monomer.sequence),
                    best_known_id=known.identifier,
                    best_known_length=len(known.sequence),
                    best_orientation=orientation,
                    shared_kmer_fraction=shared,
                    jaccard=jaccard,
                    dice=dice,
                    containment_discovered_in_known=containment_discovered,
                    containment_known_in_discovered=containment_known,
                    local_identity=local_identity,
                    local_overlap_bp=local_overlap,
                    annotation_status=status,
                    notes=notes,
                )
                if best is None or annotation_sort_key(candidate) > annotation_sort_key(best):
                    best = candidate
        assert best is not None
        annotations.append(best)
    return annotations


def annotation_sort_key(annotation: RepeatAnnotation) -> tuple[float, float, float, int]:
    return (
        annotation.dice,
        max(annotation.containment_discovered_in_known, annotation.containment_known_in_discovered),
        annotation.local_identity,
        annotation.local_overlap_bp,
    )


def write_repeat_annotation(path: Path, annotations: Sequence[RepeatAnnotation]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ANNOTATION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for annotation in annotations:
            writer.writerow(
                {
                    "family_id": annotation.family_id,
                    "monomer_length": annotation.monomer_length,
                    "best_known_id": annotation.best_known_id,
                    "best_known_length": annotation.best_known_length,
                    "best_orientation": annotation.best_orientation,
                    "shared_kmer_fraction": f"{annotation.shared_kmer_fraction:.6f}",
                    "jaccard": f"{annotation.jaccard:.6f}",
                    "dice": f"{annotation.dice:.6f}",
                    "containment_discovered_in_known": f"{annotation.containment_discovered_in_known:.6f}",
                    "containment_known_in_discovered": f"{annotation.containment_known_in_discovered:.6f}",
                    "local_identity": f"{annotation.local_identity:.6f}",
                    "local_overlap_bp": annotation.local_overlap_bp,
                    "annotation_status": annotation.annotation_status,
                    "notes": annotation.notes,
                }
            )


def annotate_repeat_catalog(
    catalog_path: Path,
    known_path: Path,
    output_path: Path,
    *,
    kmer_size: int = 11,
) -> list[RepeatAnnotation]:
    catalog = read_discovered_catalog(catalog_path)
    known = read_known_repeats(known_path)
    annotations = annotate_discovered_repeats(catalog, known, kmer_size=kmer_size)
    write_repeat_annotation(output_path, annotations)
    return annotations
