"""Post hoc similarity checks between discovered and known repeat monomers."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from tandemx.io.sequences import read_sequence_records
from tandemx.simulate.toy import reverse_complement


@dataclass(frozen=True)
class RepeatSequence:
    identifier: str
    sequence: str


@dataclass(frozen=True)
class KnownRepeatMatch:
    known_repeat_id: str
    known_repeat_length: int
    best_family_id: str
    best_monomer_length: int
    similarity_score: float
    shared_kmer_fraction: float
    orientation: str
    interpretation: str


def parse_family_id(description: str, fallback: str) -> str:
    for part in description.split(";"):
        if part.startswith("family_id="):
            return part.split("=", 1)[1]
    return fallback


def read_repeat_sequences(path: Path, *, catalog: bool) -> list[RepeatSequence]:
    records = []
    for record in read_sequence_records(path):
        identifier = parse_family_id(record.description, record.id) if catalog else record.id
        records.append(RepeatSequence(identifier=identifier, sequence=record.sequence))
    return records


def sequence_kmers(sequence: str, k: int) -> set[str]:
    return {
        sequence[index : index + k]
        for index in range(len(sequence) - k + 1)
        if "N" not in sequence[index : index + k]
    }


def compare_kmer_sets(known: set[str], candidate: set[str]) -> tuple[float, float]:
    if not known or not candidate:
        return 0.0, 0.0
    shared = len(known.intersection(candidate))
    similarity = 2.0 * shared / (len(known) + len(candidate))
    shared_fraction = shared / len(known)
    return similarity, shared_fraction


def interpretation_for_score(score: float) -> str:
    if score >= 0.8:
        return "strong_post_hoc_match"
    if score >= 0.5:
        return "possible_post_hoc_match"
    if score > 0:
        return "weak_post_hoc_match"
    return "no_kmer_match"


def match_known_repeats(
    catalog: Sequence[RepeatSequence],
    known_repeats: Iterable[RepeatSequence],
    *,
    kmer_size: int = 11,
) -> list[KnownRepeatMatch]:
    if kmer_size <= 0:
        raise ValueError("--kmer-size must be positive")
    if not catalog:
        raise ValueError("Discovered catalog contains no monomers")

    matches = []
    for known in known_repeats:
        best: tuple[float, float, RepeatSequence, str] | None = None
        for monomer in catalog:
            effective_k = min(kmer_size, len(known.sequence), len(monomer.sequence))
            known_kmers = sequence_kmers(known.sequence, effective_k)
            for orientation, sequence in (
                ("forward", monomer.sequence),
                ("reverse", reverse_complement(monomer.sequence)),
            ):
                similarity, shared_fraction = compare_kmer_sets(
                    known_kmers,
                    sequence_kmers(sequence, effective_k),
                )
                candidate = (similarity, shared_fraction, monomer, orientation)
                if best is None or candidate[:2] > best[:2]:
                    best = candidate
        assert best is not None
        similarity, shared_fraction, monomer, orientation = best
        matches.append(
            KnownRepeatMatch(
                known_repeat_id=known.identifier,
                known_repeat_length=len(known.sequence),
                best_family_id=monomer.identifier,
                best_monomer_length=len(monomer.sequence),
                similarity_score=similarity,
                shared_kmer_fraction=shared_fraction,
                orientation=orientation,
                interpretation=interpretation_for_score(similarity),
            )
        )
    return matches


def write_known_repeat_matches(path: Path, matches: Sequence[KnownRepeatMatch]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "known_repeat_id",
        "known_repeat_length",
        "best_family_id",
        "best_monomer_length",
        "similarity_score",
        "shared_kmer_fraction",
        "orientation",
        "interpretation",
    )
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for match in matches:
            writer.writerow(
                {
                    "known_repeat_id": match.known_repeat_id,
                    "known_repeat_length": match.known_repeat_length,
                    "best_family_id": match.best_family_id,
                    "best_monomer_length": match.best_monomer_length,
                    "similarity_score": f"{match.similarity_score:.6f}",
                    "shared_kmer_fraction": f"{match.shared_kmer_fraction:.6f}",
                    "orientation": match.orientation,
                    "interpretation": match.interpretation,
                }
            )


def check_known_repeats(
    catalog_path: Path,
    known_path: Path,
    output_path: Path,
    *,
    kmer_size: int = 11,
) -> list[KnownRepeatMatch]:
    catalog = read_repeat_sequences(catalog_path, catalog=True)
    known = read_repeat_sequences(known_path, catalog=False)
    matches = match_known_repeats(catalog, known, kmer_size=kmer_size)
    write_known_repeat_matches(output_path, matches)
    return matches
