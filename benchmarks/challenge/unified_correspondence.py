"""Independent native-catalogue correspondence for unified benchmark scoring.

This module only labels sequence correspondence.  It neither changes a native
motif nor estimates abundance.  It intentionally does not import TandemX
clustering or family-assignment code.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import edlib


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class CatalogueCorrespondence:
    """All truth families meeting the declared identity threshold for one motif."""

    status: str
    matches: tuple[str, ...]
    identities: tuple[tuple[str, float], ...]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _cyclic_orientations(sequence: str) -> tuple[str, ...]:
    """Return unique rotations in both orientations without external tool code."""
    variants = {
        oriented[offset:] + oriented[:offset]
        for oriented in (sequence, _reverse_complement(sequence))
        for offset in range(len(sequence))
    }
    return tuple(sorted(variants))


def _global_identity(first: str, second: str) -> float:
    """Needleman--Wunsch edit identity from independent edlib global alignment."""
    distance = edlib.align(first, second, mode="NW", task="distance")["editDistance"]
    if distance < 0:
        raise ValueError("edlib did not return a global edit distance")
    return 1 - distance / max(len(first), len(second))


def _repeat_identity(first: str, second: str, multiple: int) -> float:
    """Compare the longer sequence to all rotations/strands of the shorter unit."""
    shorter, longer = (first, second) if len(first) <= len(second) else (second, first)
    return max(_global_identity(variant * multiple, longer) for variant in _cyclic_orientations(shorter))


def _integer_multiple(first: str, second: str) -> int | None:
    shorter = min(len(first), len(second))
    longer = max(len(first), len(second))
    multiple = round(longer / shorter)
    if multiple < 1 or abs(longer - shorter * multiple) / longer > 0.10:
        return None
    return multiple


def _detectable_mixed_composite(native: str, truth: Mapping[str, str]) -> bool:
    """Flag exact multi-truth composition rather than force it to one family.

    This deliberately conservative check catches native HOR motifs that contain
    complete units from at least two truth families.  It is not a decomposition
    algorithm; motifs not caught here remain ``unmatched``.
    """
    contained = {
        family_id
        for family_id, sequence in truth.items()
        if any(variant in native for variant in _cyclic_orientations(sequence))
    }
    return len(contained) >= 2


def _validate_catalogue(catalogue: Mapping[str, str], label: str) -> None:
    for identifier, sequence in catalogue.items():
        if not isinstance(identifier, str) or not identifier:
            raise ValueError(f"{label} catalogue identifiers must be nonempty strings")
        if not isinstance(sequence, str) or not sequence or set(sequence.upper()) - set("ACGTN"):
            raise ValueError(f"Invalid {label} catalogue sequence for {identifier}")


def match_native_catalogue(
    native: Mapping[str, str],
    truth: Mapping[str, str],
    threshold: float = 0.9,
    max_multiple: int = 64,
) -> dict[str, CatalogueCorrespondence]:
    """Map native motifs to every threshold-passing truth family.

    A pure unit may be compared with an integer 1--``max_multiple`` repeat of
    the shorter sequence when their expected lengths agree within 10 percent.
    If two or more truth families pass, the result is ``ambiguous`` and no
    best-scoring family is selected.  A detectable mixed-family HOR is reported
    as ``unmatched_or_composite`` rather than assigned a spurious exact family.
    """
    if not isinstance(threshold, (int, float)) or isinstance(threshold, bool) or not 0 <= threshold <= 1:
        raise ValueError("threshold must be a finite fraction in [0,1]")
    if not isinstance(max_multiple, int) or isinstance(max_multiple, bool) or max_multiple < 1:
        raise ValueError("max_multiple must be a positive integer")
    _validate_catalogue(native, "native")
    _validate_catalogue(truth, "truth")

    result: dict[str, CatalogueCorrespondence] = {}
    for native_id, native_sequence in native.items():
        matches: list[tuple[str, float]] = []
        beyond_multiple = False
        for truth_id, truth_sequence in truth.items():
            multiple = _integer_multiple(native_sequence, truth_sequence)
            if multiple is None:
                continue
            if multiple > max_multiple:
                # Evaluate only the immediate boundary excess to distinguish an
                # otherwise matching 65x unit from an unrelated length mismatch.
                if multiple == max_multiple + 1 and _repeat_identity(native_sequence, truth_sequence, multiple) >= threshold:
                    beyond_multiple = True
                continue
            identity = _repeat_identity(native_sequence, truth_sequence, multiple)
            if identity >= threshold:
                matches.append((truth_id, identity))
        matches.sort(key=lambda item: item[0])
        if len(matches) == 1:
            status = "unique"
        elif len(matches) > 1:
            status = "ambiguous"
        elif beyond_multiple:
            status = "out_of_scope"
        elif _detectable_mixed_composite(native_sequence, truth):
            status = "unmatched_or_composite"
        else:
            status = "unmatched"
        result[native_id] = CatalogueCorrespondence(status, tuple(item[0] for item in matches), tuple(matches))
    return result
