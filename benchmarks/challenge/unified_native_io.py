"""Strict native-output readers for unified benchmark scoring.

``candidate_reads.tsv`` fields: ``read_id``, ``candidate_id``, ``read_start``,
``read_end`` and ``period_bp``.  ``monomer_membership.tsv`` fields: ``read_id``,
``candidate_id`` and ``family_id``.  PAF uses its standard mandatory twelve
columns; query coordinates are interpreted as 0-based half-open read intervals.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from benchmarks.challenge.schema import ArrayRecord, iter_table
from tandemx.quantify.mvp import parse_family_id, read_fasta
from tandemx.io.sequences import SequenceFormatError


_CANDIDATE_FIELDS = {"read_id", "candidate_id", "read_start", "read_end", "period_bp"}
_MEMBERSHIP_FIELDS = {"read_id", "candidate_id", "family_id", "status"}


def load_tandemx_catalogue(monomers: Path, *, allow_empty: bool = False) -> dict[str, str]:
    """Read machine-header ``family_id`` records from TandemX ``monomers.fa``."""
    catalogue: dict[str, str] = {}
    try:
        records = read_fasta(monomers)
        for record in records:
            if not any(part.startswith("family_id=") for part in record.description.split(";")):
                raise ValueError(f"TandemX monomer header lacks family_id: {record.description}")
            family_id = parse_family_id(record.description)
            if family_id in catalogue:
                raise ValueError(f"Duplicate TandemX family_id: {family_id}")
            catalogue[family_id] = record.sequence
    except SequenceFormatError:
        if not allow_empty or monomers.stat().st_size:
            raise
    if not catalogue and not allow_empty:
        raise ValueError("TandemX monomer catalogue contains no records")
    return catalogue


def load_tandemx_arrays(discoverdir: Path) -> list[ArrayRecord]:
    """Strictly join TandemX candidates to their assigned family monomers.

    Each candidate key ``(read_id, candidate_id)`` must have exactly one
    membership row and that row's family must exist in ``monomers.fa``.  This
    reader deliberately rejects unresolved and below-support rows rather than
    assigning them a synthetic native family.
    """
    catalogue = load_tandemx_catalogue(discoverdir / "monomers.fa", allow_empty=True)
    membership: dict[tuple[str, str], tuple[str, str]] = {}
    for row in iter_table(discoverdir / "monomer_membership.tsv", _MEMBERSHIP_FIELDS):
        key = (row["read_id"], row["candidate_id"])
        if key in membership:
            raise ValueError(f"Multiple TandemX membership rows for {key}")
        family_id, status = row["family_id"], row["status"]
        if family_id == "NA" and status in {"unresolved_sequence", "below_minimum_support"}:
            membership[key] = (f"unassigned_candidate:{key[1]}@{key[0]}", "")
            continue
        if family_id not in catalogue:
            raise ValueError(f"Unknown TandemX membership family for {key}: {family_id}")
        membership[key] = (family_id, catalogue[family_id])

    arrays: list[ArrayRecord] = []
    seen: set[tuple[str, str]] = set()
    for row in iter_table(discoverdir / "candidate_reads.tsv", _CANDIDATE_FIELDS):
        key = (row["read_id"], row["candidate_id"])
        if key in seen:
            raise ValueError(f"Duplicate TandemX candidate row for {key}")
        seen.add(key)
        if key not in membership:
            raise ValueError(f"Unmatched TandemX candidate row for {key}")
        family_id, sequence = membership[key]
        arrays.append(
            ArrayRecord(
                read_id=row["read_id"],
                start=_integer(row["read_start"], "candidate start"),
                end=_integer(row["read_end"], "candidate end"),
                period=_integer(row["period_bp"], "candidate period"),
                sequence=sequence,
                family_id=family_id,
            )
        )
    if set(membership) != seen:
        extra = sorted(set(membership) - seen)[0]
        raise ValueError(f"Unmatched TandemX membership row for {extra}")
    return arrays


def mapping_paf_to_unique_arrays(
    paf: Path,
    catalogue: dict[str, str],
    minblock: int = 100,
    identity: float = 0.9,
) -> tuple[list[ArrayRecord], int]:
    """Normalize all qualifying PAF maps without truth-assisted tie breaking.

    Primary and secondary rows are retained if their alignment block and
    identity pass the supplied gates.  Query intervals first union within each
    native family. Any query base covered by two or more native families is
    excluded from returned arrays and counted once in ``ambiguous_bp``.
    """
    if not isinstance(minblock, int) or isinstance(minblock, bool) or minblock < 1:
        raise ValueError("minblock must be a positive integer")
    if not isinstance(identity, (int, float)) or isinstance(identity, bool) or not 0 <= identity <= 1:
        raise ValueError("identity must be in [0,1]")
    if not catalogue or any(not key or not value for key, value in catalogue.items()):
        raise ValueError("catalogue must contain nonempty native family IDs and sequences")

    by_read_family: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF row {paf}:{line_number}")
            query, qlength, qstart, qend = fields[:4]
            target, tlength, tstart, tend, matches, block = fields[5:11]
            qlength_i = _integer(qlength, "PAF query length")
            qstart_i = _integer(qstart, "PAF query start")
            qend_i = _integer(qend, "PAF query end")
            tlength_i = _integer(tlength, "PAF target length")
            tstart_i = _integer(tstart, "PAF target start")
            tend_i = _integer(tend, "PAF target end")
            matches_i = _integer(matches, "PAF matches")
            block_i = _integer(block, "PAF alignment block")
            if not query or not 0 <= qstart_i < qend_i <= qlength_i:
                raise ValueError(f"Invalid PAF query coordinates {paf}:{line_number}")
            if target not in catalogue:
                raise ValueError(f"Unknown PAF native family {target!r} at {paf}:{line_number}")
            if not 0 <= tstart_i < tend_i <= tlength_i or not 0 <= matches_i <= block_i:
                raise ValueError(f"Invalid PAF target coordinates/counts {paf}:{line_number}")
            if block_i >= minblock and matches_i / block_i >= identity:
                by_read_family[query][target].append((qstart_i, qend_i))

    arrays: list[ArrayRecord] = []
    ambiguous_bp = 0
    for read_id in sorted(by_read_family):
        merged = {family_id: _union(intervals) for family_id, intervals in by_read_family[read_id].items()}
        events: dict[int, list[tuple[str, int]]] = defaultdict(list)
        for family_id, intervals in merged.items():
            for start, end in intervals:
                events[start].append((family_id, 1))
                events[end].append((family_id, -1))
        active: set[str] = set()
        previous: int | None = None
        for position in sorted(events):
            if previous is not None and previous < position:
                if len(active) == 1:
                    family_id = next(iter(active))
                    arrays.append(ArrayRecord(read_id, previous, position, len(catalogue[family_id]), catalogue[family_id], family_id))
                elif len(active) > 1:
                    ambiguous_bp += position - previous
            for family_id, direction in events[position]:
                if direction == 1:
                    active.add(family_id)
                else:
                    active.remove(family_id)
            previous = position
    return arrays, ambiguous_bp


def mapping_paf_to_all_arrays(
    paf: Path,
    catalogue: dict[str, str],
    minblock: int = 100,
    identity: float = 0.9,
) -> list[ArrayRecord]:
    """Return every qualifying primary or secondary PAF interval unfiltered.

    This separate view is for global precision scoring.  It intentionally keeps
    cross-family overlaps; callers needing non-overlapping abundance support use
    :func:`mapping_paf_to_unique_arrays` instead.
    """
    if not isinstance(minblock, int) or isinstance(minblock, bool) or minblock < 1:
        raise ValueError("minblock must be a positive integer")
    if not isinstance(identity, (int, float)) or isinstance(identity, bool) or not 0 <= identity <= 1:
        raise ValueError("identity must be in [0,1]")
    if not catalogue or any(not key or not value for key, value in catalogue.items()):
        raise ValueError("catalogue must contain nonempty native family IDs and sequences")
    arrays: list[ArrayRecord] = []
    with paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF row {paf}:{line_number}")
            query, qlength, qstart, qend = fields[:4]
            target, tlength, tstart, tend, matches, block = fields[5:11]
            qlength_i, qstart_i, qend_i = (_integer(qlength, "PAF query length"), _integer(qstart, "PAF query start"), _integer(qend, "PAF query end"))
            tlength_i, tstart_i, tend_i = (_integer(tlength, "PAF target length"), _integer(tstart, "PAF target start"), _integer(tend, "PAF target end"))
            matches_i, block_i = _integer(matches, "PAF matches"), _integer(block, "PAF alignment block")
            if not query or not 0 <= qstart_i < qend_i <= qlength_i:
                raise ValueError(f"Invalid PAF query coordinates {paf}:{line_number}")
            if target not in catalogue:
                raise ValueError(f"Unknown PAF native family {target!r} at {paf}:{line_number}")
            if not 0 <= tstart_i < tend_i <= tlength_i or not 0 <= matches_i <= block_i:
                raise ValueError(f"Invalid PAF target coordinates/counts {paf}:{line_number}")
            if block_i >= minblock and matches_i / block_i >= identity:
                arrays.append(ArrayRecord(query, qstart_i, qend_i, len(catalogue[target]), catalogue[target], target))
    return arrays


def _integer(value: str, label: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"Invalid {label}: {value!r}") from error


def _union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result
