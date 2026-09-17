"""Bounded local monomer segmentation and read--assembly label-path audit.

Inputs are pretrimmed array intervals with the same independently verified
left/right locus anchors. This module neither establishes anchoring nor
calibrates biological error probabilities. All scores are edit costs, not
posterior probabilities. The route is research-only and intentionally does not
touch TandemX's public commands.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import floor
from typing import Mapping

import edlib


_RC = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Copy:
    start: int
    end: int
    label: str
    orientation: str
    edit_cost: int
    label_margin: int | None


@dataclass(frozen=True)
class Decomposition:
    state: str
    copies: tuple[Copy, ...]
    score: int | None
    alternative_score: int | None
    reason: str


@dataclass(frozen=True)
class AuditResult:
    state: str
    reason: str
    assembly: Decomposition
    read_paths: tuple[tuple[str, tuple[str, ...]], ...]
    support_count: int
    discordant_count: int
    ambiguous_count: int
    candidate_event: str | None
    label_edit_distance: int | None
    breakpoint_label_interval: tuple[int, int] | None
    score_is_calibrated: bool = False


@dataclass(frozen=True)
class _Path:
    cost: int
    calls: tuple[Copy, ...]
    boundary_tied: bool = False

    @property
    def signature(self) -> tuple[tuple[str, str], ...]:
        return tuple((copy.label, copy.orientation) for copy in self.calls)


def _distance(left: str, right: str, limit: int) -> int | None:
    value = edlib.align(left, right, mode="NW", task="distance", k=limit)["editDistance"]
    return None if value < 0 else int(value)


def _validate(sequence: str, monomers: Mapping[str, str]) -> dict[str, str]:
    if len(sequence) > 4096 or set(sequence.upper()) - set("ACGT"):
        raise ValueError("Expected 0..4096 unambiguous array bases")
    if not monomers or len(monomers) > 16:
        raise ValueError("Expected 1..16 supplied candidate monomers")
    motifs = {str(label): str(motif).upper() for label, motif in monomers.items()}
    if any(not label or not motif or not 4 <= len(motif) <= 300 or
           set(motif) - set("ACGT") for label, motif in motifs.items()):
        raise ValueError("Monomers need labels and 4..300 ACGT bases")
    if len(set(motifs.values())) != len(motifs):
        raise ValueError("Duplicate candidate templates are non-identifiable")
    return motifs


def _insert(paths: list[_Path], proposed: _Path) -> None:
    for index, old in enumerate(paths):
        if old.signature == proposed.signature:
            if proposed.cost < old.cost:
                paths[index] = proposed
            elif proposed.cost == old.cost and proposed.calls != old.calls:
                winner = min((old, proposed), key=lambda x: tuple(
                    (c.start, c.end, c.label, c.orientation) for c in x.calls))
                paths[index] = _Path(winner.cost, winner.calls, True)
            break
    else:
        paths.append(proposed)
    paths.sort(key=lambda path: (path.cost, path.signature))
    del paths[3:]


def decompose(
    sequence: str,
    monomers: Mapping[str, str],
    *,
    max_edit_fraction: float = 0.15,
    max_indel_fraction: float = 0.08,
    minimum_label_margin: int = 2,
    max_candidate_alignments: int = 250_000,
) -> Decomposition:
    """Global tiling of a pretrimmed array by candidate monomer templates.

    Each segment uses an independent bounded Levenshtein emission. The best
    three distinct label paths are retained per sequence endpoint. A label
    whose best and next template costs differ by < ``minimum_label_margin``
    is unresolved; repeating the same error across reads cannot resolve it.
    This is an edit-cost DP, not a calibrated HMM.
    """
    motifs = _validate(sequence, monomers)
    if not sequence:
        return Decomposition("RESOLVED", (), 0, None, "empty_anchored_array_interval")
    if not 0 <= max_edit_fraction <= 0.5 or not 0 <= max_indel_fraction <= 0.3:
        raise ValueError("Invalid edit or indel fraction")
    if minimum_label_margin < 1 or max_candidate_alignments < 1:
        raise ValueError("Invalid search limit or label margin")
    sequence = sequence.upper()
    templates = [(label, strand, motif if strand == "+" else motif.translate(_RC)[::-1])
                 for label, motif in sorted(motifs.items()) for strand in ("+", "-")]
    paths: list[list[_Path]] = [[] for _ in range(len(sequence) + 1)]
    paths[0] = [_Path(0, ())]
    comparisons = 0
    for start, prefixes in enumerate(paths[:-1]):
        if not prefixes:
            continue
        for label, strand, motif in templates:
            limit = max(1, floor(len(motif) * max_edit_fraction))
            wiggle = max(1, floor(len(motif) * max_indel_fraction))
            for width in range(max(1, len(motif) - wiggle), len(motif) + wiggle + 1):
                end = start + width
                if end > len(sequence):
                    continue
                comparisons += 1
                if comparisons > max_candidate_alignments:
                    return Decomposition("AMBIGUOUS", (), None, None,
                                         "candidate_alignment_budget_exceeded")
                cost = _distance(motif, sequence[start:end], limit)
                if cost is None:
                    continue
                call = Copy(start, end, label, strand, cost, None)
                for prefix in prefixes:
                    _insert(paths[end], _Path(prefix.cost + cost, prefix.calls + (call,),
                                              prefix.boundary_tied))
    if not paths[-1]:
        return Decomposition("AMBIGUOUS", (), None, None, "no_full_monomer_tiling")
    best = paths[-1][0]
    runner = paths[-1][1].cost if len(paths[-1]) > 1 else None
    assessed: list[Copy] = []
    for call in best.calls:
        competing = []
        width = call.end - call.start
        for other_label, other_strand, other_motif in templates:
            if (other_label, other_strand) == (call.label, call.orientation):
                continue
            if abs(len(other_motif) - width) > max(1, floor(
                    len(other_motif) * max_indel_fraction)):
                continue
            comparisons += 1
            if comparisons > max_candidate_alignments:
                return Decomposition("AMBIGUOUS", (), None, None,
                                     "candidate_alignment_budget_exceeded")
            alternate = _distance(other_motif, sequence[call.start:call.end],
                                  call.edit_cost + minimum_label_margin - 1)
            if alternate is not None:
                competing.append(alternate)
        margin = min(competing) - call.edit_cost if competing else None
        assessed.append(Copy(call.start, call.end, call.label, call.orientation,
                             call.edit_cost, margin))
    calls = tuple(assessed)
    if any(c.label_margin is not None and c.label_margin < minimum_label_margin
           for c in calls):
        return Decomposition("AMBIGUOUS", calls, best.cost, runner,
                             "locally_nonidentifiable_monomer_label_or_orientation")
    if runner is not None and runner - best.cost < minimum_label_margin:
        return Decomposition("AMBIGUOUS", calls, best.cost, runner,
                             "competing_global_label_paths")
    return Decomposition("RESOLVED", calls, best.cost, runner,
                         "boundary_tied" if best.boundary_tied else "unique_label_path")


def _alignment_breakpoints(
    assembly: tuple[tuple[str, str], ...],
    read: tuple[tuple[str, str], ...],
) -> tuple[int, tuple[int, int] | None]:
    """Edit distance plus envelope of every optimal changed assembly index."""
    n, m = len(assembly), len(read)
    forward = [[0] * (m + 1) for _ in range(n + 1)]
    reverse = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        forward[i][0] = i
    for j in range(m + 1):
        forward[0][j] = j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            forward[i][j] = min(forward[i - 1][j] + 1,
                                forward[i][j - 1] + 1,
                                forward[i - 1][j - 1] + (assembly[i - 1] != read[j - 1]))
    for i in range(n - 1, -1, -1):
        reverse[i][m] = n - i
    for j in range(m - 1, -1, -1):
        reverse[n][j] = m - j
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            reverse[i][j] = min(reverse[i + 1][j] + 1,
                                reverse[i][j + 1] + 1,
                                reverse[i + 1][j + 1] + (assembly[i] != read[j]))
    optimum = forward[n][m]
    if optimum == 0:
        return 0, None
    positions: list[int] = []
    for i in range(n + 1):
        for j in range(m + 1):
            before = forward[i][j]
            if i < n and before + 1 + reverse[i + 1][j] == optimum:
                positions.append(i)
            if j < m and before + 1 + reverse[i][j + 1] == optimum:
                positions.append(i)
            if i < n and j < m and assembly[i] != read[j] and \
                    before + 1 + reverse[i + 1][j + 1] == optimum:
                positions.append(i)
    if n == 0:
        return optimum, (0, 0)
    return optimum, (min(positions), max(positions) + 1)


def _event_type(assembly: tuple[tuple[str, str], ...],
                read: tuple[tuple[str, str], ...]) -> str:
    if len(read) > len(assembly):
        return "assembly_copy_loss_or_compression"
    if len(read) < len(assembly):
        return "assembly_extra_copies_or_read_truncation"
    if [x[0] for x in assembly] == [x[0] for x in read]:
        return "orientation_disagreement"
    return "equal_length_order_or_label_disagreement"


def audit(
    assembly_sequence: str,
    anchored_reads: Mapping[str, str],
    monomers: Mapping[str, str],
    *,
    pairing_status: str,
    min_informative_reads: int = 3,
    min_consensus_fraction: float = 0.75,
    min_resolved_fraction: float = 0.75,
    minimum_label_margin: int = 2,
) -> AuditResult:
    """Audit a single same-locus array using anchored independent molecules.

    ``pairing_status='synthetic'`` only admits controlled development cases;
    biological interpretation requires ``'verified'`` external flank,
    molecule and haplotype checks. ``'unverified'`` explicitly abstains.
    """
    if min_informative_reads < 2 or not 0.5 < min_consensus_fraction <= 1 or \
            not 0.5 < min_resolved_fraction <= 1:
        raise ValueError("Invalid read support settings")
    if pairing_status not in {"verified", "synthetic", "synthetic_simulated",
                              "unverified", "absent"}:
        raise ValueError("Unknown read pairing status")
    assembly = decompose(assembly_sequence, monomers,
                         minimum_label_margin=minimum_label_margin)
    empty = ((), 0, 0, 0, None, None, None)
    if assembly.state != "RESOLVED":
        return AuditResult("AMBIGUOUS", "assembly_" + assembly.reason, assembly, *empty)
    if pairing_status in {"unverified", "absent"}:
        return AuditResult("AMBIGUOUS", "read_pairing_not_verified", assembly, *empty)
    if not anchored_reads:
        return AuditResult("INSUFFICIENT_READ_SUPPORT", "no_anchored_reads", assembly, *empty)
    assembly_path = tuple((c.label, c.orientation) for c in assembly.copies)
    read_paths: list[tuple[str, tuple[str, ...]]] = []
    groups: dict[tuple[tuple[str, str], ...], list[str]] = {}
    ambiguous = 0
    for read_id, sequence in sorted(anchored_reads.items()):
        read = decompose(sequence, monomers,
                         minimum_label_margin=minimum_label_margin)
        if read.state != "RESOLVED":
            ambiguous += 1
            continue
        path = tuple((c.label, c.orientation) for c in read.copies)
        groups.setdefault(path, []).append(read_id)
        read_paths.append((read_id, tuple(c.label + c.orientation for c in read.copies)))
    informative = len(read_paths)
    support = len(groups.get(assembly_path, []))
    discordant = informative - support
    prefix = (tuple(read_paths), support, discordant, ambiguous)
    if informative < min_informative_reads:
        return AuditResult("INSUFFICIENT_READ_SUPPORT", "too_few_resolved_anchored_reads",
                           assembly, *prefix, None, None, None)
    if informative / len(anchored_reads) < min_resolved_fraction:
        return AuditResult("AMBIGUOUS", "too_many_unresolved_read_paths",
                           assembly, *prefix, None, None, None)
    dominant_path, dominant_ids = max(groups.items(), key=lambda item: (len(item[1]), item[0]))
    if len(dominant_ids) / informative < min_consensus_fraction:
        return AuditResult("AMBIGUOUS", "mixed_read_architectures_or_haplotypes",
                           assembly, *prefix, None, None, None)
    if dominant_path == assembly_path:
        if discordant:
            return AuditResult("AMBIGUOUS", "minority_discordant_reads",
                               assembly, *prefix, None, None, None)
        return AuditResult("SUPPORTED", "all_resolved_reads_match_assembly_label_path",
                           assembly, *prefix, None, 0, None)
    # A minority assembly-like path can arise from alleles, mapping, or read
    # errors. Until separately modeled, it is not called an assembly defect.
    if support or len(groups) > 1:
        return AuditResult("AMBIGUOUS", "minority_path_or_haplotype_conflict",
                           assembly, *prefix, None, None, None)
    distance, bp = _alignment_breakpoints(assembly_path, dominant_path)
    return AuditResult("DISCORDANT", "consistent_distinct_read_label_path;"
                       "not_assembly_error_truth", assembly, *prefix,
                       _event_type(assembly_path, dominant_path), distance, bp)
