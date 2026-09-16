"""Small, explicit monomer-order prototype for pre-trimmed array intervals.

This is a development experiment, not a genome scanner or a production HOR
caller. Full-interval calls require a provided candidate monomer catalogue.
Copy-count discordance is evaluated only for reads independently established
to span both unique flanks of the same locus.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class MonomerCopy:
    start: int
    end: int
    label: str
    orientation: str
    edit_distance: int


@dataclass(frozen=True)
class Architecture:
    copies: tuple[MonomerCopy, ...]
    cyclic_unit: tuple[str, ...]
    variant_copy_indices: tuple[int, ...]
    complete_unit_count: int
    template_ambiguous: bool
    status: str


@dataclass(frozen=True)
class Discordance:
    status: str
    assembly: Architecture
    read_support: tuple[str, ...]
    discordant_reads: tuple[str, ...]
    warning: str


def _rotations(sequence: str) -> tuple[str, ...]:
    return tuple(sorted({sequence[i:] + sequence[:i] for i in range(len(sequence))}))


def _edit_distance(left: str, right: str, limit: int) -> int:
    if abs(len(left) - len(right)) > limit:
        return limit + 1
    previous = list(range(len(right) + 1))
    for i, base in enumerate(left, 1):
        current = [i]
        for j, other in enumerate(right, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1,
                               previous[j - 1] + (base != other)))
        previous = current
    return previous[-1]


def _validate(sequence: str, monomers: Mapping[str, str]) -> dict[str, str]:
    sequence = sequence.upper()
    if not sequence or len(sequence) > 4096 or set(sequence) - set("ACGT"):
        raise ValueError("Array interval must contain 1..4096 unambiguous bases")
    if not monomers or len(monomers) > 16:
        raise ValueError("Provide 1..16 candidate monomers")
    normalized = {label: motif.upper() for label, motif in monomers.items()}
    if any(not label or not motif or len(motif) > 64 or set(motif) - set("ACGT")
           for label, motif in normalized.items()):
        raise ValueError("Monomer labels and 1..64 bp ACGT motifs are required")
    canonical = [min(_rotations(motif) + _rotations(motif.translate(_COMPLEMENT)[::-1]))
                 for motif in normalized.values()]
    if len(set(canonical)) != len(canonical):
        raise ValueError("Cyclic/reverse-complement-equivalent monomers are not identifiable")
    return normalized


def _period(labels: tuple[str, ...]) -> tuple[tuple[str, ...], tuple[int, ...], int, bool]:
    if len(labels) < 2:
        return (), (), 0, False
    options: list[tuple[int, int, tuple[str, ...], tuple[int, ...], bool]] = []
    for width in range(1, len(labels) // 2 + 1):
        if len(labels) % width:
            continue
        unit_list: list[str] = []
        tied = False
        for phase in range(width):
            counts = Counter(labels[phase::width])
            maximum = max(counts.values())
            modal = sorted(label for label, count in counts.items() if count == maximum)
            unit_list.append(modal[0])
            tied |= len(modal) > 1
        unit = tuple(unit_list)
        variants = tuple(i for i, label in enumerate(labels)
                         if label != unit[i % width])
        if len(variants) / len(labels) <= 0.25:
            options.append((len(variants), width, unit, variants, tied))
    if not options:
        return (), (), 0, False
    # Identical mismatch counts prefer the smaller fundamental unit. A tied
    # modal column has no supported ancestral/standard variant direction.
    _, width, unit, variants, tied = min(options)
    canonical_unit = min(unit[i:] + unit[:i] for i in range(width))
    return canonical_unit, variants, len(labels) // width, tied


def infer_architecture(
    sequence: str,
    monomers: Mapping[str, str],
    *,
    max_edit_fraction: float = 0.20,
) -> Architecture:
    """Segment one pre-trimmed array, then infer a repeated label unit.

    Exhaustive cyclic/strand matching is deliberately limited to short motifs
    and short intervals. Ties among distinct labels are marked unresolved.
    """
    normalized = _validate(sequence, monomers)
    if not 0 <= max_edit_fraction <= 0.5:
        raise ValueError("max_edit_fraction must be between 0 and 0.5")
    sequence = sequence.upper()
    n = len(sequence)
    estimated_cells = sum(
        n * 2 * len(_rotations(motif)) *
        (2 * int(len(motif) * max_edit_fraction) + 1) *
        len(motif) * (len(motif) + int(len(motif) * max_edit_fraction))
        for motif in normalized.values()
    )
    if estimated_cells > 40_000_000:
        raise ValueError("Prototype work estimate exceeds 40 million DP cells")
    # Each path entry is (total errors, number of copies, calls, ambiguous).
    paths: list[tuple[int, int, tuple[MonomerCopy, ...], bool] | None] = [None] * (n + 1)
    paths[0] = (0, 0, (), False)
    candidates: list[tuple[str, str, str, int]] = []
    for label, motif in sorted(normalized.items()):
        limit = int(len(motif) * max_edit_fraction)
        strand_identifiable = not (set(_rotations(motif)) &
                                   set(_rotations(motif.translate(_COMPLEMENT)[::-1])))
        for orientation, template in (("+", motif),
                                      ("-", motif.translate(_COMPLEMENT)[::-1])):
            for rotation in _rotations(template):
                candidates.append((label, orientation if strand_identifiable else "?",
                                   rotation, limit))
    for start in range(n):
        prior = paths[start]
        if prior is None:
            continue
        for label, orientation, template, limit in candidates:
            for width in range(max(1, len(template) - limit), len(template) + limit + 1):
                end = start + width
                if end > n:
                    continue
                distance = _edit_distance(template, sequence[start:end], limit)
                if distance > limit:
                    continue
                call = MonomerCopy(start, end, label, orientation, distance)
                proposal = (prior[0] + distance, prior[1] + 1,
                            prior[2] + (call,), prior[3])
                incumbent = paths[end]
                if incumbent is None or proposal[:2] < incumbent[:2]:
                    paths[end] = proposal
                elif proposal[:2] == incumbent[:2]:
                    # Alternative phases of the same label are harmless; a
                    # distinct label/order or boundary makes the call ambiguous.
                    different = tuple((c.start, c.end, c.label) for c in proposal[2]) != tuple(
                        (c.start, c.end, c.label) for c in incumbent[2])
                    winner = min((proposal, incumbent), key=lambda p: tuple(
                        (c.label, c.start, c.end, c.orientation) for c in p[2]))
                    paths[end] = (winner[0], winner[1], winner[2],
                                  proposal[3] or incumbent[3] or different)
    best = paths[n]
    if best is None:
        return Architecture((), (), (), 0, False, "unresolved_no_full_decomposition")
    labels = tuple(copy.label for copy in best[2])
    unit, variants, count, template_ambiguous = _period(labels)
    status = ("unresolved_ambiguous_decomposition" if best[3] else
              "candidate_periodic_template_ambiguous" if unit and template_ambiguous else
              "candidate_periodic_with_variants" if variants else
              "candidate_periodic" if unit else "unresolved_no_repeated_unit")
    return Architecture(best[2], unit, variants, count, template_ambiguous, status)


def compare_anchored_reads(
    assembly_sequence: str,
    anchored_read_intervals: Mapping[str, str],
    monomers: Mapping[str, str],
    *,
    min_supporting_reads: int = 2,
) -> Discordance:
    """Compare pre-verified, independently anchored read intervals to assembly.

    This function does not prove that supplied reads really span both unique
    flanks. That provenance must be verified before using its output as evidence.
    """
    if min_supporting_reads < 2:
        raise ValueError("At least two independent molecules are required")
    assembly = infer_architecture(assembly_sequence, monomers)
    if not assembly.cyclic_unit or assembly.status == "unresolved_ambiguous_decomposition":
        return Discordance("unresolved", assembly, (), (), "assembly_architecture_unresolved")
    supporting: list[str] = []
    discordant: list[str] = []
    for read_id, sequence in sorted(anchored_read_intervals.items()):
        read = infer_architecture(sequence, monomers)
        if not read.cyclic_unit or read.status == "unresolved_ambiguous_decomposition":
            continue
        supporting.append(read_id)
        assembly_order = tuple((copy.label, copy.orientation) for copy in assembly.copies)
        read_order = tuple((copy.label, copy.orientation) for copy in read.copies)
        # The smallest HOR unit is cyclically normalized, but the supplied
        # intervals have the same verified left unique flank. Moving a
        # variant/order relative to that flank is a real locus difference.
        if read_order != assembly_order:
            discordant.append(read_id)
    if len(supporting) < min_supporting_reads:
        return Discordance("unresolved", assembly, tuple(supporting),
                           tuple(discordant), "insufficient_independent_anchored_reads")
    if len(discordant) >= min_supporting_reads and len(discordant) / len(supporting) > 2 / 3:
        return Discordance("candidate_discordance", assembly, tuple(supporting),
                           tuple(discordant), "requires_external_flank_and_molecule_qc")
    if discordant:
        return Discordance("unresolved_mixed_molecules", assembly, tuple(supporting),
                           tuple(discordant), "possible_haplotype_or_alignment_mixture")
    return Discordance("no_supported_discordance", assembly, tuple(supporting),
                       tuple(discordant), "does_not_prove_assembly_correctness")
