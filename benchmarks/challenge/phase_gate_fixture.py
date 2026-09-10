"""Small, development-only accuracy-gate fixture.

Truth is attached while constructing the source genome and projected to sampled
reads.  This module does not call a discovery or abundance algorithm and does
not generate future holdout data.
"""

from __future__ import annotations

import hashlib
import random
from collections import defaultdict

from .simulate import DNA, random_dna

CONDITIONS = (
    "clean", "substitutions", "indels_small", "indels_long", "partial",
    "interruptions", "background_homology", "close_families",
    "heterogeneous", "abundance_low",
)


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _mutate_unit(unit: str, rng: random.Random, *, substitutions: float = 0.0,
                 indel: int = 0, insertion: bool = False) -> tuple[str, list[bool]]:
    out: list[str] = []
    mask: list[bool] = []
    for base in unit:
        indel_event = bool(indel and rng.random() < 0.01)
        if indel_event and rng.random() < 0.5:
            continue
        if substitutions and rng.random() < substitutions:
            base = rng.choice(DNA.replace(base, ""))
        out.append(base)
        mask.append(True)
        if indel_event:
            out.append(rng.choice(DNA))
            mask.append(True)
    return "".join(out), mask


def _array(unit: str, copies: int, rng: random.Random, condition: str,
           *, heterogeneous: bool = False) -> tuple[str, list[bool]]:
    pieces: list[str] = []
    masks: list[bool] = []
    for copy in range(copies):
        if condition == "substitutions":
            seq, mask = _mutate_unit(unit, rng, substitutions=0.01)
        elif condition == "indels_small":
            seq, mask = _mutate_unit(unit, rng, indel=1, insertion=True)
        elif condition == "indels_long":
            seq = unit
            mask = [True] * len(unit)
            size = rng.randint(1, 20)
            at = rng.randrange(len(seq))
            if copy % 2:
                seq = seq[:at] + random_dna(rng, size) + seq[at:]
                mask = mask[:at] + [True] * size + mask[at:]
            else:
                seq = seq[:at] + seq[at + size:]
                mask = mask[:at] + mask[at + size:]
        elif heterogeneous:
            seq, mask = _mutate_unit(unit, rng, substitutions=0.02)
        else:
            seq, mask = unit, [True] * len(unit)
        pieces.append(seq)
        masks.extend(mask)
    sequence = "".join(pieces)
    if condition == "partial":
        trim = len(unit) // 2
        sequence, masks = sequence[trim:-trim], masks[trim:-trim]
    if condition == "interruptions":
        middle = len(sequence) // 2
        interruption = random_dna(rng, 300)
        sequence = sequence[:middle] + interruption + sequence[middle:]
        masks = masks[:middle] + [False] * 300 + masks[middle:]
    return sequence, masks


def _masked_intervals(mask: list[bool], family: str, offset: int = 0) -> list[tuple[int, int, str]]:
    rows: list[tuple[int, int, str]] = []
    start: int | None = None
    for i, is_repeat in enumerate(mask + [False]):
        if is_repeat and start is None:
            start = i
        elif not is_repeat and start is not None:
            rows.append((offset + start, offset + i, family))
            start = None
    return rows


def generate_case(seed_label: str, condition: str, coverage: int) -> dict:
    """Generate one deterministic development case with projected truth.

    ``seed_label`` must identify a development seed (for example ``dev-01``);
    labels containing ``holdout`` are refused so this helper cannot silently
    create a future holdout.  Coverage is restricted to the runner's 2/5/10x
    design points.
    """
    if "holdout" in seed_label.lower() or not seed_label.lower().startswith("dev"):
        raise ValueError("only dev seed labels are permitted")
    if condition not in CONDITIONS:
        raise ValueError(f"unknown condition: {condition}")
    if coverage not in {2, 5, 10}:
        raise ValueError("coverage must be one of 2, 5 or 10")
    seed = int.from_bytes(hashlib.sha256(seed_label.encode()).digest()[:8], "big")
    rng = random.Random(seed)
    f1 = random_dna(rng, 120)
    f2 = random_dna(rng, 120)
    if condition == "close_families":
        f2 = f1
        changed: list[int] = []
        while len(changed) < 10:
            pos = rng.randrange(120)
            if pos not in changed:
                changed.append(pos)
                f2 = f2[:pos] + rng.choice(DNA.replace(f2[pos], "")) + f2[pos + 1:]
    catalogue = {"F1": f1, "F2": f2}
    genome = list(random_dna(rng, 22_000))
    masks: dict[str, list[bool]] = {"F1": [False] * len(genome), "F2": [False] * len(genome)}

    def place(family: str, sequence: str, mask: list[bool], start: int) -> None:
        genome[start:start + len(sequence)] = sequence
        for i, value in enumerate(mask):
            masks[family][start + i] = value

    f1_copies = 2 if condition == "abundance_low" else 6
    seq1, mask1 = _array(f1, f1_copies, rng, condition,
                         heterogeneous=condition == "heterogeneous")
    seq2, mask2 = _array(f2, 12, rng, condition,
                         heterogeneous=condition == "heterogeneous")
    place("F1", seq1, mask1, 4_000)
    place("F2", seq2, mask2, 12_000)
    if condition == "background_homology":
        background = "".join(f1[:80] if i % 2 == 0 else random_dna(rng, 80) for i in range(30))
        genome[16_000:16_000 + len(background)] = background
    source = "".join(genome)
    source_intervals = {family: _masked_intervals(mask, family) for family, mask in masks.items()}
    source_truth_bp = {family: sum(end - start for start, end, _ in rows)
                       for family, rows in source_intervals.items()}

    reads: list[tuple[str, str, int]] = []
    truth_intervals: list[tuple[str, int, int, str]] = []
    read_truth_bp: dict[str, int] = defaultdict(int)
    target = coverage * len(source)
    total = 0
    read_index = 0
    while total < target:
        start = rng.randrange(0, len(source) - 1000 + 1)
        seq = source[start:start + 1000]
        read_id = f"{seed_label}_{condition}_r{read_index:06d}"
        reads.append((read_id, seq, start))
        for family, rows in source_intervals.items():
            for left, right, _ in rows:
                overlap_left, overlap_right = max(start, left), min(start + len(seq), right)
                if overlap_left < overlap_right:
                    l, r = overlap_left - start, overlap_right - start
                    truth_intervals.append((read_id, l, r, family))
                    read_truth_bp[family] += r - l
        total += len(seq)
        read_index += 1
    reads_out = [(read_id, seq) for read_id, seq, _ in reads]
    reads_digest = _digest("\n".join(f"{rid}\t{seq}" for rid, seq in reads_out))
    return {
        "catalogue": catalogue,
        "reads": reads_out,
        "truth_intervals": truth_intervals,
        "source_truth_intervals": source_intervals,
        "source_truth_bp": source_truth_bp,
        "read_truth_bp": dict(read_truth_bp),
        "read_bases": total,
        "genome_size": len(source),
        "condition": condition,
        "coverage": coverage,
        "source_hash": _digest(source),
        "reads_hash": reads_digest,
    }
