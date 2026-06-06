"""Reproducible toy dataset generator for TandemX."""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


DNA_ALPHABET = "ACGT"


@dataclass(frozen=True)
class ToyFamily:
    family_id: str
    monomer_id: str
    monomer_length: int
    read_copies: int
    assembly_copies: int
    monomer_sequence: str


@dataclass(frozen=True)
class ToySimulationConfig:
    outdir: Path
    seed: int
    num_reads: int
    read_length: int
    background_length: int
    monomer_lengths: tuple[int, ...]
    copies: tuple[int, ...]
    error_rate: float


def parse_int_list(value: str, option_name: str) -> tuple[int, ...]:
    """Parse a comma-separated list of positive integers."""
    try:
        parsed = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    except ValueError as exc:
        raise ValueError(f"{option_name} must be a comma-separated list of integers") from exc
    if not parsed:
        raise ValueError(f"{option_name} must contain at least one integer")
    if any(item <= 0 for item in parsed):
        raise ValueError(f"{option_name} values must be positive integers")
    return parsed


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(table)[::-1].upper()


def generate_toy_dataset(config: ToySimulationConfig) -> None:
    """Generate a deterministic toy dataset with truth files."""
    validate_config(config)
    rng = random.Random(config.seed)
    config.outdir.mkdir(parents=True, exist_ok=True)

    families = build_families(config, rng)
    background = random_sequence(config.background_length, rng)
    source_sequence = build_source_sequence(background, families)
    assembly_sequence, assembly_arrays = build_assembly_sequence(background, families)
    reads = sample_reads(source_sequence, config, rng)

    write_fasta(config.outdir / "reads.fa", reads)
    write_fasta(
        config.outdir / "assembly.fa",
        [("toy_chr1", assembly_sequence)],
    )
    write_fasta(
        config.outdir / "truth_monomers.fa",
        [
            (
                (
                    f"family_id={family.family_id};monomer_id={family.monomer_id};"
                    f"length_bp={family.monomer_length};source=simulated"
                ),
                family.monomer_sequence,
            )
            for family in families
        ],
    )
    write_truth_arrays(config.outdir / "truth_arrays.bed", assembly_arrays)
    write_truth_copy_number(config.outdir / "truth_copy_number.tsv", families)
    write_simulation_config(config.outdir / "simulation_config.yaml", config, families)


def validate_config(config: ToySimulationConfig) -> None:
    if config.num_reads <= 0:
        raise ValueError("--num-reads must be positive")
    if config.read_length <= 0:
        raise ValueError("--read-length must be positive")
    if config.background_length <= 0:
        raise ValueError("--background-length must be positive")
    if len(config.monomer_lengths) != len(config.copies):
        raise ValueError("--monomer-lengths and --copies must have the same number of values")
    if not 0 <= config.error_rate <= 1:
        raise ValueError("--error-rate must be between 0 and 1")


def build_families(config: ToySimulationConfig, rng: random.Random) -> list[ToyFamily]:
    families = []
    for index, (length, copies) in enumerate(zip(config.monomer_lengths, config.copies), start=1):
        assembly_copies = copies
        if index == 1:
            assembly_copies = max(1, copies // 3)
        families.append(
            ToyFamily(
                family_id=f"TXF{index:06d}",
                monomer_id=f"TXM{index:06d}",
                monomer_length=length,
                read_copies=copies,
                assembly_copies=assembly_copies,
                monomer_sequence=random_sequence(length, rng),
            )
        )
    return families


def random_sequence(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA_ALPHABET) for _ in range(length))


def build_source_sequence(background: str, families: Sequence[ToyFamily]) -> str:
    pieces = [background[: len(background) // 2]]
    for family in families:
        pieces.append(family.monomer_sequence * family.read_copies)
        pieces.append(background[:100])
    pieces.append(background[len(background) // 2 :])
    return "".join(pieces)


def build_assembly_sequence(
    background: str,
    families: Sequence[ToyFamily],
) -> tuple[str, list[tuple[str, int, int, str, int, str, int, int, str]]]:
    pieces = [background[: len(background) // 2]]
    arrays = []
    cursor = len(pieces[0])
    for family in families:
        array_sequence = family.monomer_sequence * family.assembly_copies
        start = cursor
        end = start + len(array_sequence)
        arrays.append(
            (
                "toy_chr1",
                start,
                end,
                family.family_id,
                1000,
                "+",
                family.monomer_length,
                family.assembly_copies,
                "compressed" if family.assembly_copies < family.read_copies else "truth_like",
            )
        )
        pieces.append(array_sequence)
        cursor = end
        spacer = background[:100]
        pieces.append(spacer)
        cursor += len(spacer)
    tail = background[len(background) // 2 :]
    pieces.append(tail)
    return "".join(pieces), arrays


def sample_reads(
    source_sequence: str,
    config: ToySimulationConfig,
    rng: random.Random,
) -> list[tuple[str, str]]:
    reads = []
    max_start = max(0, len(source_sequence) - config.read_length)
    for index in range(1, config.num_reads + 1):
        if index <= 2:
            start = min(config.background_length // 2, max_start)
        else:
            start = rng.randint(0, max_start) if max_start else 0
        sequence = source_sequence[start : start + config.read_length]
        if len(sequence) < config.read_length:
            sequence += source_sequence[: config.read_length - len(sequence)]
        strand = "-" if index % 2 == 0 else "+"
        if strand == "-":
            sequence = reverse_complement(sequence)
        sequence = introduce_errors(sequence, config.error_rate, rng)
        header = f"toy_read_{index:04d};source_start={start};strand={strand};error_rate={config.error_rate}"
        reads.append((header, sequence))
    return reads


def introduce_errors(sequence: str, error_rate: float, rng: random.Random) -> str:
    bases = []
    for base in sequence:
        if rng.random() < error_rate:
            alternatives = [candidate for candidate in DNA_ALPHABET if candidate != base]
            bases.append(rng.choice(alternatives))
        else:
            bases.append(base)
    return "".join(bases)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_fasta(path: Path, records: Sequence[tuple[str, str]]) -> None:
    lines = []
    for header, sequence in records:
        lines.append(f">{header}")
        lines.append(wrap_sequence(sequence))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_truth_arrays(
    path: Path,
    arrays: Sequence[tuple[str, int, int, str, int, str, int, int, str]],
) -> None:
    lines = []
    for row in arrays:
        lines.append("\t".join(str(value) for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_truth_copy_number(path: Path, families: Sequence[ToyFamily]) -> None:
    lines = [
        (
            "family_id\tmonomer_id\tmonomer_length_bp\tread_copies\tassembly_copies\t"
            "read_repeat_bp\tassembly_repeat_bp\tassembly_status"
        )
    ]
    for family in families:
        status = (
            "simulated_under_assembly"
            if family.assembly_copies < family.read_copies
            else "truth_like"
        )
        lines.append(
            "\t".join(
                str(value)
                for value in (
                    family.family_id,
                    family.monomer_id,
                    family.monomer_length,
                    family.read_copies,
                    family.assembly_copies,
                    family.monomer_length * family.read_copies,
                    family.monomer_length * family.assembly_copies,
                    status,
                )
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_simulation_config(
    path: Path,
    config: ToySimulationConfig,
    families: Sequence[ToyFamily],
) -> None:
    lines = [
        'command: "tandemx simulate toy"',
        f"seed: {config.seed}",
        f"num_reads: {config.num_reads}",
        f"read_length: {config.read_length}",
        f"background_length: {config.background_length}",
        f"error_rate: {config.error_rate}",
        "monomer_lengths:",
    ]
    for length in config.monomer_lengths:
        lines.append(f"  - {length}")
    lines.append("copies:")
    for copies in config.copies:
        lines.append(f"  - {copies}")
    lines.append("families:")
    for family in families:
        lines.extend(
            [
                f"  - family_id: \"{family.family_id}\"",
                f"    monomer_id: \"{family.monomer_id}\"",
                f"    monomer_length_bp: {family.monomer_length}",
                f"    read_copies: {family.read_copies}",
                f"    assembly_copies: {family.assembly_copies}",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
