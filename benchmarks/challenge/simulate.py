"""Streaming planted-array simulations with observed-coordinate truth.

These are controlled stress tests, not an empirical model of PacBio chemistry.
Reads have neutral IDs; truth is written separately and is never a tool input.
"""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path

from .schema import ArrayRecord, digest_file, write_table

DNA = "ACGT"


@dataclass(frozen=True)
class Scenario:
    name: str
    period: int = 171
    read_count: int = 100
    read_length: int = 5000
    families: int = 3
    copies: int = 12
    arrays_per_read: int = 1
    positive_fraction: float = 0.7
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    unit_divergence: float = 0.0
    related_families: bool = False
    negative_kind: str = "random"

    def validate(self) -> None:
        if not self.name or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789_" for c in self.name):
            raise ValueError("Scenario name must be lowercase alphanumeric/underscore")
        for field in ("period", "read_count", "read_length", "families", "copies", "arrays_per_read"):
            if getattr(self, field) < 1:
                raise ValueError(f"{field} must be positive")
        if self.copies < 2 or self.families < self.arrays_per_read:
            raise ValueError("Need >=2 copies and enough distinct families for each read")
        if self.period * self.copies * self.arrays_per_read + 200 * (self.arrays_per_read + 1) > self.read_length:
            raise ValueError("Read length cannot accommodate arrays and random flanks")
        for value in (self.positive_fraction, self.substitution_rate, self.insertion_rate,
                      self.deletion_rate, self.unit_divergence):
            if not 0 <= value <= 1:
                raise ValueError("Fractions must be in [0,1]")
        if self.deletion_rate == 1:
            raise ValueError("Deletion rate must be <1")
        if self.negative_kind not in {"random", "at_rich", "dispersed", "low_complexity"}:
            raise ValueError(f"Unknown negative kind: {self.negative_kind}")


def random_dna(rng: random.Random, size: int, *, at_rich: bool = False) -> str:
    alphabet = "AAAAAAAATTTTTTTTCG" if at_rich else DNA
    return "".join(rng.choices(alphabet, k=size))


def mutate(sequence: str, rng: random.Random, sub: float = 0.0,
           ins: float = 0.0, deletion: float = 0.0) -> str:
    """Independent per-source-base deletion, substitution and single insertion."""
    result: list[str] = []
    for base in sequence:
        if rng.random() >= deletion:
            result.append(rng.choice(DNA.replace(base, "")) if rng.random() < sub else base)
        if rng.random() < ins:
            result.append(rng.choice(DNA))
    return "".join(result)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def negative_read(config: Scenario, monomer: str, rng: random.Random) -> str:
    if config.negative_kind == "at_rich":
        return random_dna(rng, config.read_length, at_rich=True)
    if config.negative_kind == "low_complexity":
        # Genuine STRs outside the declared satellite period range, not satellites.
        return ("A" * (config.read_length // 2) + "AT" * config.read_length)[:config.read_length]
    if config.negative_kind == "dispersed":
        # Separated single copies with independently sampled, unequal long gaps.
        parts = [random_dna(rng, 300)]
        while sum(map(len, parts)) + 4 * len(monomer) < config.read_length:
            parts.extend([monomer, random_dna(rng, rng.randrange(2 * len(monomer), 3 * len(monomer)))])
        return ("".join(parts) + random_dna(rng, config.read_length))[:config.read_length]
    return random_dna(rng, config.read_length)


def generate_dataset(config: Scenario, seed: int, outdir: Path) -> dict:
    config.validate()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite dataset: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    founder = random_dna(rng, config.period)
    monomers = {
        f"truth_f{i + 1}": (founder if i == 0 else mutate(founder, rng, sub=0.12))
        if config.related_families else random_dna(rng, config.period)
        for i in range(config.families)
    }
    labels = [i < round(config.read_count * config.positive_fraction) for i in range(config.read_count)]
    rng.shuffle(labels)
    arrays: list[ArrayRecord] = []
    reads: list[dict] = []
    with (outdir / "reads.fa").open("w", encoding="utf-8") as fasta:
        for index, positive in enumerate(labels):
            read_id = f"r{index + 1:07d}"
            local: list[ArrayRecord] = []
            if positive:
                selected = rng.sample(list(monomers), config.arrays_per_read)
                flanks_bp = config.read_length - config.period * config.copies * config.arrays_per_read
                flank = flanks_bp // (config.arrays_per_read + 1)
                sequence = random_dna(rng, flank)
                for family_id in selected:
                    monomer = monomers[family_id]
                    phase = rng.randrange(len(monomer))
                    rotated = monomer[phase:] + monomer[:phase]
                    array = "".join(mutate(rotated, rng, sub=config.unit_divergence) for _ in range(config.copies))
                    array = mutate(array, rng, config.substitution_rate, config.insertion_rate, config.deletion_rate)
                    start = len(sequence)
                    sequence += array
                    local.append(ArrayRecord(read_id, start, len(sequence), config.period, monomer, family_id))
                    sequence += random_dna(rng, flank)
                sequence += random_dna(rng, flanks_bp - flank * (config.arrays_per_read + 1))
            else:
                sequence = negative_read(config, next(iter(monomers.values())), rng)
            strand = "-" if rng.random() < 0.5 else "+"
            if strand == "-":
                sequence = reverse_complement(sequence)
                local = [ArrayRecord(r.read_id, len(sequence) - r.end, len(sequence) - r.start,
                                     r.period, r.sequence, r.family_id) for r in local]
            arrays.extend(sorted(local, key=lambda r: r.start))
            fasta.write(f">{read_id}\n{sequence}\n")
            reads.append({"read_id": read_id, "length_bp": len(sequence), "truth_positive": int(positive),
                          "strand": strand, "negative_kind": "NA" if positive else config.negative_kind})
    write_table(outdir / "truth_arrays.tsv", (asdict(r) for r in arrays), list(ArrayRecord.__dataclass_fields__))
    write_table(outdir / "truth_reads.tsv", reads, list(reads[0]))
    with (outdir / "truth_monomers.fa").open("w", encoding="utf-8") as fasta:
        for family_id, monomer in monomers.items():
            fasta.write(f">{family_id}\n{monomer}\n")
    manifest = {"schema_version": 1, "scenario": asdict(config), "seed": seed,
                "read_count": len(reads), "total_bases": sum(r["length_bp"] for r in reads),
                "truth_arrays": len(arrays), "simulator": "independent_planted_arrays_v1",
                "warning": "controlled_stress_test_not_empirical_hifi_error_model",
                "files": {p.name: {"sha256": digest_file(p), "size_bytes": p.stat().st_size}
                          for p in sorted(outdir.iterdir()) if p.is_file()}}
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
