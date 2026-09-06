"""Bounded genomes with planted copy truth and uniformly sampled circular reads.

This is a sampling experiment, not a PacBio chemistry or diploid genome model.
Its known catalogue tests quantification conditional on correct monomers.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import random

from benchmarks.challenge.schema import digest_file, write_table

DNA = "ACGT"


@dataclass(frozen=True)
class GenomeSpec:
    seed: int
    periods: tuple[int, ...] = (61, 171, 421)
    copies: tuple[int, ...] = (20, 80, 200)
    flank_bp: int = 25000
    max_genome_bp: int = 2000000

    def validate(self) -> None:
        if not self.periods or len(self.periods) != len(self.copies):
            raise ValueError("Need one copy count per period")
        if any(p < 2 for p in self.periods) or any(c < 2 for c in self.copies) or self.flank_bp < 1:
            raise ValueError("Periods/copies must be >=2; flanks must be positive")
        size = (len(self.periods) + 1) * self.flank_bp + sum(p*c for p, c in zip(self.periods, self.copies))
        if size > self.max_genome_bp:
            raise ValueError("Genome exceeds the explicit bounded-simulation memory limit")


def build_genome(spec: GenomeSpec, fraction: float = 1.0) -> tuple[str, dict[str, str], list[dict]]:
    spec.validate()
    if not math.isfinite(fraction) or not 0 <= fraction <= 2:
        raise ValueError("Assembly copy fraction must be in [0,2]")
    size = (len(spec.periods)+1)*spec.flank_bp + sum(p*math.floor(c*fraction+.5) for p,c in zip(spec.periods,spec.copies))
    if size > spec.max_genome_bp:
        raise ValueError("Assembly exceeds the explicit bounded-simulation memory limit")
    rng = random.Random(spec.seed)
    monomers = {f"f{i+1}": "".join(rng.choices(DNA, k=p)) for i, p in enumerate(spec.periods)}
    # All monomers/flanks are generated before truncating arrays: assembly versions
    # have identical backgrounds and differ only in planted copy counts.
    flanks = ["".join(rng.choices(DNA, k=spec.flank_bp)) for _ in range(len(monomers) + 1)]
    parts, truth, position = [flanks[0]], [], len(flanks[0])
    for i, (family, sequence) in enumerate(monomers.items()):
        copies = math.floor(spec.copies[i] * fraction + .5)
        end = position + copies * len(sequence)
        truth.append(dict(chrom="chr_sim", family_id=family, start=position, end=end,
                          period=len(sequence), copies=copies, repeat_bp=end-position))
        parts.extend([sequence * copies, flanks[i+1]])
        position = end + len(flanks[i+1])
    return "".join(parts), monomers, truth


def write_genome(spec: GenomeSpec, outdir: Path, fractions: tuple[float, ...]) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    genome, monomers, truth = build_genome(spec)
    (outdir / "genome.fa").write_text(f">chr_sim\n{genome}\n")
    (outdir / "catalogue.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in monomers.items()))
    write_table(outdir / "truth_copy_number.tsv", truth, list(truth[0]))
    variants = []
    for index, fraction in enumerate(fractions):
        assembly, _, retained = build_genome(spec, fraction)
        name = f"assembly_{index}"
        (outdir / f"{name}.fa").write_text(f">chr_sim\n{assembly}\n")
        write_table(outdir / f"{name}.truth.tsv", retained, list(retained[0]))
        variants.append(dict(name=name, fraction=fraction, genome_bp=len(assembly)))
    manifest = {"spec": asdict(spec), "genome_bp": len(genome), "variants": variants,
                "catalogue_role": "known planted catalogue; conditional abundance/localization, not de novo recovery",
                "files": {p.name: digest_file(p) for p in outdir.iterdir()}}
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def sample_reads(genome: str, truth: list[dict], outdir: Path, *, seed: int,
                 coverage: float, read_length: int, substitution_rate: float) -> dict:
    if not genome or set(genome) - set(DNA):
        raise ValueError("Genome must contain ACGT")
    if not math.isfinite(coverage) or coverage <= 0 or not 1 <= read_length <= len(genome):
        raise ValueError("Positive finite coverage and read length <= genome length are required")
    if not math.isfinite(substitution_rate) or not 0 <= substitution_rate < 1:
        raise ValueError("Substitution rate must be in [0,1)")
    outdir.mkdir(parents=True, exist_ok=False)
    # Independent RNG streams keep starts/strands identical across error settings.
    starts_rng, errors_rng = random.Random(seed), random.Random(seed + 104729)
    count = max(1, math.floor(coverage * len(genome) / read_length + .5))
    sampled_bp = {row["family_id"]: 0 for row in truth}
    changed = 0
    with (outdir / "reads.fa").open("w") as fasta, (outdir / "sampling.tsv").open("w") as audit:
        audit.write("read_id\tgenome_start\tsource_length\tstrand\tsubstitutions\n")
        for i in range(count):
            start = starts_rng.randrange(len(genome))
            reverse = starts_rng.random() < .5
            sequence = genome[start:start+read_length]
            if len(sequence) < read_length:
                sequence += genome[:read_length-len(sequence)]
            spans = [(start, min(start+read_length, len(genome)))]
            if start + read_length > len(genome):
                spans.append((0, start+read_length-len(genome)))
            for row in truth:
                sampled_bp[row["family_id"]] += sum(max(0, min(b,row["end"])-max(a,row["start"])) for a,b in spans)
            bases, substitutions = list(sequence), 0
            for j, base in enumerate(bases):
                if errors_rng.random() < substitution_rate:
                    bases[j] = errors_rng.choice(DNA.replace(base, ""))
                    substitutions += 1
            sequence = "".join(bases)
            if reverse:
                sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            name = f"r{i+1:08d}"
            fasta.write(f">{name}\n{sequence}\n")
            audit.write(f"{name}\t{start}\t{read_length}\t{'-' if reverse else '+'}\t{substitutions}\n")
            changed += substitutions
    manifest = dict(seed=seed, requested_coverage=coverage, actual_base_coverage=count*read_length/len(genome),
                    read_count=count, read_length=read_length, genome_bp=len(genome),
                    substitution_rate=substitution_rate, observed_substitutions=changed,
                    sampled_repeat_bp=sampled_bp, source_genome_sha256=hashlib.sha256(genome.encode()).hexdigest(),
                    warning="uniform circular sampling; iid substitutions only; exact founder copies; no ploidy or empirical HiFi chemistry",
                    files={p.name:digest_file(p) for p in outdir.iterdir()})
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
