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
    unit_substitution_rate: float = 0.0
    array_fragments: int = 1
    fragment_gap_bp: int = 0

    def validate(self) -> None:
        if not self.periods or len(self.periods) != len(self.copies):
            raise ValueError("Need one copy count per period")
        if any(p < 2 for p in self.periods) or any(c < 2 for c in self.copies) or self.flank_bp < 1:
            raise ValueError("Periods/copies must be >=2; flanks must be positive")
        if (
            not math.isfinite(self.unit_substitution_rate)
            or not 0 <= self.unit_substitution_rate < 1
            or self.array_fragments < 1
            or self.fragment_gap_bp < 0
        ):
            raise ValueError("Invalid unit-divergence or array-fragmentation settings")
        gaps = len(self.periods) * (self.array_fragments - 1) * self.fragment_gap_bp
        size = ((len(self.periods) + 1) * self.flank_bp
                + sum(p*c for p, c in zip(self.periods, self.copies)) + gaps)
        if size > self.max_genome_bp:
            raise ValueError("Genome exceeds the explicit bounded-simulation memory limit")


def _copy_sequence(founder: str, rate: float, rng: random.Random) -> str:
    if rate == 0:
        return founder
    bases = list(founder)
    for index, base in enumerate(bases):
        if rng.random() < rate:
            bases[index] = rng.choice(DNA.replace(base, ""))
    return "".join(bases)


def _split_copies(copies: int, fragments: int) -> list[int]:
    if copies == 0:
        return []
    used = min(copies, fragments)
    quotient, remainder = divmod(copies, used)
    return [quotient + (index < remainder) for index in range(used)]


def challenge_scenarios(config: dict) -> list[tuple[float, int]]:
    unit_rates = config.get("unit_substitution_rates", [0.0])
    fragment_counts = config.get("array_fragment_counts", [1])
    if (
        not unit_rates or len(set(unit_rates)) != len(unit_rates)
        or any(not math.isfinite(rate) or not 0 <= rate < 1 for rate in unit_rates)
        or not fragment_counts or len(set(fragment_counts)) != len(fragment_counts)
        or any(not isinstance(count, int) or isinstance(count, bool) or count < 1
               for count in fragment_counts)
    ):
        raise ValueError("Invalid unit-divergence or array-fragmentation matrix")
    return [(float(rate), count) for rate in unit_rates for count in fragment_counts]


def scenario_directory(base: Path, scenario_index: int, scenarios: list[tuple[float, int]],
                       config: dict) -> Path:
    explicit = "unit_substitution_rates" in config or "array_fragment_counts" in config
    if len(scenarios) == 1 and scenarios[0] == (0.0, 1) and not explicit:
        return base
    return base / f"scenario_{scenario_index:03d}"


def build_genome(spec: GenomeSpec, fraction: float = 1.0) -> tuple[str, dict[str, str], list[dict]]:
    spec.validate()
    if not math.isfinite(fraction) or not 0 <= fraction <= 2:
        raise ValueError("Assembly copy fraction must be in [0,2]")
    retained_copies = [math.floor(c*fraction+.5) for c in spec.copies]
    fragment_counts = [min(copies, spec.array_fragments) if copies else 0 for copies in retained_copies]
    size = ((len(spec.periods)+1)*spec.flank_bp
            + sum(p*c for p,c in zip(spec.periods,retained_copies))
            + sum(max(0, count-1)*spec.fragment_gap_bp for count in fragment_counts))
    if size > spec.max_genome_bp:
        raise ValueError("Assembly exceeds the explicit bounded-simulation memory limit")
    rng = random.Random(spec.seed)
    monomers = {f"f{i+1}": "".join(rng.choices(DNA, k=p)) for i, p in enumerate(spec.periods)}
    # All monomers/flanks are generated before truncating arrays: assembly versions
    # have identical backgrounds and differ only in planted copy counts.
    flanks = ["".join(rng.choices(DNA, k=spec.flank_bp)) for _ in range(len(monomers) + 1)]
    fragment_gaps = {
        family: ["".join(rng.choices(DNA, k=spec.fragment_gap_bp))
                 for _ in range(spec.array_fragments - 1)]
        for family in monomers
    }
    parts, truth, position = [flanks[0]], [], len(flanks[0])
    for i, (family, sequence) in enumerate(monomers.items()):
        copies = retained_copies[i]
        copy_rng = random.Random((spec.seed + 1) * 1_000_003 + (i + 1) * 104_729)
        variants = [_copy_sequence(sequence, spec.unit_substitution_rate, copy_rng)
                    for _ in range(copies)]
        counts = _split_copies(copies, spec.array_fragments)
        if not counts:
            truth.append(dict(chrom="chr_sim", family_id=family, array_index=1,
                              start=position, end=position, period=len(sequence),
                              copies=0, repeat_bp=0))
        offset = 0
        for fragment_index, fragment_copies in enumerate(counts, 1):
            fragment = "".join(variants[offset:offset+fragment_copies])
            start, end = position, position + len(fragment)
            truth.append(dict(chrom="chr_sim", family_id=family, array_index=fragment_index,
                              start=start, end=end, period=len(sequence),
                              copies=fragment_copies, repeat_bp=end-start))
            parts.append(fragment)
            position = end
            offset += fragment_copies
            if fragment_index < len(counts):
                gap = fragment_gaps[family][fragment_index-1]
                parts.append(gap)
                position += len(gap)
        parts.append(flanks[i+1])
        position += len(flanks[i+1])
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
