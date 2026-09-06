"""Independent factorial genomes written with bounded memory and random access."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import random

from benchmarks.challenge.schema import digest_file, write_table


@dataclass(frozen=True)
class ArraySpec:
    family_id: str
    period: int
    copies: int
    gc: float = 0.5
    unit_substitution_rate: float = 0.0

    def validate(self) -> None:
        if (not self.family_id or not all(c.isalnum() or c == '_' for c in self.family_id)
                or not isinstance(self.period, int) or not isinstance(self.copies, int)
                or not 2 <= self.period <= 20_000 or self.copies < 2
                or not math.isfinite(self.gc) or not 0 < self.gc < 1
                or not math.isfinite(self.unit_substitution_rate) or not 0 <= self.unit_substitution_rate < 1):
            raise ValueError('Invalid array family, period/copies, GC or biological divergence')


class FixedFastaWriter:
    """Single-contig fixed-width FASTA; retain less than one output line."""
    def __init__(self, handle, width: int = 80):
        self.handle, self.width, self.pending = handle, width, ''
        self.bases = 0
        self.header = b'>chr_sim\n'
        handle.write(self.header)

    def add(self, sequence: str) -> None:
        self.bases += len(sequence)
        sequence = self.pending+sequence
        whole = len(sequence)//self.width*self.width
        for start in range(0, whole, self.width):
            self.handle.write(sequence[start:start+self.width].encode()+b'\n')
        self.pending = sequence[whole:]

    def finish(self) -> dict:
        if self.pending:
            self.handle.write(self.pending.encode()+b'\n')
        return dict(contig='chr_sim', length=self.bases, offset=len(self.header),
                    line_bases=self.width, line_bytes=self.width+1)


class FixedFastaReader:
    """Random access to a verified generated fixed-width single-contig file."""
    def __init__(self, handle, index: dict):
        self.handle, self.index = handle, index

    def get(self, start: int, end: int) -> str:
        meta = self.index
        if not 0 <= start <= end <= meta['length']:
            raise ValueError('Random-access bounds exceed generated genome')
        if start == end:
            return ''
        width = meta['line_bases']
        offset = meta['offset'] + start//width*meta['line_bytes'] + start % width
        self.handle.seek(offset)
        size = end-start + (end-1)//width-start//width
        sequence = self.handle.read(size).replace(b'\n', b'').decode('ascii')
        if len(sequence) != end-start or set(sequence)-set('ACGT'):
            raise ValueError('Generated FASTA/index mismatch or invalid sequence')
        return sequence


def generate(specs: list[ArraySpec], genome_bp: int, seed: int, outdir: Path,
             background_gc: float = .45) -> dict:
    """Uniform IID backgrounds and independently diverged copies, not plant chemistry."""
    if not specs or len({s.family_id for s in specs}) != len(specs):
        raise ValueError('Require distinct nonempty array families')
    for spec in specs:
        spec.validate()
    repeat_bp = sum(s.period*s.copies for s in specs)
    if (not isinstance(genome_bp, int) or not 0 < background_gc < 1
            or not repeat_bp+len(specs)+1 <= genome_bp <= 1_000_000_000):
        raise ValueError('Genome must fit repeats plus flanks, with an explicit 1 Gb generator limit')
    outdir.mkdir(parents=True, exist_ok=False)
    # Separate streams retain identical founders/background at different divergence.
    motif_rng, background_rng, variant_rng = [random.Random(seed+n) for n in (0, 104729, 209759)]
    monomers = {s.family_id: ''.join(motif_rng.choices('ACGT',
                 weights=[(1-s.gc)/2, s.gc/2, s.gc/2, (1-s.gc)/2], k=s.period)) for s in specs}
    flank, remainder = divmod(genome_bp-repeat_bp, len(specs)+1)
    truth = []
    with (outdir/'genome.fa').open('wb') as handle:
        writer = FixedFastaWriter(handle)
        for i in range(len(specs)+1):
            remaining = flank+(i < remainder)
            while remaining:
                amount = min(remaining, 65_536)
                writer.add(''.join(background_rng.choices('ACGT',
                    weights=[(1-background_gc)/2, background_gc/2, background_gc/2, (1-background_gc)/2], k=amount)))
                remaining -= amount
            if i == len(specs):
                break
            spec = specs[i]
            monomer = monomers[spec.family_id]
            start, substitutions = writer.bases, 0
            for _ in range(spec.copies):
                unit = []
                for base in monomer:
                    change = variant_rng.random() < spec.unit_substitution_rate
                    unit.append(variant_rng.choice('ACGT'.replace(base, '')) if change else base)
                    substitutions += change
                writer.add(''.join(unit))
            truth.append(dict(chrom='chr_sim', family_id=spec.family_id, start=start, end=writer.bases,
                              period=spec.period, copies=spec.copies, repeat_bp=writer.bases-start,
                              requested_gc=spec.gc, founder_gc=(monomer.count('G')+monomer.count('C'))/len(monomer),
                              unit_substitution_rate=spec.unit_substitution_rate, observed_unit_substitutions=substitutions,
                              founder_sha256=hashlib.sha256(monomer.encode()).hexdigest()))
        index = writer.finish()
    assert index['length'] == genome_bp
    (outdir/'genome_index.json').write_text(json.dumps(index, indent=2)+'\n')
    (outdir/'catalogue.fa').write_text(''.join(f'>{f}\n{s}\n' for f, s in monomers.items()))
    write_table(outdir/'truth_copy_number.tsv', truth, list(truth[0]))
    manifest = dict(schema_version=1, generator='streamed_factorial_genome_v1', seed=seed,
                    genome_bp=genome_bp, background_gc=background_gc, array_specs=[asdict(s) for s in specs],
                    warning='IID_background_haploid_independent_unit_substitutions_not_real_plant_genome',
                    script_sha256=digest_file(Path(__file__)),
                    files={p.name: digest_file(p) for p in outdir.iterdir() if p.is_file()})
    (outdir/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest
