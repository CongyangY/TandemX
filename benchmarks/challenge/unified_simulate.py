"""Formal comparator fixtures; historical simulator source remains frozen."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import random

from .adapters import read_fasta
from .schema import digest_file, read_table, write_table
from .simulate import Scenario, generate_dataset, random_dna, reverse_complement


def shared_fragment_read(length: int, monomer: str, rng: random.Random) -> str:
    if len(monomer) <= 100:
        raise ValueError('shared_fragment requires period >100 bp')
    parts = [random_dna(rng, 200)]
    size = 200
    while size < length:
        gap = random_dna(rng, rng.randint(30, 80))
        parts.extend([monomer[:100], gap])
        size += 100 + len(gap)
    return ''.join(parts)[:length]


def generate_unified_dataset(config: Scenario, seed: int, outdir: Path) -> dict:
    if config.negative_kind != 'shared_fragment':
        return generate_dataset(config, seed, outdir)
    if config.period <= 100:
        raise ValueError('shared_fragment requires period >100 bp')
    manifest = generate_dataset(replace(config, negative_kind='random'), seed, outdir)
    monomer = next(iter(read_fasta(outdir / 'truth_monomers.fa').values()))
    reads = read_fasta(outdir / 'reads.fa')
    labels = read_table(outdir / 'truth_reads.tsv')
    for index, row in enumerate(labels):
        if int(row['truth_positive']) == 0:
            # Independent declared stream; positive arrays remain unchanged.
            seq = shared_fragment_read(config.read_length, monomer, random.Random(seed * 1000 + index))
            reads[row['read_id']] = reverse_complement(seq) if row['strand'] == '-' else seq
            row['negative_kind'] = 'shared_fragment'
    with (outdir / 'reads.fa').open('w') as handle:
        for name, seq in reads.items():
            handle.write(f'>{name}\n{seq}\n')
    write_table(outdir / 'truth_reads.tsv', labels, list(labels[0]))
    manifest.update(scenario=asdict(config), simulator='independent_planted_arrays_v1_shared_fragment_negative_overlay',
                    negative_rng='Random(seed*1000+zero_based_read_index)')
    manifest['files'] = {p.name: {'sha256': digest_file(p), 'size_bytes': p.stat().st_size}
                         for p in sorted(outdir.iterdir()) if p.is_file() and p.name != 'manifest.json'}
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest
