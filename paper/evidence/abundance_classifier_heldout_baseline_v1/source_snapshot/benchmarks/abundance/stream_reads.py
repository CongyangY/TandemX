"""Uniform source sampling with empirical lengths and observed indel coordinates."""
from __future__ import annotations

from bisect import bisect_right
from collections import Counter
import csv
import json
import math
from pathlib import Path
import random

from benchmarks.abundance.stream_genome import FixedFastaReader
from benchmarks.challenge.schema import digest_file, read_table
from benchmarks.challenge.simulate import reverse_complement


def alter_segment(sequence: str, rng: random.Random, sub: float, ins: float, deletion: float) -> tuple[str, tuple[int, int, int]]:
    """Per-source-base events, including insertions after deleted bases."""
    if sub == ins == deletion == 0:
        return sequence, (0, 0, 0)
    output, substitutions, insertions, deletions = [], 0, 0, 0
    for base in sequence:
        if rng.random() >= deletion:
            changed = rng.random() < sub
            output.append(rng.choice('ACGT'.replace(base, '')) if changed else base)
            substitutions += changed
        else:
            deletions += 1
        if rng.random() < ins:
            output.append(rng.choice('ACGT'))
            insertions += 1
    return ''.join(output), (substitutions, insertions, deletions)


def length_distribution(histogram: Path, maximum: int) -> tuple[list[int], list[int]]:
    rows = read_table(histogram, {'length_bp', 'read_count'})
    counts = {}
    for row in rows:
        length, count = int(row['length_bp']), int(row['read_count'])
        if not 1 <= length <= maximum or count <= 0 or length in counts:
            raise ValueError('Length histogram has invalid, oversized or duplicate bins')
        counts[length] = count
    lengths, cumulative, total = sorted(counts), [], 0
    for length in lengths:
        total += counts[length]
        cumulative.append(total)
    if not total:
        raise ValueError('Empty length distribution')
    return lengths, cumulative


def sample(genome_dir: Path, outdir: Path, histogram: Path, *, seed: int, coverage: float,
           substitution_rate: float, insertion_rate: float = 0, deletion_rate: float = 0) -> dict:
    if (not math.isfinite(coverage) or not 0 < coverage <= 100
            or any(not math.isfinite(x) or not 0 <= x < 1 for x in (substitution_rate, insertion_rate, deletion_rate))):
        raise ValueError('Invalid finite coverage or error probabilities')
    genome = json.loads((genome_dir/'manifest.json').read_text())
    for name in ('genome.fa', 'genome_index.json', 'truth_copy_number.tsv', 'catalogue.fa'):
        if digest_file(genome_dir/name) != genome['files'][name]:
            raise ValueError('Generated source/hash mismatch')
    index = json.loads((genome_dir/'genome_index.json').read_text())
    if index['length'] != genome['genome_bp']:
        raise ValueError('Genome length differs from index')
    lengths, cumulative = length_distribution(histogram, min(index['length'], 200_000))
    truth = read_table(genome_dir/'truth_copy_number.tsv')
    for row in truth:
        for field in ('start', 'end', 'period'):
            row[field] = int(row[field])
    starts = [r['start'] for r in truth]
    if any(a['end'] > b['start'] for a, b in zip(truth, truth[1:])):
        raise ValueError('Source arrays overlap or are unsorted')
    outdir.mkdir(parents=True, exist_ok=False)
    design_rng, errors_rng = random.Random(seed), random.Random(seed+314159)
    target = math.ceil(coverage*genome['genome_bp'])
    source_bases = observed_bases = read_count = 0
    sampled_repeat_bp = Counter()
    read_array_count = 0
    event_totals = [0, 0, 0]

    def alter(sequence: str) -> str:
        observed, events = alter_segment(sequence, errors_rng, substitution_rate, insertion_rate, deletion_rate)
        for i, number in enumerate(events):
            event_totals[i] += number
        return observed
    with (genome_dir/'genome.fa').open('rb') as source, (outdir/'reads.fa').open('w') as output, \
            (outdir/'sampling.tsv').open('w', newline='') as audit, (outdir/'truth_read_segments.tsv').open('w', newline='') as truth_out:
        reader = FixedFastaReader(source, index)
        audit_writer = csv.writer(audit, delimiter='\t')
        audit_writer.writerow(['read_id','genome_start','source_length','observed_length','strand'])
        truth_writer = csv.writer(truth_out, delimiter='\t')
        truth_writer.writerow(['read_id','start','end','period','family_id','sampled_source_repeat_bp','at_least_two_source_units'])
        while source_bases < target:
            length = lengths[bisect_right(cumulative, design_rng.randrange(cumulative[-1]))]
            start = design_rng.randrange(genome['genome_bp'])
            reverse = design_rng.random() < .5
            spans = [(start, min(start+length, genome['genome_bp']))]
            if start+length > genome['genome_bp']:
                spans.append((0, start+length-genome['genome_bp']))
            parts, segments, query_position = [], [], 0
            for begin, end in spans:
                position = begin
                first = max(0, bisect_right(starts, begin)-1)
                for row in truth[first:]:
                    if row['start'] >= end:
                        break
                    lo, hi = max(begin, row['start']), min(end, row['end'])
                    if lo >= hi:
                        continue
                    background = alter(reader.get(position, lo))
                    parts.append(background); query_position += len(background)
                    array = alter(reader.get(lo, hi))
                    sampled_repeat_bp[row['family_id']] += hi-lo
                    if array:
                        segments.append([query_position, query_position+len(array), row['period'], row['family_id'], hi-lo])
                    parts.append(array); query_position += len(array)
                    position = hi
                background = alter(reader.get(position, end))
                parts.append(background); query_position += len(background)
            sequence = ''.join(parts)
            if not sequence:
                raise ValueError('Error model deleted an entire read; input is not usable')
            if reverse:
                sequence = reverse_complement(sequence)
                segments = [[len(sequence)-b, len(sequence)-a, p, f, n] for a, b, p, f, n in segments]
            read_count += 1
            name = f'r{read_count:09d}'
            output.write(f'>{name}\n{sequence}\n')
            audit_writer.writerow([name,start,length,len(sequence),'-' if reverse else '+'])
            for a,b,p,f,n in sorted(segments):
                truth_writer.writerow([name,a,b,p,f,n,n >= 2*p])
                read_array_count += 1
            source_bases += length; observed_bases += len(sequence)
    if observed_bases != source_bases+event_totals[1]-event_totals[2]:
        raise ArithmeticError('Observed sequence length is inconsistent with sampled source and error events')
    result = dict(complete=True, seed=seed, requested_coverage=coverage,
                  actual_source_coverage=source_bases/genome['genome_bp'], actual_observed_coverage=observed_bases/genome['genome_bp'],
                  source_bases=source_bases, total_bases=observed_bases, read_count=read_count, truth_read_segments=read_array_count,
                  sampled_repeat_bp=dict(sampled_repeat_bp), substitution_rate=substitution_rate,
                  insertion_rate=insertion_rate, deletion_rate=deletion_rate,
                  observed_substitutions=event_totals[0], observed_insertions=event_totals[1], observed_deletions=event_totals[2],
                  length_histogram_sha256=digest_file(histogram), genome_manifest_sha256=digest_file(genome_dir/'manifest.json'),
                  script_sha256=digest_file(Path(__file__)), strand_helper_sha256=digest_file(Path(reverse_complement.__code__.co_filename)),
                  warning='uniform_circular_sampling_IID_errors_empirical_lengths_only_not_empirical_HiFi_chemistry;partial_repeat_segments_not_automatically_detectable_arrays',
                  files={p.name: digest_file(p) for p in outdir.iterdir()})
    (outdir/'manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
