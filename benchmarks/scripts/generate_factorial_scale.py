"""Generate one development factorial genome and independent read conditions."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
import shutil

from benchmarks.abundance.stream_genome import ArraySpec, generate
from benchmarks.abundance.stream_reads import length_distribution, sample
from benchmarks.challenge.schema import digest_file


def run(
    config_path: Path,
    histogram: Path,
    outdir: Path,
    seed: int,
    max_output_bases: int,
    split: str = 'development',
) -> None:
    config = json.loads(config_path.read_text())
    seed_groups = config['seeds']
    declared = [value for values in seed_groups.values() for value in values]
    if (
        split not in {'development', 'validation'}
        or split not in seed_groups
        or seed not in seed_groups[split]
        or len(set(declared)) != len(declared)
    ):
        raise ValueError(
            'Only distinct predeclared development or validation seeds are allowed'
        )
    for name in ('periods', 'copies', 'gc_fractions', 'unit_substitution_rates', 'coverages'):
        if not config[name] or len(set(config[name])) != len(config[name]):
            raise ValueError('Factorial axes must be nonempty and contain unique values')
    grid = list(product(config['periods'], config['copies'], config['gc_fractions'], config['unit_substitution_rates']))
    specs = [ArraySpec(f'f{i:03d}', p, c, gc, d) for i,(p,c,gc,d) in enumerate(grid,1)]
    specs += [ArraySpec(**entry) for entry in config.get('additional_arrays', [])]
    for spec in specs:
        spec.validate()
    labels = [model['label'] for model in config['read_error_models']]
    if not labels or any(not isinstance(label, str) or not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError('Error-model labels must be nonempty and distinct')
    lengths, _ = length_distribution(histogram, min(config['genome_bp'], 200_000))
    conditions = []
    for coverage in config['coverages']:
        if not math.isfinite(coverage) or not 0 < coverage <= 100:
            raise ValueError('Invalid factorial coverage')
        for error in config['read_error_models']:
            rates = [error[k] for k in ('substitution_rate','insertion_rate','deletion_rate')]
            if any(not math.isfinite(v) or not 0 <= v < 1 for v in rates):
                raise ValueError('Invalid factorial error model')
            conditions.append(dict(coverage=coverage, **error))
    # Each source base emits at most one retained base and one insertion.
    bound = config['genome_bp'] + sum((2 if c['insertion_rate'] else 1)*
            (math.ceil(config['genome_bp']*c['coverage'])+max(lengths)) for c in conditions)
    if not conditions or bound > max_output_bases or max_output_bases <= 0:
        raise ValueError('Requested simulation exceeds its worst-case sequence-base output budget')
    outdir.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(config_path, outdir/'run_config.json')
    shutil.copyfile(histogram, outdir/'length_histogram.tsv')
    source = outdir/'source_snapshot'
    root = Path(__file__).resolve().parents[2]
    hashes = {}
    for rel in ['benchmarks/scripts/generate_factorial_scale.py','benchmarks/abundance/stream_genome.py',
                'benchmarks/abundance/stream_reads.py','benchmarks/challenge/simulate.py','benchmarks/challenge/schema.py']:
        path = root/rel; target=source/rel; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target); hashes[rel]=digest_file(target)
    receipt = dict(complete=False, seed=seed, split=split, source_sha256=hashes,
                   config_sha256=digest_file(config_path), histogram_sha256=digest_file(histogram),
                   worst_case_sequence_base_budget=bound, family_count=len(specs), genome_bp=config['genome_bp'],
                   conditions_completed=[], validation_used=split == 'validation', heldout_used=False,
                   warning=('one_simulated_genome_not_an_independent_plant;factorial_conditions_share_genome_and_read_starts'
                            if split == 'development' else
                            'frozen_validation_genome_not_an_independent_plant;validation_consumed_when_scored;factorial_conditions_share_genome_and_read_starts'))
    receipt_path = outdir/'generation_receipt.json'
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    generate(specs, config['genome_bp'], seed, outdir/'genome', config['background_gc'])
    for i, condition in enumerate(conditions, 1):
        name = f'condition_{i:03d}'
        result = sample(outdir/'genome', outdir/'reads'/name, histogram, seed=seed+1_000_003,
                        coverage=condition['coverage'], substitution_rate=condition['substitution_rate'],
                        insertion_rate=condition['insertion_rate'], deletion_rate=condition['deletion_rate'])
        receipt['conditions_completed'].append(dict(condition_id=name, **condition, read_count=result['read_count'],
            total_bases=result['total_bases'], manifest_sha256=digest_file(outdir/'reads'/name/'manifest.json')))
        receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        with (outdir/'generation.log').open('a') as log:
            log.write(f'{name}\treads={result["read_count"]}\tbases={result["total_bases"]}\n')
    if any(digest_file(root/rel) != digest for rel,digest in hashes.items()):
        raise ValueError('Generator source changed during execution; retain outputs as incomplete')
    receipt['complete'] = True
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--length-histogram', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--split', choices=('development', 'validation'), default='development')
    parser.add_argument('--max-output-bases', type=int, required=True, help='Worst-case emitted sequence bases, excluding headers/TSV/source code')
    args = parser.parse_args()
    run(args.config, args.length_histogram, args.outdir, args.seed, args.max_output_bases, args.split)
