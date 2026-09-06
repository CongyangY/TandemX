"""Replay a completed real-input discovery command with frozen current source."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file

PRODUCTS = ('candidate_reads.tsv', 'candidate_monomers.fa', 'monomers.fa', 'families.tsv',
            'monomer_membership.tsv', 'family_similarity.tsv', 'family_audit_summary.json')


def command_with_threads(command: list[str], threads: int | None) -> tuple[list[str], int, int]:
    """Return a copied command with an optional, validated discover thread budget."""
    updated = command[:]
    if updated.count('--threads') != 1:
        raise ValueError('Baseline must declare exactly one --threads value')
    index = updated.index('--threads') + 1
    try:
        baseline_threads = int(updated[index])
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError('Baseline has an invalid --threads value') from exc
    if not 1 <= baseline_threads <= 64:
        raise ValueError('Baseline has an invalid --threads value')
    replay_threads = baseline_threads if threads is None else threads
    if not isinstance(replay_threads, int) or not 1 <= replay_threads <= 64:
        raise ValueError('Threads must be an integer in [1,64]')
    updated[index] = str(replay_threads)
    return updated, baseline_threads, replay_threads


def replay(previous_run: Path, outdir: Path, timeout: float = 1800, profile: bool = False,
           threads: int | None = None) -> dict:
    root = Path(__file__).resolve().parents[2]
    previous_run, outdir = previous_run.resolve(), outdir.resolve()
    old_environment = json.loads((previous_run/'environment.json').read_text())
    old_execution = json.loads((previous_run/'tandemx/execution.json').read_text())
    command = old_execution['command'][:]
    if old_execution['exit_code'] or old_execution['timed_out'] or command[1:4] != ['-m', 'tandemx.cli', 'discover']:
        raise ValueError('Require a successful TandemX discovery baseline')
    if command.count('--reads') != 1 or command.count('--outdir') != 1:
        raise ValueError('Baseline must declare exactly one read input and output')
    command, baseline_threads, replay_threads = command_with_threads(command, threads)
    fasta = Path(command[command.index('--reads')+1])
    input_hash = digest_file(fasta)
    if input_hash != old_environment['input']['fasta_sha256']:
        raise ValueError('Baseline FASTA no longer matches its input receipt')
    baseline = previous_run/'tandemx/discover'
    expected = {name: digest_file(baseline/name) for name in PRODUCTS}
    outdir.mkdir(parents=True, exist_ok=False)
    snapshot = outdir/'source_snapshot'
    environment = source_manifest(root, snapshot)
    shutil.copyfile(Path(__file__), outdir/'replay_discovery.py')
    command[0] = sys.executable
    command[command.index('--outdir')+1] = str(outdir/'discover')
    if profile:
        command[1:1] = ['-m', 'cProfile', '-o', str(outdir/'discovery.prof')]
    environment.update(previous_run=str(previous_run), baseline_environment_sha256=digest_file(previous_run/'environment.json'),
                       baseline_execution_sha256=digest_file(previous_run/'tandemx/execution.json'),
                       input_sha256=input_hash, expected_product_hashes=expected, profile=profile,
                       baseline_threads=baseline_threads, replay_threads=replay_threads,
                       helper_sha256=digest_file(outdir/'replay_discovery.py'),
                       scope='full live-candidate pipeline parity; one run; optional thread-scaling replay; '
                             'concurrent jobs; profiler adds overhead if enabled')
    (outdir/'environment.json').write_text(json.dumps(environment, indent=2)+'\n')
    result = dict(complete=False, command=command, profile=profile)
    try:
        result['execution'] = run_process(command, outdir/'stdout.log', outdir/'stderr.log', timeout,
                                          {**os.environ, 'PYTHONPATH': str(snapshot)}, snapshot)
        if result['execution']['exit_code'] or result['execution']['timed_out']:
            raise RuntimeError('Discovery replay failed; logs retained')
        result['products'] = {name: dict(expected_sha256=value,
                             observed_sha256=digest_file(outdir/'discover'/name),
                             byte_identical=value == digest_file(outdir/'discover'/name))
                             for name, value in expected.items()}
        if digest_file(fasta) != input_hash or any(digest_file(baseline/n) != h for n, h in expected.items()):
            raise ValueError('Baseline input/product changed during replay')
        if not all(p['byte_identical'] for p in result['products'].values()):
            raise ValueError('Full pipeline products differ from baseline')
        result['complete'] = True
        return result
    except Exception as exc:
        result['error'] = str(exc)
        raise
    finally:
        (outdir/'validation.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-run', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--profile', action='store_true', help='Collect cProfile hotspots; timings then include profiling overhead')
    parser.add_argument('--threads', type=int, help='Optional discover thread budget for output-parity scaling')
    args = parser.parse_args()
    replay(args.previous_run, args.outdir, args.timeout, args.profile, args.threads)
