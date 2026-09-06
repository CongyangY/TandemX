"""Run two clustering implementations in fresh children with exact output checks."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.replay_clustering import load_candidates


def worker(candidates: Path, module_path: Path, output: Path) -> dict:
    rows = load_candidates(candidates)
    spec = importlib.util.spec_from_file_location('isolated_clustering', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started = time.perf_counter()
    families, membership = module.cluster_monomers(rows, 1, .95, 'rust')
    elapsed = time.perf_counter()-started
    # Stream JSON encoding so a second entire serialized copy does not dominate RSS.
    with output.open('x') as handle:
        json.dump(dict(families=[asdict(f) for f in families], membership=membership),
                  handle, sort_keys=True, separators=(',', ':'))
    result = dict(candidate_count=len(rows), family_count=len(families), clustering_seconds=elapsed,
                  output_sha256=digest_file(output), clustering_source_sha256=digest_file(module_path))
    output.with_suffix('.receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def replay(candidate_run: Path, baseline_run: Path, outdir: Path, timeout: float = 3600) -> dict:
    root = Path(__file__).resolve().parents[2]
    candidate_run, baseline_run, outdir = candidate_run.resolve(), baseline_run.resolve(), outdir.resolve()
    old = json.loads((baseline_run/'environment.json').read_text())
    snapshot = baseline_run/'source_snapshot'
    # Clustering is the only core source allowed to differ. This also fixes native
    # kernels, distance rules, dataclasses, orientation and low-complexity helpers.
    for relative, expected in old['file_hashes'].items():
        if not relative.startswith(('tandemx/', 'rust-core/')):
            continue
        if digest_file(snapshot/relative) != expected:
            raise ValueError('Baseline source snapshot changed: '+relative)
        if relative != 'tandemx/discover/clustering.py' and digest_file(root/relative) != expected:
            raise ValueError('A different core module changed: '+relative)
    outdir.mkdir(parents=True, exist_ok=False)
    environment = source_manifest(root, outdir/'source_snapshot')
    helpers = {}
    for name in ('replay_clustering_isolated.py', 'replay_clustering.py'):
        src = Path(__file__).with_name(name)
        dst = outdir/'source_snapshot/benchmarks/scripts'/name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        helpers[name] = digest_file(dst)
        if helpers[name] != digest_file(src):
            raise ValueError('Replay helper changed during snapshot')
    package = root/'benchmarks/__init__.py'
    shutil.copyfile(package, outdir/'source_snapshot/benchmarks/__init__.py')
    helpers['benchmarks/__init__.py'] = digest_file(package)
    folder = candidate_run/'tandemx/discover'
    environment.update(baseline_run=str(baseline_run), candidate_run=str(candidate_run),
                       baseline_clustering_sha256=digest_file(snapshot/'tandemx/discover/clustering.py'),
                       helper_hashes=helpers,
                       input_hashes={name: digest_file(folder/name) for name in ('candidate_reads.tsv', 'candidate_monomers.fa')},
                       resources='fresh-child wait4 peak includes import/input loading/clustering/output encoding; concurrent jobs possible',
                       scope='serialized rounded candidates; not full-pipeline parity; one repetition in fixed order')
    (outdir/'environment.json').write_text(json.dumps(environment, indent=2)+'\n')
    result = dict(complete=False, measurements=[], exact_output_parity=False)
    try:
        variants = [('baseline', snapshot/'tandemx/discover/clustering.py'),
                    ('packed', outdir/'source_snapshot/tandemx/discover/clustering.py')]
        for label, module_path in variants:
            output = outdir/(label+'.json')
            command = [sys.executable, '-m', 'benchmarks.scripts.replay_clustering_isolated',
                       '--worker-candidates', str(folder), '--worker-module', str(module_path), '--worker-output', str(output)]
            measured = run_process(command, outdir/(label+'.stdout.log'), outdir/(label+'.stderr.log'), timeout,
                                   {**os.environ, 'PYTHONPATH': str(outdir/'source_snapshot')}, outdir/'source_snapshot')
            entry = dict(label=label, command=command, **measured)
            result['measurements'].append(entry)
            if measured['exit_code'] or measured['timed_out']:
                raise RuntimeError('Clustering child failed; logs and outputs retained')
            entry.update(json.loads(output.with_suffix('.receipt.json').read_text()))
        result['exact_output_parity'] = result['measurements'][0]['output_sha256'] == result['measurements'][1]['output_sha256']
        if not result['exact_output_parity']:
            raise ValueError('Clustering outputs differ')
        result['complete'] = True
        return result
    except Exception as exc:
        result['error'] = str(exc)
        raise
    finally:
        (outdir/'validation.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-run', type=Path)
    parser.add_argument('--baseline-run', type=Path)
    parser.add_argument('--outdir', type=Path)
    parser.add_argument('--timeout', type=float, default=3600)
    parser.add_argument('--worker-candidates', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--worker-module', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--worker-output', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_candidates and args.worker_module and args.worker_output:
        worker(args.worker_candidates, args.worker_module, args.worker_output)
    elif args.candidate_run and args.baseline_run and args.outdir:
        replay(args.candidate_run, args.baseline_run, args.outdir, args.timeout)
    else:
        parser.error('Require candidate run, baseline run and a new output directory')
