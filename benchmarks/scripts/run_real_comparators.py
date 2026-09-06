"""Task-matched discovery pilot on a verified whole-library random sample.

Observed calls/coverage are descriptive, not accuracy without external truth.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import random
import shutil
import sys

from benchmarks.challenge.adapters import build_command, parse_arrays
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.fastq_stream import hashed_fastq, records


def prepare_input(receipt: Path, sample_id: str, output: Path) -> tuple[dict, dict[str, int]]:
    sampling = json.loads(receipt.read_text())
    rows = [r for r in sampling['samples'] if r['sample_id'] == sample_id]
    if not sampling.get('complete') or len(rows) != 1 or rows[0]['status'] != 'ok':
        raise ValueError('Require one completed nonempty whole-library sample')
    sample = rows[0]
    if sample['read_count'] > 100_000:
        raise ValueError('This pilot has a 100,000-read controller cap; larger scale requires a disk-backed evaluator')
    path = receipt.parent/(sample_id+'.fastq.gz')
    lengths = {}
    with hashed_fastq(path) as (handle, digest), output.open('wb') as fasta:
        for record in records(handle):
            identifier = record.identifier.decode('ascii')
            if identifier in lengths:
                raise ValueError('Duplicate read ID in selected input')
            lengths[identifier] = len(record.sequence)
            fasta.write(b'>'+record.header[1:]+b'\n'+record.sequence+b'\n')
    if digest.hexdigest() != sample['fastq_sha256'] or len(lengths) != sample['read_count'] or sum(lengths.values()) != sample['total_bases']:
        raise ValueError('Selected FASTQ does not match its sampling receipt')
    return dict(**sample, sampling_receipt_sha256=digest_file(receipt), fasta_sha256=digest_file(output)), lengths


def describe_arrays(arrays, lengths: dict[str, int]) -> dict:
    grouped = defaultdict(list)
    for array in arrays:
        if array.read_id not in lengths or not 0 <= array.start < array.end <= lengths[array.read_id]:
            raise ValueError('Native prediction has an unknown read or invalid coordinates')
        grouped[array.read_id].append((array.start, array.end))
    covered = 0
    for intervals in grouped.values():
        end = -1
        for start, stop in sorted(intervals):
            covered += max(0, stop-max(start, end))
            end = max(end, stop)
    return dict(observed_in_scope_calls=len(arrays), observed_positive_reads=len(grouped),
                observed_union_bp=covered, observed_union_base_fraction=covered/sum(lengths.values()))


def run(receipt: Path, sample_id: str, outdir: Path, trf: Path, tidehunter: Path, timeout: float) -> None:
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    sample, lengths = prepare_input(receipt.resolve(), sample_id, outdir/'reads.fa')
    snapshot = outdir/'source_snapshot'
    root = Path(__file__).resolve().parents[2]
    manifest = source_manifest(root, snapshot)
    for name in ('run_real_comparators.py', 'fastq_stream.py'):
        source = Path(__file__).with_name(name)
        target = snapshot/'benchmarks/scripts'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    tools = dict(tandemx=sys.executable, trf=str(trf.resolve()), tidehunter=str(tidehunter.resolve()))
    if any(not Path(path).is_file() for path in tools.values()) or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Missing executable or invalid timeout')
    order = list(tools)
    random.Random(6101).shuffle(order)
    manifest.update(input=sample, tool_paths=tools, tool_hashes={t: digest_file(Path(p)) for t, p in tools.items()},
                    tool_order=order, scope=dict(min_period=30, max_period=1000, min_span=100),
                    repetitions=1, threads=1, timeout_per_tool_seconds=timeout,
                    runner_sha256=digest_file(Path(__file__)), parser_sha256=digest_file(Path(__file__).with_name('fastq_stream.py')),
                    accuracy='not_assessed_without_curated_independent_truth',
                    resources='direct-child wait4; excludes controller; acquisition may overlap; not publication ranking')
    (outdir/'environment.json').write_text(json.dumps(manifest, indent=2)+'\n')
    summaries = []
    for tool in order:
        folder = outdir/tool
        folder.mkdir()
        command, output = build_command(tool, tools[tool], outdir/'reads.fa', folder, 30, 1000, 100)
        if tool == 'tandemx':
            command[0:1] = [sys.executable, '-m', 'tandemx.cli']
            command += ['--discovery-method', 'elastic', '--clustering-method', 'sequence', '--cluster-identity', '.95']
        measured = run_process(command, output if tool == 'trf' else folder/'stdout.log', folder/'stderr.log', timeout,
                               {**os.environ, 'PYTHONPATH': str(snapshot)} if tool == 'tandemx' else None,
                               snapshot if tool == 'tandemx' else folder)
        (folder/'execution.json').write_text(json.dumps(dict(command=command, **measured), indent=2)+'\n')
        row = dict(tool=tool, **measured, observed_in_scope_calls=None, observed_positive_reads=None,
                   observed_union_bp=None, observed_union_base_fraction=None, normalization='not_run',
                   warning='descriptive_real_input_pilot_no_accuracy_truth_no_resource_ranking')
        if measured['exit_code'] == 0 and not measured['timed_out']:
            try:
                arrays = parse_arrays(tool, output, 30, 1000, 100)
                row.update(describe_arrays(arrays, lengths), normalization='ok')
                write_table(folder/'normalized_arrays.tsv', [asdict(a) for a in arrays],
                            ['read_id', 'start', 'end', 'period', 'sequence', 'family_id'])
            except Exception as exc:
                row['normalization'] = 'failed: '+str(exc)
        summaries.append(row)
        write_table(outdir/'summary.tsv', summaries, list(row))
    if any(r['normalization'] != 'ok' for r in summaries):
        raise RuntimeError('Real pilot has a failed execution or normalization; inspect retained receipts')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sampling-receipt', type=Path, required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--trf', type=Path, required=True)
    parser.add_argument('--tidehunter', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=900)
    args = parser.parse_args()
    run(args.sampling_receipt, args.sample_id, args.outdir, args.trf, args.tidehunter, args.timeout)
