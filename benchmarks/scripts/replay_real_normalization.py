"""Check disk normalization against retained real outputs, or stress input preparation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import time

from benchmarks.challenge.run import source_manifest
from benchmarks.challenge.schema import digest_file, read_table
from benchmarks.scripts.real_disk import normalize_disk, prepare_disk_input


def replay(receipt: Path, sample_id: str, outdir: Path, previous: Path | None = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    environment = source_manifest(Path(__file__).resolve().parents[2], outdir/'source_snapshot')
    environment['evaluator_hashes'] = {}
    for name in ('replay_real_normalization.py', 'real_disk.py', 'fastq_stream.py'):
        source = Path(__file__).with_name(name)
        target = outdir/'source_snapshot/benchmarks/scripts'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        checksum = digest_file(source)
        shutil.copyfile(source, target)
        environment['evaluator_hashes'][name] = checksum
        if digest_file(target) != checksum:
            raise ValueError('Source changed during snapshot')
    (outdir/'environment.json').write_text(json.dumps(environment, indent=2)+'\n')
    start = time.perf_counter()
    result = dict(complete=False, sample_id=sample_id, sampling_receipt=str(receipt.resolve()),
                  previous_run=str(previous.resolve()) if previous else None, normalized_tools=[],
                  warning='evaluator_check_only_no_new_discovery_or_accuracy_result')
    try:
        sample = prepare_disk_input(receipt, sample_id, outdir/'reads.fa', outdir/'evaluation.sqlite')
        result.update(input=sample, input_preparation_seconds=time.perf_counter()-start)
        if previous is not None:
            prior = json.loads((previous/'environment.json').read_text())['input']
            if any(sample[key] != prior[key] for key in ('sample_id', 'fastq_sha256', 'fasta_sha256', 'read_count', 'total_bases')):
                raise ValueError('Prior real benchmark used a different selected input')
            if digest_file(previous/'reads.fa') != sample['fasta_sha256']:
                raise ValueError('Prior input FASTA does not match its receipt')
            result['previous_environment_sha256'] = digest_file(previous/'environment.json')
            rows = read_table(previous/'summary.tsv')
            if len(rows) != 3 or {r['tool'] for r in rows} != {'tandemx', 'trf', 'tidehunter'}:
                raise ValueError('Require exactly three completed comparator summaries')
            for row in rows:
                if row['normalization'] != 'ok':
                    raise ValueError('Prior normalization is not complete')
                tool = row['tool']
                native = previous/tool/{'tandemx':'discover/candidate_reads.tsv', 'trf':'trf.txt', 'tidehunter':'tidehunter.tsv'}[tool]
                normalized = outdir/(tool+'.normalized.tsv')
                before = time.perf_counter()
                metrics = normalize_disk(tool, native, outdir/'evaluation.sqlite', normalized)
                byte_equal = digest_file(normalized) == digest_file(previous/tool/'normalized_arrays.tsv')
                metric_equal = all(value == float(row[key]) for key, value in metrics.items())
                result['normalized_tools'].append(dict(tool=tool, **metrics,
                    normalization_seconds=time.perf_counter()-before, byte_identical=byte_equal,
                    metrics_identical=metric_equal, native_sha256=digest_file(native),
                    normalized_sha256=digest_file(normalized)))
                if not byte_equal or not metric_equal:
                    raise ValueError('Disk normalization differs from retained previous outputs')
        result['complete'] = True
        return result
    except Exception as exc:
        result['error'] = str(exc)
        raise
    finally:
        result['total_seconds'] = time.perf_counter()-start
        (outdir/'replay_receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sampling-receipt', type=Path, required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--previous-run', type=Path)
    args = parser.parse_args()
    replay(args.sampling_receipt, args.sample_id, args.outdir, args.previous_run)
