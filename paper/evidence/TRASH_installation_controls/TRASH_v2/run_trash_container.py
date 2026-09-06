"""Execute pinned TRASH assembly comparators offline; retain native outputs/failures."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import uuid

from benchmarks.challenge.run import run_process
from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.qc_reference_fasta import reference_qc


def tool_command(tool: str) -> list[str]:
    if tool == 'trash':
        return ['bash', '/opt/TRASH/TRASH_run.sh', '/input/assembly.fa', '--def',
                '--o', '/output', '--par', '1', '--randomseed', '6101']
    if tool == 'trash2':
        return ['Rscript', '--vanilla', '/opt/TRASH_2/src/TRASH.R', '-f', '/input/assembly.fa',
                '-o', '/output', '-p', '1']
    raise ValueError('Tool must be trash or trash2')


MEASURE = r'''printf 'elapsed_seconds\tmax_rss_kib\tuser_seconds\tsystem_seconds\texit_code\n' > /output/linux_time.tsv
/usr/bin/time -a -o /output/linux_time.tsv -f '%e\t%M\t%U\t%S\t%x' "$@"
code=$?
cat /sys/fs/cgroup/memory.peak > /output/cgroup_memory_peak_bytes.txt
exit "$code"
'''


def linux_time(path: Path) -> dict:
    lines = path.read_text().splitlines()
    expected = ['elapsed_seconds', 'max_rss_kib', 'user_seconds', 'system_seconds', 'exit_code']
    if not lines or lines[0].split('\t') != expected:
        raise ValueError('Invalid GNU time table header')
    # GNU time writes an extra status line for nonzero exits. Keep it in the
    # native file and parse exactly one numeric row; never turn failure into zero.
    values = [line.split('\t') for line in lines[1:] if len(line.split('\t')) == 5]
    if len(values) != 1:
        raise ValueError('Missing or ambiguous GNU time record')
    result = dict(zip(expected, map(float, values[0])))
    if any(not math.isfinite(v) or v < 0 for v in result.values()):
        raise ValueError('Invalid GNU time numeric measurement')
    if not result['exit_code'].is_integer() or not result['max_rss_kib'].is_integer():
        raise ValueError('GNU time exit/RSS must be integers')
    result['peak_rss_mib'] = result['max_rss_kib']/1024
    return result


def native_products(tool: str) -> tuple[str, tuple[str, ...]]:
    if tool == 'trash':
        return ('TRASH finished, exiting', ('all.repeats.from.assembly.fa.csv',
                'Summary.of.repetitive.regions.assembly.fa.csv'))
    if tool == 'trash2':
        return ('TRASH exiting correctly', ('assembly.fa_arrays.csv', 'assembly.fa_repeats.csv'))
    raise ValueError('Unsupported comparator')


def checked(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}): {command!r}\n{result.stderr}')
    return result.stdout


def run(fasta: Path, outdir: Path, tool: str, image: str, timeout: float = 1800) -> dict:
    command = tool_command(tool)
    fasta, outdir = fasta.resolve(), outdir.resolve()
    if fasta.suffix not in {'.fa', '.fna', '.fasta'} or timeout <= 0:
        raise ValueError('Require uncompressed FASTA and positive timeout')
    if any(',' in str(p) or '\n' in str(p) for p in (fasta, outdir)):
        raise ValueError('Docker bind paths cannot contain commas or newlines')
    # Resolve tags before inference; all execution uses the immutable image ID.
    info = json.loads(checked(['docker', 'image', 'inspect', image]))[0]
    if info['Os'] != 'linux' or info['Architecture'] not in {'arm64', 'amd64'}:
        raise ValueError('Expected a Linux comparator image with declared architecture')
    qc = reference_qc(fasta)
    outdir.mkdir(parents=True, exist_ok=False)
    native = outdir/'native'
    native.mkdir()
    name = 'tandemx-benchmark-'+uuid.uuid4().hex
    call = ['docker', 'run', '--name', name, '--network', 'none', '--cpus', '1', '--memory', '8g',
            '--mount', f'type=bind,source={fasta},target=/input/assembly.fa,readonly',
            '--mount', f'type=bind,source={native},target=/output',
            '-w', '/opt/TRASH_2/src' if tool == 'trash2' else '/work', info['Id'],
            'bash', '-c', MEASURE, 'tandemx-measure', *command]
    (outdir/'image_inspect.json').write_text(json.dumps(info, indent=2)+'\n')
    (outdir/'input_qc.json').write_text(json.dumps(qc, indent=2)+'\n')
    shutil.copyfile(Path(__file__), outdir/'run_trash_container.py')
    receipt = dict(complete=False, tool=tool, command=call, image_id=info['Id'],
                   script_sha256=digest_file(outdir/'run_trash_container.py'),
                   input_sha256=qc['input_sha256'], platform=info['Os']+'/'+info['Architecture'],
                   scope='native_default_assembly_annotation;one_core;no_HOR;no_accuracy_scoring',
                   resource_note='Linux GNU time max RSS is not sum of concurrent processes; cgroup peak includes cache; host Docker client is separate')
    try:
        receipt['host_docker_client'] = run_process(call, outdir/'stdout.log', outdir/'stderr.log', timeout)
        state = json.loads(checked(['docker', 'inspect', name]))[0]['State']
        if state['Running']:
            checked(['docker', 'kill', name])
            state = json.loads(checked(['docker', 'inspect', name]))[0]['State']
        receipt['container_state'] = state
        checked(['docker', 'cp', name+':/opt/provenance', str(outdir/'container_provenance')])
        measured = native/'linux_time.tsv'
        if measured.exists():
            receipt['linux_command'] = linux_time(measured)
        peak = native/'cgroup_memory_peak_bytes.txt'
        receipt['cgroup_memory_peak_bytes'] = int(peak.read_text().strip()) if peak.exists() else None
        if (receipt['host_docker_client']['timed_out'] or receipt['host_docker_client']['exit_code']
                or state['ExitCode'] or state['OOMKilled']):
            raise RuntimeError('Comparator failed or timed out; native logs/resources retained')
        if receipt.get('linux_command', {}).get('exit_code') != 0:
            raise ValueError('No successful in-container resource receipt')
        marker, products = native_products(tool)
        if marker not in (outdir/'stdout.log').read_text():
            raise ValueError('Missing native completion marker')
        for filename in products:
            with (native/filename).open(newline='') as handle:
                if not next(csv.reader(handle), None):
                    raise ValueError('Missing native output header: '+filename)
        if digest_file(fasta) != qc['input_sha256']:
            raise ValueError('Input FASTA changed during inference')
        receipt['native_products'] = {f: digest_file(native/f) for f in products}
        receipt['complete'] = True
        return receipt
    except Exception as exc:
        receipt['error'] = str(exc)
        raise
    finally:
        # This name was generated for this run; do not touch unrelated containers.
        cleanup = subprocess.run(['docker', 'rm', '-f', name], capture_output=True, text=True, check=False)
        receipt['cleanup'] = dict(exit_code=cleanup.returncode, stdout=cleanup.stdout, stderr=cleanup.stderr)
        (outdir/'execution.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--tool', choices=('trash', 'trash2'), required=True)
    parser.add_argument('--image', required=True)
    parser.add_argument('--timeout', type=float, default=1800)
    args = parser.parse_args()
    run(args.fasta, args.outdir, args.tool, args.image, args.timeout)
