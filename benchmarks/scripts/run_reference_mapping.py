"""Run a source-pinned HiFi mapping QC on a completed whole-library sample."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.qc_reference_fasta import reference_qc
from benchmarks.scripts.reference_mapping_qc import prepare_queries, load_alignments, summarize_mapping


def run(reference: Path, reference_sha256: str, receipt: Path, sample_id: str,
        minimap2: Path, organelles: set[str], outdir: Path, timeout: int = 3600) -> None:
    reference, receipt, minimap2, outdir = [p.resolve() for p in (reference, receipt, minimap2, outdir)]
    if len(reference_sha256) != 64 or set(reference_sha256)-set('0123456789abcdef') or timeout <= 0:
        raise ValueError('Require exact reference SHA-256 and positive timeout')
    sampling = json.loads(receipt.read_text())
    selected = [s for s in sampling['samples'] if s['sample_id'] == sample_id]
    if not sampling.get('complete') or len(selected) != 1 or selected[0]['status'] != 'ok':
        raise ValueError('Require completed nonempty whole-library sample')
    fastq = receipt.parent/(sample_id+'.fastq.gz')
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir/'source_snapshot'
    source = source_manifest(root, snapshot)
    helpers = [Path(__file__), *[Path(__file__).with_name(name) for name in
               ('reference_mapping_qc.py', 'qc_reference_fasta.py', 'fastq_stream.py', 'real_disk.py')]]
    helper_hashes = {}
    for path in helpers:
        target = snapshot/path.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        helper_hashes[str(path)] = digest_file(target)
    inputs = {str(receipt): digest_file(receipt), str(reference): reference_sha256,
              str(fastq): selected[0]['fastq_sha256'], str(minimap2): digest_file(minimap2)}
    version = subprocess.run([str(minimap2), '--version'], check=True, text=True, capture_output=True).stdout.strip()
    # One reference batch is required for meaningful MAPQ; do not allow a
    # hidden split index. Larger assemblies need a separately profiled design.
    command = [str(minimap2), '-x', 'map-hifi', '-c', '--cs=short', '--secondary=yes',
               '-N', '5', '-t', '1', '-K', '50M', '-I', '3G', str(reference), str(fastq)]
    provenance = dict(source=source, helper_sha256=helper_hashes, input_sha256=inputs,
        minimap2_version=version, command=command, organelle_contigs=sorted(organelles),
        source_documentation='https://lh3.github.io/minimap2/minimap2.html',
        warning='conditional_reference_concordance_not_biological_truth;one_reference_batch;concurrent_diagnostic_resources')
    (outdir/'environment.json').write_text(json.dumps(provenance, indent=2)+'\n')
    validation = dict(complete=False)
    try:
        ref = reference_qc(reference)
        if ref['input_sha256'] != reference_sha256 or ref['total_bases'] > 3_000_000_000:
            raise ValueError('Reference hash mismatch or exceeds single-index 3-Gb QC limit')
        lengths = {r['contig']: r['length_bp'] for r in ref['contigs']}
        if organelles-set(lengths):
            raise ValueError('Unknown organellar contig')
        (outdir/'reference_qc.json').write_text(json.dumps(ref, indent=2)+'\n')
        database = outdir/'mapping.sqlite'
        observed = prepare_queries(fastq, database, selected[0])
        measured = run_process(command, outdir/'alignments.paf', outdir/'minimap2.stderr.log', timeout)
        (outdir/'execution.json').write_text(json.dumps(dict(command=command, **measured), indent=2)+'\n')
        if measured['exit_code'] != 0 or measured['timed_out']:
            raise RuntimeError('Mapping failed; partial PAF must not be interpreted as low concordance')
        paf = load_alignments(outdir/'alignments.paf', database, lengths)
        summary = summarize_mapping(database, lengths, organelles, outdir)
        if summary['read_count'] != observed['read_count'] or summary['total_bases'] != observed['total_bases']:
            raise ValueError('Mapping summary lost query denominators')
        if any(digest_file(Path(name)) != value for name, value in {**inputs, **helper_hashes}.items()):
            raise ValueError('Inputs or QC source changed during run')
        validation.update(complete=True, **observed, **paf, summary=summary,
            outputs={p.name: digest_file(p) for p in outdir.iterdir() if p.suffix in {'.tsv', '.paf'}},
            evidence_limit='unmapped_not_proven_contamination;reference_collapse_and_donor_divergence_unresolved')
    except Exception as exc:
        validation['error'] = str(exc)
        raise
    finally:
        (outdir/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--reference-sha256', required=True)
    parser.add_argument('--sampling-receipt', type=Path, required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--minimap2', type=Path, required=True)
    parser.add_argument('--organelle-contig', action='append', default=[])
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--timeout', type=int, default=3600)
    args = parser.parse_args()
    run(args.reference, args.reference_sha256, args.sampling_receipt, args.sample_id,
        args.minimap2, set(args.organelle_contig), args.outdir, args.timeout)
