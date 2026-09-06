"""Observed-coordinate truth views for independent factorial discovery scoring."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays, score_base_coverage
from benchmarks.challenge.schema import ArrayRecord, digest_file, read_table


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str,int]={}
    identifier=None
    with path.open() as handle:
        for line in handle:
            line=line.strip()
            if not line:
                continue
            if line.startswith('>'):
                fields=line[1:].split()
                if not fields or fields[0] in lengths:
                    raise ValueError('Empty or duplicate FASTA identifier')
                identifier=fields[0];lengths[identifier]=0
                if len(lengths)>100000:
                    raise ValueError('Pilot evaluator cap is 100000 reads')
            else:
                if identifier is None or set(line.upper())-set('ACGTN'):
                    raise ValueError('Invalid FASTA sequence')
                lengths[identifier]+=len(line)
    if not lengths or any(n<=0 for n in lengths.values()) or sum(lengths.values())>500000000:
        raise ValueError('Require nonempty reads within the 500-Mb evaluator cap')
    return lengths


@dataclass(frozen=True)
class DiscoveryTruth:
    lengths: dict[str,int]
    all_segments: list[ArrayRecord]
    eligible_segments: list[ArrayRecord]
    excluded_reads: set[str]
    founders: dict[str,str]
    metadata: dict


def load_truth(dataset: Path, condition_id: str) -> DiscoveryTruth:
    generation=json.loads((dataset/'generation_receipt.json').read_text())
    if not generation['complete'] or generation['split']!='development' or generation['heldout_used']:
        raise ValueError('Only complete development inputs are allowed')
    conditions=[r for r in generation['conditions_completed'] if r['condition_id']==condition_id]
    if len(conditions)!=1:
        raise ValueError('Unknown or duplicate condition')
    condition=conditions[0];folder=dataset/'reads'/condition_id
    if digest_file(folder/'manifest.json')!=condition['manifest_sha256']:
        raise ValueError('Read manifest hash mismatch')
    manifest=json.loads((folder/'manifest.json').read_text())
    for name in ('reads.fa','truth_read_segments.tsv'):
        if digest_file(folder/name)!=manifest['files'][name]:
            raise ValueError('Observed input/truth hash mismatch')
    genome=json.loads((dataset/'genome/manifest.json').read_text())
    if digest_file(dataset/'genome/catalogue.fa')!=genome['files']['catalogue.fa']:
        raise ValueError('Founder catalogue hash mismatch')
    founders=read_fasta(dataset/'genome/catalogue.fa')
    lengths=fasta_lengths(folder/'reads.fa')
    if len(lengths)!=manifest['read_count'] or sum(lengths.values())!=manifest['total_bases']:
        raise ValueError('Observed read/base totals differ')
    all_segments,eligible,excluded=[],[],set()
    for row in read_table(folder/'truth_read_segments.tsv'):
        family=row['family_id'];period=int(row['period'])
        if family not in founders or len(founders[family])!=period or row['at_least_two_source_units'] not in {'True','False'}:
            raise ValueError('Unknown founder or invalid truth period/eligibility')
        record=ArrayRecord(row['read_id'],int(row['start']),int(row['end']),period,family_id=family)
        if record.read_id not in lengths or record.end>lengths[record.read_id]:
            raise ValueError('Out-of-bounds observed truth')
        source_eligible=int(row['sampled_source_repeat_bp'])>=2*period
        if source_eligible!=(row['at_least_two_source_units']=='True'):
            raise ValueError('Inconsistent source-unit eligibility')
        all_segments.append(record)
        if source_eligible and record.end-record.start>=100 and 30<=period<=1000:
            eligible.append(record)
        else:
            excluded.add(record.read_id)
    if len(all_segments)!=manifest['truth_read_segments']:
        raise ValueError('Truth segment count mismatch')
    metadata=dict(seed=generation['seed'],condition_id=condition_id,coverage=condition['coverage'],
        error_model=condition['label'],generation_sha256=digest_file(dataset/'generation_receipt.json'),
        condition_manifest_sha256=condition['manifest_sha256'],reads_sha256=manifest['files']['reads.fa'],
        truth_sha256=manifest['files']['truth_read_segments.tsv'],catalogue_sha256=genome['files']['catalogue.fa'],
        observed_read_count=len(lengths),observed_bases=sum(lengths.values()),
        all_truth_segments=len(all_segments),eligible_truth_segments=len(eligible),
        excluded_partial_or_out_of_scope_reads=len(excluded),all_planted_families=len(founders),
        observed_eligible_families=len({r.family_id for r in eligible}))
    return DiscoveryTruth(lengths,all_segments,eligible,excluded,founders,metadata)


def score_predictions(predicted: list[ArrayRecord], truth: DiscoveryTruth) -> tuple[dict,list[dict]]:
    for record in predicted:
        if record.read_id not in truth.lengths or record.end>truth.lengths[record.read_id]:
            raise ValueError('Prediction outside observed input')
    safe_lengths={r:n for r,n in truth.lengths.items() if r not in truth.excluded_reads}
    safe_truth=[r for r in truth.eligible_segments if r.read_id in safe_lengths]
    safe_predictions=[r for r in predicted if r.read_id in safe_lengths]
    if not safe_lengths:
        raise ValueError('No reads remain for the declared eligible-array endpoint')
    array_metrics,details=score_arrays(safe_predictions,safe_truth,safe_lengths,.5)
    # Every planted observed base, including short fragments, stays in the
    # full-input endpoint. Array scoring excludes reads by truth-only criteria.
    metrics={f'eligible_read_{k}':v for k,v in array_metrics.items()}
    metrics.update({f'all_input_{k}':v for k,v in score_base_coverage(predicted,truth.all_segments).items()})
    metrics.update(all_input_prediction_count=len(predicted),
                   predictions_on_excluded_reads=len(predicted)-len(safe_predictions))
    return metrics,details
