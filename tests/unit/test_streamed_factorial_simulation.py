import io
import json
from pathlib import Path
import random

import pytest

from benchmarks.abundance.stream_genome import ArraySpec, FixedFastaReader, FixedFastaWriter, generate
from benchmarks.abundance.stream_reads import alter_segment, length_distribution, sample
from benchmarks.challenge.schema import read_table
from benchmarks.challenge.simulate import mutate, reverse_complement
from benchmarks.scripts.generate_factorial_scale import run


def fasta(path: Path) -> dict[str, str]:
    result = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:]; result[name] = ''
        else:
            result[name] += line
    return result


def test_fixed_width_random_access_at_every_boundary_with_chunked_writes():
    sequence = 'ACGT'*53+'G'
    output = io.BytesIO()
    writer = FixedFastaWriter(output, width=13)
    for start in range(0, len(sequence), 7):
        writer.add(sequence[start:start+7])
    index = writer.finish()
    reader = FixedFastaReader(output, index)
    for start in range(len(sequence)+1):
        for length in [0,1,12,13,14,100]:
            end = min(len(sequence), start+length)
            assert reader.get(start,end) == sequence[start:end]
    with pytest.raises(ValueError, match='bounds'):
        reader.get(0,len(sequence)+1)


def test_streamed_genome_truth_founders_background_and_divergence(tmp_path):
    specs = [ArraySpec('a',11,8,.3,0), ArraySpec('b',23,10,.7,.1)]
    for name in ['one','same']:
        generate(specs, 3000, 6311, tmp_path/name)
    assert (tmp_path/'one/genome.fa').read_bytes() == (tmp_path/'same/genome.fa').read_bytes()
    sequence = fasta(tmp_path/'one/genome.fa')['chr_sim']
    motifs = fasta(tmp_path/'one/catalogue.fa')
    truth = read_table(tmp_path/'one/truth_copy_number.tsv')
    assert len(sequence) == 3000
    for row in truth:
        observed = sequence[int(row['start']):int(row['end'])]
        original = motifs[row['family_id']]*int(row['copies'])
        assert len(observed) == int(row['repeat_bp']) == len(original)
        assert sum(a != b for a,b in zip(observed,original)) == int(row['observed_unit_substitutions'])
    assert truth[0]['observed_unit_substitutions'] == '0'
    assert int(truth[1]['observed_unit_substitutions']) > 0
    generate([ArraySpec('a',11,8,.3,.2), ArraySpec('b',23,10,.7,.2)],3000,6311,tmp_path/'changed')
    assert (tmp_path/'one/catalogue.fa').read_bytes() == (tmp_path/'changed/catalogue.fa').read_bytes()
    other = fasta(tmp_path/'changed/genome.fa')['chr_sim']
    excluded = {i for row in truth for i in range(int(row['start']),int(row['end']))}
    assert all(a == b for i,(a,b) in enumerate(zip(sequence,other)) if i not in excluded)


def test_error_events_match_independent_mutator_and_sequence_length_identity():
    sequence = 'ACGT'*300
    for sub,ins,delete in [(0,0,0),(.03,0,0),(.01,.03,.05),(0,.9,.9)]:
        observed, events = alter_segment(sequence, random.Random(6321), sub,ins,delete)
        assert observed == mutate(sequence, random.Random(6321), sub,ins,delete)
        assert len(observed) == len(sequence)+events[1]-events[2]
        assert 0 <= events[0] <= len(sequence)-events[2]


def test_sampled_source_occupancy_and_observed_coordinates_have_independent_oracles(tmp_path):
    generate([ArraySpec('a',11,8),ArraySpec('b',23,7)],1000,6331,tmp_path/'genome')
    hist = tmp_path/'hist.tsv';hist.write_text('length_bp\tread_count\n150\t2\n280\t1\n')
    source = fasta(tmp_path/'genome/genome.fa')['chr_sim']
    truth = read_table(tmp_path/'genome/truth_copy_number.tsv')
    receipts = []
    for name,rates in [('clean',(0,0,0)),('repeat',(0,0,0)),('error',(.01,.03,.05))]:
        receipts.append(sample(tmp_path/'genome',tmp_path/name,hist,seed=6341,coverage=30,
                        substitution_rate=rates[0],insertion_rate=rates[1],deletion_rate=rates[2]))
    assert (tmp_path/'clean/reads.fa').read_bytes() == (tmp_path/'repeat/reads.fa').read_bytes()
    coords = read_table(tmp_path/'clean/sampling.tsv')
    error_coords = read_table(tmp_path/'error/sampling.tsv')
    assert [(r['genome_start'],r['source_length'],r['strand']) for r in coords] == [(r['genome_start'],r['source_length'],r['strand']) for r in error_coords]
    observed = fasta(tmp_path/'clean/reads.fa')
    for row in coords:
        start,length = int(row['genome_start']),int(row['source_length'])
        expected = ''.join(source[(start+i)%len(source)] for i in range(length))
        if row['strand'] == '-': expected = reverse_complement(expected)
        assert observed[row['read_id']] == expected
    assert {r['strand'] for r in coords} == {'+','-'}
    assert any(int(r['genome_start'])+int(r['source_length']) > len(source) for r in coords)
    segments = read_table(tmp_path/'clean/truth_read_segments.tsv')
    for row in truth:
        positions = set(range(int(row['start']),int(row['end'])))
        expected = sum((int(r['genome_start'])+i)%len(source) in positions for r in coords for i in range(int(r['source_length'])))
        assert receipts[0]['sampled_repeat_bp'][row['family_id']] == expected
        assert sum(int(s['end'])-int(s['start']) for s in segments if s['family_id']==row['family_id']) == expected
    result = receipts[2]
    assert result['total_bases'] == result['source_bases']+result['observed_insertions']-result['observed_deletions']
    error_reads = fasta(tmp_path/'error/reads.fa')
    for row in read_table(tmp_path/'error/truth_read_segments.tsv'):
        assert 0 <= int(row['start']) < int(row['end']) <= len(error_reads[row['read_id']])
    # An independent per-query-base oracle verifies strand-transformed label masks.
    for row in coords:
        start,length = int(row['genome_start']),int(row['source_length'])
        mask = [any(int(t['start']) <= (start+i)%len(source) < int(t['end']) for t in truth) for i in range(length)]
        if row['strand']=='-': mask.reverse()
        predicted = {i for s in segments if s['read_id']==row['read_id'] for i in range(int(s['start']),int(s['end']))}
        assert predicted == {i for i,v in enumerate(mask) if v}


def test_controller_budgets_reserved_seeds_and_actual_small_generation(tmp_path):
    hist=tmp_path/'hist.tsv';hist.write_text('length_bp\tread_count\n50\t1\n')
    config=dict(seeds={'development':[6351],'heldout':[7351]},genome_bp=1000,background_gc=.4,
        periods=[11],copies=[3,7],gc_fractions=[.3,.7],unit_substitution_rates=[0,.1],
        coverages=[1],read_error_models=[dict(label='clean',substitution_rate=0,insertion_rate=0,deletion_rate=0)])
    path=tmp_path/'config.json';path.write_text(json.dumps(config))
    for seed,budget,message in [(7351,100000,'development'),(6351,1,'budget')]:
        with pytest.raises(ValueError,match=message):run(path,hist,tmp_path/'invalid',seed,budget)
        assert not (tmp_path/'invalid').exists()
    run(path,hist,tmp_path/'ok',6351,100000)
    receipt=json.loads((tmp_path/'ok/generation_receipt.json').read_text())
    assert receipt['complete'] and receipt['family_count']==8 and not receipt['heldout_used']
    assert len(receipt['conditions_completed'])==1
    hist.write_text('length_bp\tread_count\n50\t1\n50\t2\n')
    with pytest.raises(ValueError,match='duplicate'):length_distribution(hist,100)
