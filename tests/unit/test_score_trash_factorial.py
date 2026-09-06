"""Construct parser fixtures only; never present them as comparator results."""
import csv
import json

import pytest

from benchmarks.abundance.stream_genome import ArraySpec, generate
from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.schema import digest_file, read_table
from benchmarks.challenge.trash_adapters import ARRAY_FILES, UNIT_FILES
from benchmarks.scripts.score_trash_factorial import evaluate


def csv_file(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.parametrize('tool', ['trash', 'trash2'])
def test_independent_fixture_scoring_coordinates_and_tamper_rejection(tmp_path, tool):
    genome = tmp_path/'genome'
    manifest = generate([ArraySpec('f01', 31, 6)], 1000, 98001, genome)
    truth = read_table(genome/'truth_copy_number.tsv')[0]
    start, end = int(truth['start']), int(truth['end'])
    sequence = read_fasta(genome/'catalogue.fa')['f01']
    run = tmp_path/'run'
    native = run/'native'
    native.mkdir(parents=True)
    if tool == 'trash':
        arrays = [dict(start=start+1, end=end, **{'fasta.name': 'chr_sim',
                  'most.freq.value.N': 31, 'consensus.primary': sequence})]
        units = [dict(start=start+i*31+2, end=start+(i+1)*31+1, width=31,
                      strand='+', seq=sequence, **{'seq.name': 'chr_sim'}) for i in range(6)]
    else:
        arrays = [dict(start=start+1, end=end, seqID='chr_sim', top_N=31, representative=sequence.lower())]
        units = [dict(start=start+i*31+1, end=start+(i+1)*31, width=31,
                      strand='-', sequence=sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1],
                      seqID='chr_sim') for i in range(6)]
    csv_file(native/ARRAY_FILES[tool], arrays)
    csv_file(native/UNIT_FILES[tool], units)
    if tool == 'trash2':
        csv_file(native/'assembly.fa_repeats_with_seq.csv', units)
    receipt = dict(tool=tool, complete=True, synthetic_parser_test_only=True,
                   input_sha256=manifest['files']['genome.fa'],
                   native_products={name: digest_file(native/name) for name in (ARRAY_FILES[tool], UNIT_FILES[tool])})
    (run/'execution.json').write_text(json.dumps(receipt))
    result = evaluate(genome, run, tmp_path/'evaluation')
    assert result['complete']
    assert all(row['cyclic_monomer_recall'] == 1 for row in result['region_metrics'])
    exact = [row for row in result['unit_coverage_metrics'] if row['offset_from_native_1based'] == (-1 if tool == 'trash' else 0)][0]
    assert exact['base_union_recall'] == exact['base_union_precision'] == 1
    audit = read_table(tmp_path/'evaluation/unit_coordinate_audit.tsv')
    assert int(audit[1 if tool == 'trash' else 2]['strand_adjusted_matches']) == 6
    with (native/ARRAY_FILES[tool]).open('a') as handle:
        handle.write('\n')
    with pytest.raises(ValueError, match='changed'):
        evaluate(genome, run, tmp_path/'tampered')
