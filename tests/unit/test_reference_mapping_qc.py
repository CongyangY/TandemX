import csv
import hashlib

import pytest

from benchmarks.scripts.reference_mapping_qc import prepare_queries, cigar_blocks, load_alignments, summarize_mapping


def indexed_reads(tmp_path):
    path = tmp_path/'reads.fq'
    path.write_text('@a\n'+'ACGT'*5+'\n+\n'+'I'*20+'\n@b\n'+'G'*20+'\n+\n'+'I'*20+'\n@c\n'+'A'*20+'\n+\n'+'I'*20+'\n')
    db = tmp_path/'mapping.sqlite'
    prepare_queries(path, db, dict(read_count=3, total_bases=60, fastq_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return db


def test_cigar_target_blocks_exclude_deletions_and_check_both_axes():
    assert cigar_blocks('5M2I3M4D2M', 12, 20, 34) == [(20, 25), (25, 28), (32, 34)]
    assert cigar_blocks('3=2X', 5, 0, 5) == [(0, 3), (3, 5)]
    for cigar in ['5M0I5M', '5S5M', '', '10Mtrailing']:
        with pytest.raises(ValueError):
            cigar_blocks(cigar, 10, 0, 10)
    with pytest.raises(ValueError, match='consumption'):
        cigar_blocks('10M', 9, 0, 10)


def test_paf_unions_include_unmapped_denominators_and_missing_mapq_is_not_high(tmp_path):
    db = indexed_reads(tmp_path)
    paf = tmp_path/'out.paf'
    paf.write_text('a\t20\t0\t12\t+\tchr\t100\t0\t14\t10\t16\t30\ttp:A:P\tcg:Z:5M2I3M4D2M\n'
                   'a\t20\t8\t18\t-\tchr\t100\t30\t40\t10\t10\t255\ttp:A:I\tcg:Z:10M\n'
                   'a\t20\t0\t20\t+\torg\t50\t0\t20\t20\t20\t0\ttp:A:S\tcg:Z:20M\n'
                   'b\t20\t5\t15\t+\torg\t50\t20\t30\t10\t10\t0\ttp:A:P\tcg:Z:10M\n')
    assert load_alignments(paf, db, {'chr': 100, 'org': 50})['alignment_rows'] == 4
    summary = summarize_mapping(db, {'chr': 100, 'org': 50}, {'org'}, tmp_path)
    assert summary['read_count'] == 3 and summary['total_bases'] == 60
    assert summary['any_alignment_span_bp'] == 30
    assert summary['primary_span_bp'] == 28
    assert summary['primary_mapq20_span_bp'] == 12
    assert summary['organelle_primary_reads'] == 1
    assert summary['primary_organelle_span_bp'] == 10
    assert summary['primary_other_reference_span_bp'] == 18
    assert summary['primary_compartment_overlap_bp'] == 0
    refs = {r['contig']: r for r in csv.DictReader((tmp_path/'reference_mapping.tsv').open(), delimiter='\t')}
    assert refs['chr']['primary_aligned_target_bases'] == '20'
    assert refs['chr']['primary_reference_union_bp'] == '20'
    assert refs['chr']['mapq20_aligned_target_bases'] == '10'
    assert refs['org']['primary_aligned_target_bases'] == '10'
    reads = list(csv.DictReader((tmp_path/'read_mapping.tsv').open(), delimiter='\t'))
    assert reads[-1]['read_id'] == 'c' and reads[-1]['primary_rows'] == '0'


@pytest.mark.parametrize('bad', [
    'unknown\t20\t0\t10\t+\tchr\t100\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:10M',
    'a\t20\t0\t10\t\tchr\t100\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:10M',
    'a\t20\t0\t10\t+\tchr\t100\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:9M',
    'a\t20\t0\t10\t+\tchr\t100\t0\t10\t10\t10\t60\ttp:A:P\ttp:A:S\tcg:Z:10M',
])
def test_late_invalid_native_alignment_fails_instead_of_reducing_concordance(tmp_path, bad):
    db = indexed_reads(tmp_path)
    paf = tmp_path/'out.paf'
    paf.write_text('a\t20\t0\t10\t+\tchr\t100\t0\t10\t10\t10\t60\ttp:A:P\tcg:Z:10M\n'+bad+'\n')
    with pytest.raises(ValueError):
        load_alignments(paf, db, {'chr': 100})


def test_successful_empty_paf_has_zero_mapping_not_missing_queries(tmp_path):
    db = indexed_reads(tmp_path)
    paf = tmp_path/'out.paf'
    paf.touch()
    assert load_alignments(paf, db, {'chr': 100})['alignment_rows'] == 0
    summary = summarize_mapping(db, {'chr': 100}, set(), tmp_path)
    assert summary['any_mapped_read_fraction'] == 0
    assert summary['read_count'] == 3


def test_overlapping_primary_compartments_are_explicit_not_double_counted(tmp_path):
    db = indexed_reads(tmp_path)
    paf = tmp_path/'out.paf'
    paf.write_text('a\t20\t0\t12\t+\tchr\t100\t0\t12\t12\t12\t30\ttp:A:P\tcg:Z:12M\n'
                   'a\t20\t8\t18\t+\torg\t50\t0\t10\t10\t10\t30\ttp:A:P\tcg:Z:10M\n')
    load_alignments(paf, db, {'chr': 100, 'org': 50})
    summary = summarize_mapping(db, {'chr': 100, 'org': 50}, {'org'}, tmp_path)
    assert summary['primary_span_bp'] == 18
    assert summary['primary_organelle_span_bp'] == 10
    assert summary['primary_other_reference_span_bp'] == 12
    assert summary['primary_compartment_overlap_bp'] == 4
