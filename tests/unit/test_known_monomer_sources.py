from pathlib import Path
import shutil

import pytest

from benchmarks.challenge.schema import read_table
from benchmarks.scripts.curate_known_monomers import SOURCE_UNITS, curate, validate_record


RECORDS = Path(__file__).resolve().parents[2]/'paper/evidence/known_repeat_sources/records'


def test_real_source_bank_preserves_materials_and_excludes_multiunit_clone(tmp_path):
    result = curate(RECORDS, tmp_path/'curated')
    assert result['complete'] and result['source_records'] == 4 and result['included_monomers'] == 3
    rows = read_table(tmp_path/'curated/source_units.tsv')
    assert [(r['material'], r['monomer_length']) for r in rows] == [
        ('Seneca 60', '156'), ('IR-BB21', 'NA'), ('c.v. Donetsky A.', '118'), ('Cigalon variety', '358')]
    assert rows[1]['included_monomer'] == 'False'
    assert 'CentO_RCS2_AF058902_1' not in (tmp_path/'curated/known_monomers.fa').read_text()


@pytest.mark.parametrize('fault', ['base', 'version', 'material', 'feature'])
def test_changed_source_content_fails_independent_metadata_and_sequence_checks(tmp_path, fault):
    expected = SOURCE_UNITS[2]
    embl = tmp_path/'record.embl'
    text = (RECORDS/'X16095.1.embl').read_text()
    if fault == 'base': text = text.replace('cgaaactcgc', 'agaaactcgc')
    if fault == 'version': text = text.replace('SV 1;', 'SV 2;')
    if fault == 'material': text = text.replace('c.v. Donetsky A.', 'unresolved')
    if fault == 'feature': text = text.replace('misc_feature    1..118', 'misc_feature    1..117')
    embl.write_text(text)
    with pytest.raises(ValueError):
        validate_record(embl, RECORDS/'X16095.1.fasta', expected)


def test_retrieval_hash_failure_produces_no_completed_bank(tmp_path):
    source = tmp_path/'records'
    shutil.copytree(RECORDS, source)
    (source/'AF078922.1.fasta').write_text('changed')
    with pytest.raises(ValueError, match='retrieval'):
        curate(source, tmp_path/'curated')
    assert not (tmp_path/'curated').exists()
