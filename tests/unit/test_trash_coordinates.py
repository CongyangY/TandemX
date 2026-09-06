import csv

import pytest

from benchmarks.scripts.audit_trash_coordinates import audit


@pytest.mark.parametrize('tool,shift', [('trash', 1), ('trash2', 0)])
def test_audit_recovers_deliberate_offset_on_both_strands(tmp_path, tool, shift):
    fasta = tmp_path/'control.fa'
    fasta.write_text('>chr1\nTACGATTCGGATCCGTTAGCA\n')
    name, seqkey, idkey = (('all.repeats.from.assembly.fa.csv', 'seq', 'seq.name') if tool == 'trash'
                          else ('assembly.fa_repeats_with_seq.csv', 'sequence', 'seqID'))
    # True inclusive positions 3..8 and 12..17, second sequence reverse-complemented.
    rows = [{idkey:'chr1', 'start':3+shift, 'end':8+shift, 'width':6, 'strand':'+', seqkey:'CGATTC'},
            {idkey:'chr1', 'start':12+shift, 'end':17+shift, 'width':6, 'strand':'-', seqkey:'AACGGA'}]
    with (tmp_path/name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    result = audit(fasta, tmp_path, tool, tmp_path/'audit')
    row = next(r for r in result['rows'] if r['offset_from_native_1based'] == -shift)
    assert row['strand_adjusted_matches'] == row['monomer_rows'] == 2
    assert row['forward_matches'] == 1 and row['width_discrepancies'] == 0
    assert all(r['strand_adjusted_matches'] == 0 for r in result['rows'] if r is not row)


def test_small_control_audit_rejects_large_input_before_loading(tmp_path):
    path = tmp_path/'large.fa'
    with path.open('wb') as handle:
        handle.truncate(1_100_001)
    with pytest.raises(ValueError, match='restricted'):
        audit(path, tmp_path, 'trash', tmp_path/'out')
