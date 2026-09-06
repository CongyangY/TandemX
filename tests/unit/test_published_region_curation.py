from zipfile import ZipFile

import pytest

from benchmarks.scripts.curate_mo17_regions import read_cells, extract_regions, chromosome_map, curate


def test_literal_workbook_extraction_and_explicit_sheet_formula_boundary(tmp_path):
    path = tmp_path/'source.xlsx'
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    with ZipFile(path, 'w') as z:
        z.writestr('xl/workbook.xml', f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="6" r:id="r1"/><sheet name="14" r:id="r2"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels', '<Relationships><Relationship Id="r1" Target="worksheets/sheet1.xml"/><Relationship Id="r2" Target="/xl/worksheets/sheet2.xml"/></Relationships>')
        z.writestr('xl/sharedStrings.xml', f'<sst xmlns="{ns}"><si><t>chr1</t></si></sst>')
        z.writestr('xl/worksheets/sheet1.xml', f'<worksheet xmlns="{ns}"><sheetData><row><c r="B4" t="s"><v>0</v></c><c r="C4"><v>100</v></c><c r="A1" t="inlineStr"><is><t>Title</t></is></c><c r="Z99"/></row></sheetData></worksheet>')
        z.writestr('xl/worksheets/sheet2.xml', f'<worksheet xmlns="{ns}"><sheetData><row><c r="A1"><f>2+2</f><v>4</v></c></row></sheetData></worksheet>')
    assert read_cells(path, {'6'}) == {'6': {'B4': 'chr1', 'C4': '100', 'A1': 'Title'}}
    with pytest.raises(ValueError, match='Formula'):
        read_cells(path)


def test_regions_keep_source_cells_and_mixed_array_composition():
    sheets = {n: {} for n in ('6','7','8','9','10','11','12')}
    sheets['6'] = {'A4': '1', 'B4': 'chr1', 'C4': '10', 'D4': '110', 'E4': '100', 'F4': '18.84'}
    sheets['11'] = {'A4': 'chr1', 'B4': '200', 'C4': '400', 'D4': '200', 'E4': '52.34'}
    sheets['12'] = {'A3': 'chr1S', 'B3': '0', 'C3': '10', 'D3': '10', 'E3': '99.55'}
    rows = extract_regions(sheets, {'chr1': ('CM001.1', 1000)})
    assert [r['region_type'] for r in rows] == ['satellite_region','ChIP_defined_centromere','telomere_region']
    assert rows[0]['source_row'] == 4 and rows[0]['coordinate_cells'] == 'C4:E4'
    assert rows[0]['component_percent'] == 18.84
    assert rows[1]['annotated_component'] == 'CentC' and rows[1]['component_cell'] == 'E4'
    assert rows[2]['chromosome_label'] == 'chr1S'
    assert rows[0]['start_raw'] == 10 and 'inferred' in rows[0]['coordinate_note']
    sheets['6']['E4'] = '101'
    with pytest.raises(ValueError, match='end-minus-start'):
        extract_regions(sheets, {'chr1': ('CM001.1', 1000)})
    sheets['6']['E4'] = '100'
    with pytest.raises(ValueError, match='bounds'):
        extract_regions(sheets, {'chr1': ('CM001.1', 50)})


def test_source_checksum_and_wrong_reference_fail_before_writing(tmp_path):
    source = tmp_path/'source.xlsx'
    source.write_bytes(b'not the inspected source')
    with pytest.raises(ValueError, match='checksum'):
        curate(source, tmp_path/'missing', tmp_path/'out')
    assert not (tmp_path/'out').exists()
    report = tmp_path/'report.txt'
    report.write_text('chr1\tassembled-molecule\t1\tChromosome\tCM001.1\t=\tna\tna\t1000\n')
    with pytest.raises(ValueError, match='expected Mo17'):
        chromosome_map(report)
