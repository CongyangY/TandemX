"""Extract published Mo17 regions with source cells and explicit truth limits.

Read-only XLSX extraction uses standard-library XML, with no Excel dependency.
This is not a monomer library or a base-level satellite truth set.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path, PurePosixPath
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from benchmarks.challenge.schema import digest_file, write_table

XLSX_SHA256 = '863da22bbc8fa3170035f966f6ced723e418e81541fda96de0e5cc9355569e1d'
SOURCE_URL = ('https://media.springernature.com/original/springer-static/esm/'
              'art%3A10.1038%2Fs41588-023-01419-6/MediaObjects/41588_2023_1419_MOESM4_ESM.xlsx')
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def read_cells(path: Path, sheet_names: set[str] | None = None) -> dict[str, dict[str, str]]:
    """Extract literal worksheet cells; refuse formulas in source input."""
    with ZipFile(path) as archive:
        if sum(e.file_size for e in archive.infolist()) > 50_000_000:
            raise ValueError('Workbook exceeds bounded extraction size')
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            shared = [''.join(s.itertext()) for s in
                      ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', NS)]
        relationships = {r.attrib['Id']: r.attrib['Target'] for r in
                         ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
        sheets = {}
        for sheet in ET.fromstring(archive.read('xl/workbook.xml')).findall('s:sheets/s:sheet', NS):
            if sheet_names is not None and sheet.attrib['name'] not in sheet_names:
                continue
            rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target = relationships[rid]
            source = str(PurePosixPath('xl') / target) if not target.startswith('/') else target.lstrip('/')
            cells = {}
            for cell in ET.fromstring(archive.read(source)).findall('.//s:sheetData/s:row/s:c', NS):
                if cell.find('s:f', NS) is not None:
                    raise ValueError('Formula source requires an explicit cached-value audit')
                value = cell.find('s:v', NS)
                if cell.attrib.get('t') == 'inlineStr':
                    text = ''.join(cell.find('s:is', NS).itertext())
                elif value is None:
                    continue
                elif cell.attrib.get('t') == 's':
                    text = shared[int(value.text)]
                else:
                    text = value.text
                cells[cell.attrib['r']] = text
            sheets[sheet.attrib['name']] = cells
    return sheets


def chromosome_map(report: Path) -> dict[str, tuple[str, int]]:
    mapping = {}
    for line in report.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) < 9 or fields[1] != 'assembled-molecule' or not fields[2].isdigit():
            raise ValueError('Require numbered assembled chromosomes in the reference report')
        chromosome = 'chr'+fields[2]
        if chromosome in mapping:
            raise ValueError('Duplicate chromosome in NCBI report')
        mapping[chromosome] = (fields[4], int(fields[8]))
    if len(mapping) != 10 or sum(length for _, length in mapping.values()) != 2_178_604_320:
        raise ValueError('Reference report is not the expected Mo17 chromosome set')
    return mapping


def extract_regions(sheets: dict[str, dict[str, str]], mapping: dict[str, tuple[str, int]]) -> list[dict]:
    """Preserve raw coordinates; identify half-open inference and component purity."""
    result = []
    families = {'6': 'TR-1', '7': 'knob180', '8': 'CentC'}
    for sheet in ('6', '7', '8', '9', '10', '11', '12'):
        cells = sheets[sheet]
        # Data are recognized by a chromosome-labelled cell, not arbitrary blanks.
        chr_column = 'A' if sheet in ('11', '12') else 'B'
        data_rows = sorted(int(key[1:]) for key, value in cells.items()
                           if key.startswith(chr_column) and re.fullmatch(r'chr\d+(?:[SL])?', value))
        for row in data_rows:
            get = lambda col: cells[f'{col}{row}']
            chromosome_label = get(chr_column)
            chromosome = re.sub('[SL]$', '', chromosome_label)
            accession, length = mapping[chromosome]
            offset = 0 if sheet in ('11', '12') else 1
            coords = [get(chr(ord('B')+offset+i)) for i in range(3)]
            if any(not re.fullmatch(r'\d+', v) for v in coords):
                raise ValueError('Published coordinates/size must be literal nonnegative integers')
            start, end, size = map(int, coords)
            if not 0 <= start < end <= length or size != end-start:
                raise ValueError('Published region violates chromosome bounds or end-minus-start size')
            if sheet in families:
                family, component_column, kind = families[sheet], 'F', 'satellite_region'
            elif sheet == '9':
                family = get('A').removesuffix(' array')
                component_column = {'sat112': 'F', 'sat261': 'G', 'sat268': 'H', 'Cent4': 'I', 'tRNAsat': 'J'}[family]
                kind = 'satellite_region'
            elif sheet == '10':
                family, kind = get('A'), 'rDNA_region'
                component_column = {'5S rDNA': 'F', '45S rDNA': 'G'}[family]
            elif sheet == '11':
                family, component_column, kind = 'CENH3_centromere', 'E', 'ChIP_defined_centromere'
            else:
                family, component_column, kind = 'telomeric_repeat', 'E', 'telomere_region'
            component = 'CentC' if sheet == '11' else family
            percentage = float(get(component_column))
            if not 0 <= percentage <= 100:
                raise ValueError('Invalid published sequence-component percentage')
            result.append(dict(region_id=f'Mo17_S{sheet}_R{row}', region_type=kind,
                               published_label=family, chromosome_label=chromosome_label,
                               reference_contig=accession, reference_length_bp=length,
                               start_raw=start, end_raw=end, published_size_bp=size,
                               coordinate_note='end_minus_start_consistent_origin_inferred_0_based_not_explicit',
                               annotated_component=component, component_percent=percentage,
                               source_sheet=sheet, source_row=row, coordinate_cells=f'{chr(ord("B")+offset)}{row}:{chr(ord("D")+offset)}{row}',
                               component_cell=f'{component_column}{row}', source_url=SOURCE_URL,
                               source_sha256=XLSX_SHA256,
                               warning='published_region_not_base_level_repeat_truth;not_independent_of_assembly;component_percent_rounded'))
    return result


def curate(source: Path, report: Path, outdir: Path) -> None:
    if digest_file(source) != XLSX_SHA256:
        raise ValueError('Published source checksum differs from the inspected workbook')
    mapping = chromosome_map(report)
    sheets = read_cells(source, {'4', '6', '7', '8', '9', '10', '11', '12'})
    regions = extract_regions(sheets, mapping)
    counts = Counter(row['region_type'] for row in regions)
    if counts != {'satellite_region': 64, 'rDNA_region': 2, 'ChIP_defined_centromere': 10, 'telomere_region': 20}:
        raise ValueError('Published region counts differ from inspected source tables')
    if int(sheets['4']['B28']) != sum(length for _, length in mapping.values()):
        raise ValueError('Published total genome size differs from NCBI chromosome total')
    telomeres = [r for r in regions if r['region_type'] == 'telomere_region']
    for r in telomeres:
        if r['chromosome_label'].endswith('S') and r['start_raw'] != 0:
            raise ValueError('Published chromosome start disagrees with zero-based inference')
        if r['chromosome_label'].endswith('L') and r['end_raw'] != r['reference_length_bp']:
            raise ValueError('Published chromosome end differs from reference length')
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir/'published_regions.tsv', regions, list(regions[0]))
    receipt = dict(complete=True, source_sha256=XLSX_SHA256, source_url=SOURCE_URL,
                   report_sha256=digest_file(report), script_sha256=digest_file(Path(__file__)),
                   counts=dict(counts), source_total_bases=2_178_604_320,
                   all_ten_chromosome_end_coordinates_match_reference=True,
                   source_start_coordinate_zero_in_table12=True,
                   coordinate_origin='inferred from size arithmetic and telomeres; retain raw values',
                   excluded='Table13 has mixed size arithmetic; not normalized',
                   warning='Numeric coordinate compatibility does not establish assembly byte identity or donor identity',
                   output_sha256=digest_file(outdir/'published_regions.tsv'))
    (outdir/'curation_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-xlsx', type=Path, required=True)
    parser.add_argument('--ncbi-report', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    curate(args.source_xlsx, args.ncbi_report, args.outdir)
