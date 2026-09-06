"""Validate a small source-backed monomer bank for post hoc biological checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.schema import digest_file, write_table


SOURCE_UNITS = [
    dict(known_id='CentC_AF078922_1', accession='AF078922.1', source_length=156, start1=1, end1=156,
         sequence_md5='81f1c7ff0d4424fafd5cab475cb9a438', species='Zea mays', material='Seneca 60',
         role='source_repeat_unit', evidence='https://doi.org/10.1073/pnas.95.22.13073',
         boundary='one_historical_CentC_variant_not_a_Mo17_truth_catalogue'),
    dict(known_id='CentO_RCS2_AF058902_1', accession='AF058902.1', source_length=639, start1=None, end1=None,
         sequence_md5='6a0b60272c10e2daf0ac0286418e2d9a', species='Oryza sativa Indica Group', material='IR-BB21',
         role='multiunit_clone_excluded_from_monomer_scoring', evidence='https://doi.org/10.1073/pnas.95.14.8135',
         boundary='639bp_clone_not_a_155bp_monomer;indica_not_Nipponbare'),
    dict(known_id='HvT01_X16095_1_1_118', accession='X16095.1', source_length=208, start1=1, end1=118,
         sequence_md5='3a1d13f9b033234d3889a1d7837c8d9d', species='Hordeum vulgare', material='c.v. Donetsky A.',
         role='explicitly_annotated_monomer_segment', evidence='https://doi.org/10.1111/pbi.13816',
         boundary='source_monomer_feature_1_to_118_not_full208bp_record;not_Morex_donor'),
    dict(known_id='Rice358_X55642_1', accession='X55642.1', source_length=358, start1=1, end1=358,
         sequence_md5='1b271c520a06efd80d99806cdf18a89a', species='Oryza sativa', material='Cigalon variety',
         role='deposited_repeat_unit', evidence='https://www.ebi.ac.uk/ena/browser/view/X55642.1',
         boundary='database_record_lists_manuscript_unpublished;not_Nipponbare_donor;exact_native_unit_boundaries_unvalidated'),
]


def validate_record(embl: Path, fasta: Path, expected: dict) -> str:
    text = embl.read_text()
    identifier = re.search(r'^ID\s+(\w+); SV (\d+);.*; (\d+) BP\.$', text, re.M)
    if (identifier is None or f'{identifier[1]}.{identifier[2]}' != expected['accession']
            or int(identifier[3]) != expected['source_length']):
        raise ValueError('EMBL accession/version/length differs from curated source')
    md5 = re.search(r'^DR\s+MD5; ([0-9a-f]{32})\.', text, re.M)
    if md5 is None or md5[1] != expected['sequence_md5']:
        raise ValueError('EMBL sequence checksum differs from pinned record')
    if expected['material'] not in text or expected['species'] not in text:
        raise ValueError('Expected source species/material missing from record')
    if '\nSQ ' not in text or not text.rstrip().endswith('//'):
        raise ValueError('Missing EMBL sequence or terminal marker')
    block = text.split('\nSQ ', 1)[1].split('\n', 1)[1].rsplit('//', 1)[0]
    sequence = re.sub(r'[\s\d]', '', block).upper()
    if len(sequence) != expected['source_length'] or set(sequence)-set('ACGTN'):
        raise ValueError('EMBL sequence length/alphabet invalid')
    records = read_fasta(fasta)
    if len(records) != 1 or expected['accession'] not in next(iter(records)) or next(iter(records.values())) != sequence:
        raise ValueError('EMBL and FASTA sequence/version disagree')
    if hashlib.md5(sequence.encode()).hexdigest() != expected['sequence_md5']:
        raise ValueError('Sequence fails official MD5')
    if expected['role'] == 'explicitly_annotated_monomer_segment':
        feature = f"misc_feature    {expected['start1']}..{expected['end1']}"
        if feature not in text or '/note="repeat monomer"' not in text:
            raise ValueError('Explicit monomer feature missing')
    if expected['start1'] is not None and not 1 <= expected['start1'] <= expected['end1'] <= len(sequence):
        raise ValueError('Monomer extraction coordinates invalid')
    return sequence


def curate(records_dir: Path, outdir: Path) -> dict:
    retrieval = json.loads((records_dir/'retrieval_receipt.json').read_text())
    sources = {r.get('file'): r for r in retrieval['records'] if 'file' in r}
    rows, monomers = [], []
    for expected in SOURCE_UNITS:
        paths = [records_dir/(expected['accession']+'.'+suffix) for suffix in ('embl', 'fasta')]
        for path in paths:
            if path.name not in sources or digest_file(path) != sources[path.name]['sha256']:
                raise ValueError('Missing or changed original retrieval: '+path.name)
        sequence = validate_record(*paths, expected)
        unit = sequence[expected['start1']-1:expected['end1']] if expected['start1'] is not None else None
        rows.append(dict(**expected, source_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                         embl_sha256=digest_file(paths[0]), fasta_sha256=digest_file(paths[1]),
                         included_monomer=unit is not None, monomer_length=len(unit) if unit else None,
                         monomer_sha256=hashlib.sha256(unit.encode()).hexdigest() if unit else None,
                         warning='posthoc_reference_not_discovery_input_or_complete_family_truth'))
        if unit is not None:
            monomers.append(f">{expected['known_id']} accession={expected['accession']};start1={expected['start1']};end1={expected['end1']}\n{unit}\n")
    outdir.mkdir(parents=True, exist_ok=False)
    (outdir/'known_monomers.fa').write_text(''.join(monomers))
    write_table(outdir/'source_units.tsv', rows, list(rows[0]))
    result = dict(complete=True, source_records=len(rows), included_monomers=len(monomers),
                  script_sha256=digest_file(Path(__file__)),
                  retrieval_receipt_sha256=digest_file(records_dir/'retrieval_receipt.json'),
                  files={name: digest_file(outdir/name) for name in ('known_monomers.fa', 'source_units.tsv')},
                  boundary='three_source_units_not_complete_plant_repeat_truth;source_cultivars_differ_from_test_materials')
    (outdir/'curation_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    curate(args.records_dir, args.outdir)
