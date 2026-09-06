"""Small-input coordinate diagnosis; no automatic correction or accuracy ranking."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.schema import digest_file, write_table


def audit(fasta: Path, native: Path, tool: str, outdir: Path) -> dict:
    if fasta.stat().st_size > 1_100_000:
        raise ValueError('This adapter audit is restricted to <=1 Mb controls')
    references = {k: v.upper() for k, v in read_fasta(fasta).items()}
    if sum(map(len, references.values())) > 1_000_000:
        raise ValueError('This adapter audit is restricted to <=1 Mb controls')
    if tool == 'trash':
        name, sequence_key, id_key = 'all.repeats.from.assembly.fa.csv', 'seq', 'seq.name'
    elif tool == 'trash2':
        name, sequence_key, id_key = 'assembly.fa_repeats_with_seq.csv', 'sequence', 'seqID'
    else:
        raise ValueError('Require trash or trash2')
    rows = [dict(offset_from_native_1based=offset, monomer_rows=0, width_discrepancies=0,
                 forward_matches=0, strand_adjusted_matches=0, outside_reference=0)
            for offset in range(-2, 3)]
    with (native/name).open(newline='') as handle:
        reader = csv.DictReader(handle)
        if not {sequence_key, id_key, 'start', 'end', 'width', 'strand'} <= set(reader.fieldnames or []):
            raise ValueError('Unexpected native monomer schema')
        for record in reader:
            sequence = record[sequence_key].upper()
            start, end, width = (int(record[k]) for k in ('start', 'end', 'width'))
            if (record[id_key] not in references or record['strand'] not in {'+', '-'}
                    or start < 1 or end < start or not sequence or set(sequence)-set('ACGTN')):
                raise ValueError('Invalid native monomer record')
            for row in rows:
                row['monomer_rows'] += 1
                row['width_discrepancies'] += int(width != end-start+1 or width != len(sequence))
                a, b = start-1+row['offset_from_native_1based'], end+row['offset_from_native_1based']
                reference = references[record[id_key]]
                if not 0 <= a < b <= len(reference):
                    row['outside_reference'] += 1
                    continue
                fragment = reference[a:b]
                row['forward_matches'] += int(sequence == fragment)
                if record['strand'] == '-':
                    fragment = fragment.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
                row['strand_adjusted_matches'] += int(sequence == fragment)
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir/'coordinate_offsets.tsv', rows, list(rows[0]))
    result = dict(complete=True, tool=tool, input_sha256=digest_file(fasta),
                  native_file=name, native_sha256=digest_file(native/name),
                  script_sha256=digest_file(Path(__file__)), table_sha256=digest_file(outdir/'coordinate_offsets.tsv'),
                  rows=rows, scope='small_control_coordinate_diagnosis;no_automatic_correction_or_biological_accuracy')
    (outdir/'audit_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--tool', choices=('trash', 'trash2'), required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    audit(args.fasta, args.native, args.tool, args.outdir)
