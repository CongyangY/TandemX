# Ey15-2 original-read controlled edits, 2026-09-17

The input is the frozen Ey15-2 (9994) native-source package in
`benchmarks/controlled_collapse/ey15_native_pair_audit_20260917/`.
Its selected Chr1 interval is 3,244 bp inside a 9,244-bp context with
3,000-bp natural flanks. Seven distinct original HiFi molecules span the
interval with at least 1,000 bp of natural flank on each side; the source
archive records their whole-assembly alignments. The same published sample
is supported, but exact DNA extraction is not separately established.

`build_native_bp_collapse.py` produced nine development-only assemblies:
an intact control; 25%, 50%, 75%, and 100% terminal deletions; 25%, 50%, and
75% internal deletions; and a 25% left-boundary deletion. Exact deleted
lengths for the 3,244-bp interval are 0, 811, 1,622, 2,433, and 3,244 bp.
The source reads were never edited. `receipt.json` lists every source and
edited coordinate, sequence hash, FASTA hash, and the original FASTQ hash.

`verify_native_bp_collapse.py` independently reconstructed all 9/9 edited
FASTA sequences from the frozen source sequence and ledger. A clean
generator replay was byte-identical for the 10 generated files (nine FASTA
files and the receipt); the archived directory additionally contains the
independent `verification.json`. Two focused tests passed, including a
corrupted-FASTA rejection test.

The only exact truth is the **injected bp deletion in the assembly**. The
420-bp periodicity is operational, the 3,244-bp interval is not an integer
multiple of it, and monomer boundaries and physical copy count are unknown.
The frozen M2 route has a 300-bp monomer cap and is ineligible here. This
lineage is part of prior development evidence, not a final held-out donor.

The read-span score protocol was frozen in `score_protocol.json` before
running read-span scoring. Its 1,024-bp flank and 5% threshold are inherited
from the earlier Col-CEN technical pilot; the result will be reported as a
source-guided development diagnostic.
