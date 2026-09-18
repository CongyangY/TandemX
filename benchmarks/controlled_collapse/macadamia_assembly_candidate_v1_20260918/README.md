# Macadamia assembly-only interval freeze (2026-09-18)

This development selection follows the **pre-existing, committed** rule in
`../macadamia_native_pair_audit_20260917/protocol.json` (SHA-256
`2b986f3f0d3ff658d51304b58b63a3a8c13552bb82dac2eb7491de5e4178c415`).
It uses the GigaDB100906 primary assembly and TandemX's earlier `locate_new`
BED plus `discover` family periods. No original HiFi FASTQ, mapping, or read
support was inspected to select the interval.

The rule admits arrays 3,000–5,000 bp, family period 30–300 bp and 3,000 bp
of natural assembly flank on each side. Among admitted arrays it chooses the
greatest direct exact period-shift identity, then longer interval, then
lexicographic contig and ascending start. The catalogue has 28,368 BED rows;
seven pass the length range and two also pass the period range. All seven
length-filter decisions are in `length_filter_trace.tsv`; both sequence
scores are in `all_candidate_scores.tsv`.

The selected 0-based half-open interval is
**ctg.000105F:[118439,121594)**, family `TXF000219`, operational period
**144 bp**. Its direct period-shift identity is **2796/3011 = 0.9285951511**.
`selected_context.fa` contains ctg.000105F:[115439,124594), a 9,155-bp
assembly sequence with 3,000 bp of left and right flank. It contains no
non-ACGT bases. The context sequence SHA-256 (without FASTA header/newline)
is `583a95f41d7132baa76b90f6b32d3619435496dabbee2a3c5e980801ea203b03`.
`selection_receipt.json` records current source SHA-256 values, counts,
all output checksums and the selected locus.

Replay on a verified source copy with:

```bash
conda run -n tandemx-dev python benchmarks/scripts/screen_macadamia_assembly_candidate.py \
  --protocol benchmarks/controlled_collapse/macadamia_native_pair_audit_20260917/protocol.json \
  --output-dir /path/to/new-empty-output-directory
```

This is an **assembly-only development nomination**, not a native-read
eligible interval or biological copy-number truth. The operational period
comes from an earlier TandemX family and may not represent a complete HOR.
No global flank uniqueness, donor identity, molecule support, or physical
under-representation claim has been established by this package. If the
predeclared native-read screen fails, retain that negative result; do not
substitute a different locus based on read outcomes.
