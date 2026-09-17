# Ey15-2 original HiFi interval source audit (2026-09-17)

This is the **second independent development assembly lineage** for bounded
controlled edits, separate from Col-CEN. It is not final held-out evidence.
The [Ey15-2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757041/) reports
HiFi and CLR assemblies of sample 9994; the original HiFi run is ERR8666125.
The newer assembly shares HiFi evidence with these reads, so it is a
high-quality reference proxy, not independent physical copy-number truth.

## Pre-read interval selection

The prior archived `locate_new/arrays.bed` and `discover/families.tsv` were
screened before mapping any reads in this audit. Sixteen located intervals were
3,000–5,000 bp long and had a family period in the frozen catalogue. For each,
the direct fraction of matching bases at that period shift was measured on the
assembly sequence. `candidate_periodicity_screen.tsv` retains all sixteen,
including low-purity candidates. The top interval was TXF000260:

| Input | Value |
| --- | --- |
| Source assembly | Ey15-2 9994 HiFi Hifiasm, five chromosomes, 136,162,473 bases |
| Source FASTA SHA-256, current full readback | `48d67e8d26c5d9ebc1fbd772e877e2f957437b5a3b27bbdde4a1d3c77923208c` |
| Located array, 0-based half-open | Chr1 `[12829234,12832478)`; 3,244 bp |
| Operational family period | 420 bp, direct shift identity 0.9946883853 |
| Natural-flank context | Chr1 `[12826234,12835478)`; 9,244 bp, array `[3000,6244)` |
| Archived context FASTA SHA-256 | `35b10d849c6a97bb6a89640f2e54832c5b92bc7eb966e06f2153a7e9452a5ace` |

Current full FASTA readback also confirmed five unique nonempty Chr1–Chr5
records, 136,162,473 A/C/G/T/N bases, and the source receipt's contig lengths.

The 420-bp shift is operational periodicity. The 3,244-bp interval is not an
integer number of 420-bp units; individual monomer boundaries and physical
copy count are unresolved. Exact **injected bp changes** can be recorded in a
controlled edit ledger. A copy-count or naturally missing-base truth label
cannot be inferred from this interval. The frozen M2 route's 300-bp period cap
makes this interval ineligible for that route without changing its rules.

## Original-read support

The existing seeded `sample_003` is 50,267 whole original HiFi reads selected
from ERR8666125 at 0.06 inclusion probability. Its 1,125,515,371-byte gzip
passed current complete SHA-256 readback
(`1f05335279f80e604d0f451a72741d45fb53e31d35e139023e0afbe528b26baf`),
gzip trailer and FASTQ structure checks. Seven records with different PacBio
ZMW identifiers were extracted **without changing their headers, sequences or
qualities**. `raw_read_extract_receipt.json` records each exact four-line
record SHA-256. Their compressed selected FASTQ SHA-256 is
`e38b4bb7726a9975992137f681c56ccd1adc8ee15dde1401ad447dcc27a736d0`.
These are seven molecules from one movie and one material, not biological
replicates. The full 18.6-Gb source FASTQ was not reread in this audit.

Minimap2 2.31-r1302 `-x map-hifi -c --secondary=yes -N 10` was run first on
the context against the bounded 50,267-read subset, then on the exact seven
selected records against the **entire five-chromosome 136-Mb assembly**. All
seven have a main context alignment spanning the array with at least 1 kb of
natural sequence on both sides and >=0.99 nmatch/alignment-column identity.
All seven have one reported whole-assembly alignment, primary on Chr1 across
the selected array (identity 0.994559–0.999133, MAPQ 60). Three extra short
local-context alignments arise from repeated sequence and remain in the PAF.
`alignment_summary.tsv` preserves per-read coordinates, spans and identity;
both PAF files preserve the actual aligner output. An aligner's sole reported
hit on one assembly is strong locus support, subject to mapping heuristics and
the representation of that assembly; it does not prove an absent alternative
haplotype or absolute molecule origin.

`source_eligibility_manifest.json` gives all input and output hashes, source
accessions, read-pairing level, coordinates, and eligibility decisions. The
archive contains small original-read and assembly extracts only. No T7 write,
new download, or 11-GB-class read scan was performed. The earlier Col-CEN
archive was left unchanged.

Run the artifact guard without T7 access:

```bash
conda run -n tandemx-dev pytest -q tests/unit/test_ey15_native_pair_audit.py
```
