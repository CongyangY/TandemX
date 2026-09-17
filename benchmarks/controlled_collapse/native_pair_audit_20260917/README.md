# Original-read source audit for controlled edits (2026-09-17)

This is a **development input**. The three 13,560-bp contexts in
`colcen_native_contexts.fa` are exact extracts of the internally rechecked
Col-CEN v1.2 assembly (compressed FASTA SHA-256
`b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8`).
Each frozen B1 v3 array occupies context `[5000,8560)`; 5,000 bp of **natural**
assembly sequence was retained on each side. The array subsequences match
`../v3/source_arrays.fa` byte for byte. The three contexts remain one assembly
lineage, not three biological donors.

| Context | Col-CEN chromosome and 0-based context | Frozen array | Original reads with a primary local alignment spanning the array and at least 1 kb on each side |
| --- | --- | --- | ---: |
| C1 | Chr1 `[15580019,15593579)` | `[15585019,15588579)` | 2 |
| C2 | Chr2 `[3960542,3974102)` | `[3965542,3969102)` | 1 |
| C3 | Chr3 `[15055000,15068560)` | `[15060000,15063560)` | 3 |

The six unmodified FASTQ records in `colcen_native_spanners.fastq.gz` were
extracted from the existing seeded ERR6210723 `sample_003` (74,600 reads;
compressed SHA-256
`3460277f5659bf362f2b32050ae0333cd8762ee58db3e8c65dc338e60978b8a9`).
The current readback checked the source SHA-256, all FASTQ records and gzip
trailer, then stored each selected record's SHA-256 in
`raw_read_extract_receipt.json`. The six archive IDs have six different PacBio
ZMW IDs from one movie. They are separate molecules, not biological replicates.

The 8.9x nested subset was screened for at least two exact 41-mers in each
natural flank, unique **within the three selected contexts**. Twenty-one
candidate reads were aligned to all three contexts using minimap2 2.31-r1302
with `-x map-hifi -c --secondary=yes -N 10`. The six retained alignments are
primary, have local-context MAPQ >=20, nmatch/alignment-column identity >=0.99,
and extend >=1 kb into each natural flank. `alignment_summary.tsv` gives each
read's coordinates, identity and best reported alignment to another selected
context; `colcen_native_spanners.paf` preserves all reported alignments for the
six records, including secondary hits. C3's three alignments have 0.99835 to
0.99976 identity; its observed other-context alignments are <=0.931689. This
is **local three-context discrimination**, not a whole-genome uniqueness test.
MAPQ 60 must not be read as whole-genome locus proof here.

The [Col-CEN study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10164409/)
links the assembly to Col-0 ONT/HiFi data and lists HiFi project PRJEB46164.
It does not establish that this specific read library and the assembled
chromosome came from an identical plant/extraction. The source pairing is
therefore limited to study and accession. Prior use of the Col-CEN lineage
keeps these data in development; it cannot serve as independent final held-out
evidence. The exact edited-assembly ledger can establish injected edit truth,
but not native biological missing-copy truth.

`eligibility_manifest.json` also records three distinct assembly lineages with
existing genuine HiFi reads. [Ey15-2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757041/)
is a published same-sample 9994 source pair, but this audit did not recruit
native reads to a chosen interval. The [Macadamia jansenii update](https://www.gigabytejournal.com/articles/24)
reports a higher-contiguity HiFi IPA assembly from the same sample as an
earlier CLR comparison; archival BioSample identifiers differ and the same DNA
extraction is unverified. Its 284-contig, 4.49-Mb-N50 assembly is a reference
proxy, with no native interval selected in this audit. [Mo17 SRR15447419](https://www.ncbi.nlm.nih.gov/sra/SRX11746829)
and [GCA_022117705.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022117705.1/)
share a project/cultivar; their BioSamples differ, so exact donor pairing is
pending. The Col-CEN, Ey15-2 and Mo17 small 11.7-Mb read subsets were copied to internal `tmp/` and
rechecked for source/copy SHA-256, gzip, and FASTQ structure; those subsets do
not establish array support. No new public download or large T7 write was made.

Source and alignment SHA-256 values, exact 0-based coordinates, and per-read
metrics are in `source_and_alignment_receipt.json`. Full raw FASTQs and the
Ey15/Mo17 whole assemblies retain their historical receipts; they were not
fully reread during this audit. C3 has three local-context spanning molecules;
C1/C2 remain below a three-molecule M2 threshold. Native reads are 11.8-19.6
kb long and require a prespecified original-subsequence window to meet a
<=4,096-bp route input cap. Such clipping must be labeled explicitly and must
retain real sequence, quality, and distinct ZMW provenance.

Focused verification:

```bash
conda run -n tandemx-dev pytest -q tests/unit/test_prepare_native_pair_audit.py
shasum -a 256 benchmarks/controlled_collapse/native_pair_audit_20260917/colcen_native_spanners.fastq.gz
```
