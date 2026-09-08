# unitFinder source and real-validation scope audit v1

This audit fixes the evidence boundary for adding unitFinder as a serious
assembly-side comparator.  It does not report a successful unitFinder run or an
accuracy result.

The official Genome Biology study applies unitFinder to three telomere-to-
telomere soybean assemblies (Jack, ZH13 and WM82).  The workflow first predicts
potential centromeric regions de novo from tandem-repeat structure and then uses
potential CentGm sequences in a reference-assisted pass.  Cross-chromosome
grouping is partly manual in the released workflow, so end-to-end timing must
include those steps rather than timing only `bin/unitFinder.py`.

The publisher's supplementary workbook was downloaded unchanged.  Its SHA-256
is recorded in `audit.json`; the first PMC download returned a 1,816-byte HTML
interstitial and is retained outside Git as a failed acquisition.  Table S1 has
60 chromosome rows (20 per variety) and two interval sets.  Every reported
length equals `end - start + 1`, supporting a 1-based inclusive interpretation.
The reported TRS interval totals are 40,512,042 bp (Jack), 43,336,133 bp (ZH13)
and 45,793,097 bp (WM82).  These are published workflow outputs and are suitable
for reproducibility checks, not independent truth.

The highest-value real experiment is therefore staged as follows:

1. build the pinned container and pass the already frozen interface smoke;
2. enroll and checksum the exact ZH13 assembly used in the paper;
3. freeze per-chromosome commands, the manual cross-chromosome steps and all
   interval-conversion rules before inspecting new unitFinder output;
4. measure reproduction of Table S1 separately from biological accuracy;
5. derive independent coordinate evidence from the matching CENH3 ChIP-seq
   runs only after their donor/assembly mapping and peak-calling plan is frozen;
6. retain FISH as qualitative family-level support, not base-pair truth.

The current source audit resolves the ZH13 T2T assembly as Genome Warehouse
`GWHBWDJ00000000.1` (BioProject `PRJCA015269`, BioSample `SAMC1127443`).  The
paper's data-availability statement assigns ZH13 ChIP-seq to
`CRR638211`--`CRR638216`.  Those accessions establish an available route to an
orthogonal assay, but per-run input/control roles, replicate structure and the
exact biological relationship to the assembly donor still require enrollment
before the reads can be treated as an accuracy standard.

This comparator does not estimate read-based copy number and must not be placed
in the copy-number table as though its prediction unit were equivalent to
TandemX `quantify`.
