# Mo17 complete-input QC and nested sampling

One included CCS-labelled run, SRR15447419, from PRJNA751841. Full source size and
MD5 verified before QC: 407,670 reads / 5,624,644,958 bases. N50 13,968 bp; median
13,455 bp; GC 45.9219%; no N bases or duplicate archive IDs. Mean reported base
error probability .00288144 (Phred transform 25.4039) is not empirical accuracy.
All read-quality bins start at >=20; no additional filtering was applied.

Seed 6101, namespace SRR15447419, independent inclusion probabilities:

| Fraction | Selected reads | Observed bases | Nominal bases / reference size |
| --- | ---: | ---: | ---: |
| .002 | 840 | 11,680,888 | .00536164 |
| .02 | 8,084 | 111,505,681 | .05118216 |
| .2 | 81,775 | 1,128,793,699 | .51812699 |

These are nested samples of one technical batch, not independent materials or
the full published Mo17 HiFi dataset. Whole reads were selected after full-file
QC; archive order did not determine membership. Nuclear depth/contaminants and
library bias are not resolved by total bases or these distributions.

The paper-linked NCBI reference GCA_022117705.1 was downloaded separately and
MD5/full-FASTA checked: 10 chromosomes, 2,178,604,320 bp, no N or other ambiguous
bases. Archive `reference/` contains source reports, sequence-to-chromosome mapping,
checksums and QC; raw FASTA is at T7. Reference and reads share the source project
and Mo17 cultivar context, but exact donor/extraction identity is not established.
The NCBI release is not assumed byte-identical to the alternate MaizeGDB file.
CyVerse presented a browser verification page and MaizeGDB returned HTTP403 to
the local downloader; the paper's independently public NCBI record was used.

`figures_checked/input_qc.pdf` and SVG are the visually checked four-panel
distribution/data-volume figure. A: read-length CDF; B: fraction per 1% GC bin;
C: CDF of lower bounds of 5-Phred mean-read-quality bins; D: observed/expected
base amounts on a log axis. All probabilities are relative to each included
dataset's own read count. SVG has 70 editable text elements and zero images.
`panel_source.tsv` and figure/source hashes support every plotted value.

Reproduce with `benchmarks/scripts/plot_complete_qc.py` and the archived qc/sampling
folders. Full data acquisition/sampling commands are in `docs/complete_data_qc.md`.
The raw FASTQ, exact selected-ID audits and selected FASTQs stay under
`/Volumes/T7/Codex/TandemX/data/`; they are reproducible from source hashes and seed.
No de novo accuracy, abundance calibration, biological replicate count or final
resource comparison is established by this QC figure.
