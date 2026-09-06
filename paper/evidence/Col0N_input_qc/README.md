# Complete Col-0N input QC

ERR6210723 / PRJEB46164: complete ENA FASTQ size and source MD5 verified, then
the full gzip stream parsed. Observed 933,904 reads / 14,646,601,458 bases,
median length 15,514 bp, N50 15,663 bp, GC 0.3676189867, N=0 and zero duplicate
archive read IDs. SHA256:
`79dd524760b6db4eb9d44db5a4776a25dd4b72e80b24186213e92ea83f264c25`.

Reported quality gives mean error probability 0.0015386005; this is not measured
empirical read accuracy. Archive ID uniqueness is not proof of unique molecules.
The source is Col-0N pooled plants, not a single plant or a biological replicate
of Col-0R. Biological identity, organellar content, coverage bias and reference
uncertainty require downstream evidence. Full-library hash sampling is separate
from this completed file QC. See `archive_manifest.json` for byte-preserved source
receipts and histograms. Raw reads remain under the T7 data root.

Seed6101 nested full-library sampling completed: fractions .0008/.008/.08/.36
gave 752 / 7,557 / 74,600 / 336,613 reads, with 11,765,868 / 118,496,530 /
1,170,032,165 / 5,280,163,009 bases. Exact FASTQ/ID hashes and distributions are
in `subsets/`; raw sequences and ID lists remain on T7. The final subset is about
40 nominal nuclear genome equivalents using 131,559,676 bp. Organellar content
has not been removed, so this is not measured nuclear depth.

Four-panel `figures/input_qc.pdf` and editable SVG show length, GC, reported
quality and observed/expected sampling volume. PNG visually inspected. Every
point/curve is linked to source histograms; these are not independent specimens.
