# Rice v2 bounded original-read source gate (2026-09-17)

This is a separate **development** screen after the frozen rice v1 rule returned
no eligible interval. The v1 archive and its negative outcome remain unchanged.
The v2 candidate protocol and code were committed as `07c5a98`; the selected
context and all 235,500 assembly-only candidate scores were committed as
`2020a72` **before any original read was aligned to the selected context**.

## Assembly-only selection

The source is the verified AGIS1.0 assembly SHA-256
`cbe51d6615380bb441b7cf531ebc0ae58f498bc7982d2db07c549376a7638270`.
The frozen v1 top-1,000 assembly-only candidate table was ranked by its
within-context flank 31-mer uniqueness, then shift identity; its top 100
windows were used as fixed centers. Five array lengths (1–5 kb) and every
integer period from 30 to 500 bp produced **235,500** evaluations. Entropy and
31-mer flank uniqueness were computed once per 500 distinct array contexts,
then period-shift identity was evaluated for each period. All scores, including
ineligible ones, are in `all_candidate_scores.tsv.gz`; `region_summary.tsv`
has the 500 context metrics.

Of the 235,500 scores, 920 passed the predeclared period-shift identity 0.85
and sequence entropy 1.6 filters. The deterministic ranking prioritized flank
31-mer uniqueness and then period-shift identity. It selected CP132242.1
`[12995500,13000500)` (0-based half-open), a 5,000-bp array with a 155-bp
operational period, shift identity 0.874097 and flank uniqueness 0.717172.
The selected 11,000-bp assembly context has 3,000 bp of natural sequence on
each side and SHA-256
`48d27907c4d51d9f2017db28dee38eebabdca84808674c5d16165cd2832fced8`.
The within-context uniqueness is a ranking statistic, **not** proof of
genome-wide uniqueness. It is below v1's hard 0.8 threshold, so mapping
ambiguity remains material.

## Fixed original-read gate

The only read input was the previously sampled, currently validated
SRR25241090 `sample_003` gzip: 1,071,470,526 bytes, 62,345 original complete
FASTQ records, SHA-256
`88e05ea73487dffa014ca3f8bf87881eb58071635fcde65a75349264d20a821f`.
The complete checksum/gzip/FASTQ readback is in the v1 source archive. No
full-run scan, new download, or T7 write was performed.

Minimap2 2.31-r1302 with the frozen `-x map-hifi -c --secondary=yes -N 20 -t 2`
settings reported 9,068 context alignments, including 3,041 primary alignment
segments across 1,265 distinct read IDs. Exactly **one** distinct original
read met all fixed criteria: primary alignment, at least 99% matching bases
per alignment column, and target coverage extending at least 1,000 bp into
both natural flanks. Its ID and coordinates are in `qualifying_spanners.tsv`.
The required denominator was at least **three distinct reads**. The source
therefore failed its technical spanner gate. The full PAF, minimap2 stderr,
and machine-readable result are preserved. The mapper reported 46.252 s wall,
90.499 s CPU and 2.098 GB peak RSS for this **one 11-kb-context versus
62,345-read mapping**; those numbers are not TandemX genome-scale resource
measurements.

The predeclared full-reference anchoring step was **not run** because the
spanner gate failed. No alternative interval was selected after read mapping,
and no controlled-collapse assembly edit was generated. A second independent
calculation from the frozen candidate score archive and raw PAF reproduced
the 235,500/920 candidate counts, selected coordinates and one qualifying
read.

The assembly BioSample SAMN36344332 and read BioSample SAMN36368305 have the
same Nipponbare/AGIS-1.0 label but different archive records and collection
dates. Exact individual and DNA extraction are unverified. The result is a
technical development failure for this selected interval; it has no natural
copy-number truth or final-heldout status.
