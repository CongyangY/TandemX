# Mo17 native-read bounded source screen (2026-09-17)

**Outcome:** the prespecified CentC array has zero qualifying original HiFi
spanners in the existing 81,775-read nested subset. It is **ineligible** for a
three-read controlled-collapse development case under this screen. No other
interval was selected after inspecting the read alignments.

The Mo17 T2T assembly `GCA_022117705.1` and `SRR15447419` share cultivar and
BioProject `PRJNA751841`, but their archive BioSamples differ
(`SAMN20854702` versus `SAMN20604744`). Exact donor/extraction pairing and
physical array copy count remain unresolved. The assembly provides a reference
proxy, never biological missing-base truth. This is development material, not
an independently held-out test case.

## Frozen selection before read mapping

Six hand-selected rows from the published Mo17 satellite-region supplement were
screened on assembly sequence for 3-kb tandem windows, using a 250-bp window
step and a nominal-period shift identity. The exact six extracted source
sequences are in `pre_read_six_candidate_regions.fa`, and the fully reproducible
limited screen is recorded in `region_periodicity_screen.tsv`; it is not an
exhaustive 64-region search.
The highest-scoring window was `Mo17_S8_R7`, a published CentC region on
`CM039151.1` `[98993717, 99023934)` (published coordinate origin inferred by
size arithmetic). The chosen 3,000-bp array is `[98994467, 98997467)`, with
155-bp shift identity 0.934974. The `mo17_native_context.fa` record is the
unmodified assembly `[98991467, 99000467)`, with 3,000 bp of assembly sequence
on each side and array positions `[3000, 6000)` within it. The earlier
`Mo17_S7_R9` in the six-region source FASTA is a pre-read knob180 candidate that
failed the periodicity screen (best shifted identity 0.559532) and was never
mapped to reads. Neither published component percentage nor shift identity
establishes locus uniqueness or an exact monomer boundary.

## Bounded input integrity and read result

The complete compressed assembly on T7 has current exact SHA-256
`2ac74f8a6dc9f3651cd54fcfb4ba644c012a299aa2c6d3ebcd23d756192d15e0`;
its gzip trailer and all 10 FASTA records / 2,178,604,320 ACGT bases were
read and checked. The existing nested `sample_003.fastq.gz` has current SHA-256
`58dfb8812b54f536dcf3ca76a1b41d75511aeab6f0617f58181884acb6f2bb81`,
valid gzip trailer and FASTQ structure, 81,775 reads / 1,128,793,699 bases
(nominal 0.518x assembly bases). The complete 5.35-GB read source was not
read back under the current T7 USB fault; its hash is an archived acquisition
receipt only. No new T7 data were written.

Minimap2 `2.31-r1302`, `-x map-hifi -c --secondary=yes -N 10 -t 4`, mapped the
bounded source to the frozen 9-kb context. `sample_003_context_map_hifi.paf`
contains 1,433 local alignments from 421 distinct reads. A primary alignment
had to reach at least 1 kb into both natural reference flanks and have
whole-alignment exact identity at least 0.95. Six different original reads
spanned geometrically, but their identities were 0.689264–0.900597; zero
qualified. Their exact unedited FASTQ records and per-record SHA-256 hashes
are in `geometric_spanners_original.fastq.gz` and
`read_extraction_receipt.json` as **diagnostic failures**, not supporting
native molecules. `local_alignment_audit.json` reports their left flank,
array and right flank identities; even the best read has 0.909253, 0.901296
and 0.890277 across those parts. The low identity is not confined to the
array. Repeat cross-mapping and donor divergence cannot be distinguished from
this local mapping. Local MAPQ 60 is not genome-wide uniqueness.

Because no read met the frozen identity threshold, full-reference mapping and
controlled edits were not run. The result cannot establish a true physical
array length or classify collapse. Scanning another locus or the full read
source would require a separately frozen selection rule and safe input path.

## Reproduce the bounded analysis

Source paths, receipt references, SHA-256 values, coordinates and eligibility
are in `source_eligibility_manifest_v1.json`. With the documented T7 inputs
available and the `tandemx-dev` environment active:

```bash
python benchmarks/scripts/screen_mo17_native_candidates.py \
  --regions pre_read_six_candidate_regions.fa \
  --output region_periodicity_screen.tsv
minimap2 -x map-hifi -c --secondary=yes -N 10 -t 4 mo17_native_context.fa \
  /Volumes/T7/Codex/TandemX/data/subsets/SRR15447419_seed6101_v1/sample_003.fastq.gz \
  > sample_003_context_map_hifi.paf
python benchmarks/scripts/audit_native_spanners.py \
  --paf sample_003_context_map_hifi.paf --context mo17_native_context.fa \
  --diagnostic-reads geometric_spanners_original.fastq.gz \
  --array-start 3000 --array-end 6000 --flank-bp 1000 \
  --minimum-identity 0.95 --output local_alignment_audit.json
pytest -q tests/unit/test_audit_native_spanners.py \
  tests/unit/test_screen_mo17_native_candidates.py
```

All coordinates in this source package are 0-based, half-open. The diagnostic
FASTQ is an exact-record extraction from the nested original HiFi source,
without sequence or quality edits; its gzip envelope is a new deterministic
archive. It is not a biological replicate.
