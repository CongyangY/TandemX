# YSD56 tandem-repeat pilot: retained results and evidence boundary — 2026-09-10

## Scope

This records the completed bounded YSD56 pilot. It is not a formal abundance
run, a catalogue-completeness result, or evidence that any family is novel.
All execution products are retained under
`/Volumes/T7/Codex/TandemX/results/legume_ysd56_novel_tr_pilot_v1_20260910/`.

The input was the seed-6101 nested HiFi sample `sample_003.fastq.gz` from
*Glycine soja* YSD56: 65,106 reads, 1,099,009,791 bases and nominal 1.0897x
total-base coverage. Its source snapshot binds the donor-matched assembly
`GCA_040083835.1` / `ASM4008383v1`, BioSample `SAMN40909152`, and the input
hashes. The small nested sample is not an independent biological replicate.

## Completed discovery and validation

Discovery used cascade discovery, sequence clustering at 0.95, related-family
audit, 30--1,000 bp periods, minimum span 100 bp, minimum support five, Rust
k-mer backend and four threads. It reported 77,559 candidate reads and 1,474
operational families. `tandemx validate` then validated seven discovery output
files containing 241,099 records.

The retained profiler reports 435.15 s discovery wall time and 1,963.88 MiB
peak aggregate live-process-tree RSS; validation took 1.82 s. These are host-
and input-specific execution observations. The catalogue contains operational
representatives, related-family warnings and candidate hierarchy edges; it does
not establish biological family identity or a validated HOR.

## Whole-assembly localization

The first localization attempt supplied the NCBI `.fna.gz` filename and stopped
at input validation because that suffix was not then accepted. It produced no
localization result. After the input-extension repair, the identical recorded
assembly and catalogue completed successfully: `tandemx locate` wrote 18,131
arrays and 520 compatibility-comparison rows, and `tandemx validate` validated
the three localization outputs containing 119,330 records. The successful
profile reports 606.65 s wall time and 141.97 MiB peak aggregate
live-process-tree RSS.

This is a descriptive full-assembly localization of all 1,474 operational
families. It is not read--assembly copy-number comparison: no real-data
quantification was run, and the localization compatibility file records missing
read estimates as zero with an explicit `missing_read_estimate` warning. Array
absence under this exact-k-mer profile, or an `assembly_only` row, is not proof
of biological absence.

## Known-repeat annotation and bounded exclusion

An initial `tandemx annotate-repeats` run against the recorded soybean GenBank
clone library was stopped before completion because its exhaustive Python local
alignment did not finish promptly. It is retained as
`run/annotate_known_genbank_attempt1_naive_stopped/`; it has no annotation table
and must not be used as negative evidence.

The replacement bounded exclusion used a fixed 13-record accessioned clone
library and edlib HW alignment of the complete shorter sequence against a
doubled longer sequence in both orientations. It retained its input hashes,
script hash, thresholds (strong >=0.90, possible >=0.80 and minimum aligned
length 50 bp) and pairwise output. Across the 1,474 families it reported:

| State | Families | Meaning |
| --- | ---: | --- |
| `strong_known_clone_match` | 53 | Exclude from a follow-up novelty screen against this limited library. |
| `possible_known_family_match` | 65 | Retain as related/ambiguous, rather than calling novel. |
| `weak_known_similarity` | 9 | Weak similarity only. |
| `no_match_in_limited_library` | 962 | No qualifying match in these 13 records; not proof of novelty. |
| `insufficient_length_for_exclusion` | 385 | The family was too short for this exclusion decision. |

The library has known coverage gaps, including assembly-derived labels without
independently accessioned monomers. Thus neither a no-match nor an
insufficient-length state can support a novelty claim.

## Evidence-integrated post hoc candidate triage

A machine-readable post hoc rule was locked before its final successful output
in `benchmarks/configs/legume_ysd56_posthoc_candidate_evidence_v3.json`. It
integrates family confidence and complexity, hierarchy/redundancy relations,
the 13-record soybean clone exclusion, read support span, whole-assembly array
span, both centromere definitions from the source paper, distance from
chromosome ends and representative gene annotation. The rule requires high
confidence, `low_complexity_flag=false`, no emitted hierarchy or redundancy
edge, `no_match_in_limited_library`, at least 50 kb read support span, a longest
assembly locus spanning at least ten monomers after merging intervals separated
by at most 100 bp, no overlap with either published centromere definition, at
least 100 kb from either chromosome end and no gene overlap.

Two technical attempts were retained rather than hidden. Version 1 terminated
before output because the analysis reader did not parse TandemX's structured
`family_id=...` FASTA headers. Version 2 computed rows but terminated during
output because three internal gate names did not match the declared TSV schema;
its partial directory is retained as
`run/candidate_evidence_v2_attempt2_output_schema_failure/`. Regression tests
were added before the successful v3 execution.

The successful v3 table contains all 1,474 families and exactly two passing
rows:

| Family | Monomer | Pilot read evidence | Assembly locus (BED) | Located repeat support | Context | Initial state |
| --- | ---: | --- | --- | ---: | --- | --- |
| `TXF000708` | 785 bp | 8 unique reads; 94,397 bp span; mean identity 0.9933 | chr13 / `CP154579.1:13,966,881-14,090,910` | 124,029 bp | Outside both centromere definitions; 13.97 Mb from the nearest end; zero gene overlap; flanked by `SoyYSD56_chr13G04070` (906 bp) and `SoyYSD56_chr13G04080` (947 bp) | preliminary follow-up candidate |
| `TXF000367` | 335 bp | 16 reads; 217,921 bp span; mean identity 0.9945 | chr19 / `CP154585.1:13,468,027-13,538,217` | 70,132 bp in two exact-k-mer intervals separated by 58 bp | Outside both centromere definitions; 13.47 Mb from the nearest end; zero gene overlap | preliminary follow-up candidate before rDNA exclusion |

The full evidence table, shortlist and receipt are in
`run/candidate_evidence_v3/`. The merged chr19 locus spans 70,190 bp including
the 58-bp internal gap; the repeat-supported intervals themselves total 70,132
bp. These coordinates are exact-k-mer localization outputs and do not establish
complete biological array boundaries.

## rDNA and organelle exclusion changes the shortlist

A fixed 117-record TideCluster rDNA library was applied locally. The first run
terminated before output because this public FASTA contains IUPAC `W/Y` bases,
which the strict production reader rejects. The corrected benchmark-only run
explicitly masked 11 non-ACGT library bases as `N` and recorded this operation.
It found that `TXF000367` contains a reverse-orientation, zero-edit 119/119-bp
match to a *Lupinus luteus* 5S rDNA sequence. `TXF000367` is therefore excluded
from the previously-unreported queue and assigned `rDNA_related_known_repeat`.

`TXF000708` had no qualifying rDNA match; its best rDNA identity was 0.559006.
It also had no qualifying match to the public *Glycine soja* chloroplast
`NC_022868.1` or mitochondrial `NC_039768.1` references; the best organelle
identity was 0.542675 to the chloroplast record. A separate fixed comparison
against the five published soybean retrotransposon-associated minisatellites
`Gm_ms_a`--`Gm_ms_e` also found no qualifying match; the highest identity was
0.653846 to the 26-bp `Gm_ms_a` motif. These are bounded exclusions against
fixed libraries, not proof of non-rDNA, non-organelle or non-TE origin. After
the rDNA refinement, the retained FASTA and evidence row are in
`run/candidate_refined_after_rdna_v1/`, and only `TXF000708` remains.

The historical SoyTEdb bulk FASTA endpoint failed with a server-side PHP
memory-exhaustion error, and its small public reproducibility repository does
not contain the generated bulk library. A separately frozen local screen used
the Dfam 4.0 API collection of 1,697 raw or curated *Glycine max* consensus
families (3,540,456 bp). `TXF000708` again had no qualifying match: its best
pair was a 72-bp Dfam record at 0.694444 glocal edit identity, below the locked
0.80 possible-match threshold. This adds a species-specific TE exclusion layer
but does not replace the unavailable SoyTEdb library or establish novelty.

An additional source audit found seven directly published soybean or
wild-soybean tandem-repeat sequences that were absent from the original
13-clone library: `Z26334.1`, `AF297983.1`--`AF297985.1`, and published
`SBRS1`--`SBRS3` sequences. A separately locked local screen retained
`TXF000708` as `no_match_in_limited_library`; its best comparison was the
92-bp `SBRS1` sequence at 0.586957 glocal edit identity, below the 0.80
possible-match threshold. This closes those seven accession/source omissions,
but `CentGm-2`, `CentGm273` and `CentGm444` remain unresolved because no
independently located public monomer record was obtained.

The YSD56 source paper searched assembly tandem repeats with TRF periods from
30 to 500 bp and highlighted trf91, trf92, trf182, trf183, trf184, trf273 and
trf276. A 785-bp operational monomer was therefore outside that declared
period range. This explains why the present candidate is not a direct recovery
of one of the seven reported length labels. The official 1,172,429-byte
supplementary DOCX (SHA-256
`37926502c611eab322b91ffb2b484bb4ccee1d632ed4596b73e5f2b887d3333d`)
also names only those seven families and contains no `785` text occurrence.
This supports `not_reported_in_the_YSD56_article_or_supplement`, but does not
show that the sequence or array was absent from all previous literature or
databases.

## Nested 10.9570x depth recurrence

A second configuration was committed before inspecting output and applied to
the 10.9570x nested sample from the same HiFi run: 654,450 reads and
11,050,418,089 bases. Discovery completed with 778,295 candidate intervals and
36,847 operational families, and validation accepted seven output files with
2,694,384 records. The profiler recorded 5,335.51 s wall time, 13,248.39 s
user CPU, 1,517.28 s system CPU, 5,509.84 MiB peak aggregate process-tree RSS
and 460,330,868 peak scratch bytes. These measurements expose a real high-depth
catalogue and clustering cost; the 36,847 operational representatives must not
be interpreted as 36,847 biological repeat families.

The fixed 1.0897x candidates were then compared with all 36,847 deeper
catalogue representatives under the predeclared full-shorter circular glocal
rules. Both had exact direct recurrence. `TXF000708` matched deeper family
`TXF000598` at 785/785 bp, zero edits and identity 1.000000. The deeper family
has 97 supporting reads, 1,397,702 bp supporting span, mean identity 0.9905,
high confidence and one emitted possible higher-order/partial relationship.
The 335-bp rDNA-related control likewise recurred exactly as `TXF000248`.
Because the 10.9570x reads contain the 1.0897x sample, this is same-run depth
recurrence rather than independent biological or platform replication.

## Independent detector confirmation of the assembly array

The locked assembly interval `CP154579.1:13,966,881-14,090,910` (zero-based,
half-open) was extracted from the hash-bound YSD56 assembly. Its length is
124,029 bp. TRF 4.10.0-rc.2 reported 461 records whose merged coverage spans
124,029 bp and whose dominant period is exactly 785 bp. TideHunter 1.5.5
reported two records whose merged coverage also spans 124,029 bp and whose
dominant period is exactly 785 bp. Both tools therefore passed the predeclared
minimum 0.90 locus-coverage and period-concordance rules and returned
`confirmed_tandem_structure`.

This confirms that the assembly locus is a continuous 785-bp tandem array; it
does not establish an HOR, centromeric function, physical array size, donor-
independent support or sequence novelty. Together with exact deeper-catalogue
recurrence, the retained status is
`stable_previously_unreported_in_source_article_candidate`, not `novel_TR`.

## What remains unestablished

The pilot does not provide any of the following:

1. abundance estimates from the complete YSD56 HiFi collection or a justified
   real-data normalization path;
2. independent ONT or Illumina support, an independent specimen, or a
   catalogue-completeness estimate; the 10.9570x recurrence is nested within
   the same HiFi run;
3. exclusion against a comprehensive soybean/legume TE/repeat resource, a
   full nucleotide-database search, or complete repeat annotation of the array
   context; the current soybean-clone, seven-sequence supplemental soybean,
   rDNA, organelle, five-minisatellite and Dfam *G. max* consensus libraries
   are bounded screens only, `CentGm-2`, `CentGm273` and `CentGm444` remain
   sequence-unresolved, and the SoyTEdb bulk FASTA remains unavailable;
4. evidence that the 18,131 localized arrays are complete, or that unlocalized
   families are absent from the donor genome; or
5. a validated HOR, physical copy number, assembly deficit, or biological
   novelty result.

The next evidence gate is a predeclared full-data contract: freeze the
complete-read input hashes, real-data normalization rationale, family selection
rule, library versions and search thresholds, independent-evidence plan,
resource stops and a new output directory before running quantification,
formal localization or comparison. It must preserve all failure and empty
states and retain assembly representation as a proxy rather than physical truth.
