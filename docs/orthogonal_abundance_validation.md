# Orthogonal abundance validation gate

Frozen 2026-09-09 after the Ey15-2 and Macadamia donor-matched
assembly-reference-proxy analyses. This is the final method-science gate before
manuscript and release packaging. It is not a new algorithm or comparator.

Completed 2026-09-10 without changing families or thresholds. Ey15
`TXF000002`/`TXF000154` and Macadamia `TXF000496` satisfy the frozen residual-
under-representation rule; Ey15 `TXF001517` and Macadamia
`TXF000563`/`TXF000695` are unresolved. No candidate satisfies a stable cross-k
quantification-bias rule. The compact evidence is
`paper/evidence/orthogonal_abundance_validation_v1`; method development is now
stopped for this release.

## Scientific question

For each already source-eligible repeat family, determine whether the TandemX
HiFi abundance estimate is supported by a sequencing platform that did not
generate the newer assembly. The two competing explanations are:

1. the newer assembly remains under-represented for the family, consistent with
   residual tandem-repeat collapse; or
2. the HiFi read estimator is biased upward.

The newer assembly is an assembly measurement, not absolute copy-number truth.
The primary continuous estimand is therefore the **read--assembly abundance
deficit**, not prediction error against an assembly truth.

## Frozen materials

The validation reuses the complete, unchanged TandemX catalogues and the
5/15/50-kb source-eligibility denominators from:

- `ey15_donor_matched_collapse_v1` (2,133 total; 19 primary eligible);
- `macadamia_jansenii_donor_matched_collapse_v1` (1,227 total; 43 primary
  eligible).

The primary validation denominator is every family eligible at the already
frozen 15-kb newer-assembly threshold, not only the previously positive
families. Family definitions, eligibility, k-mer exclusivity rules, assembly
localization, and the 0.6 binary threshold are not retuned after inspecting
orthogonal results.

Before orthogonal result inspection, the frozen HiFi/newer-assembly ratio marks
three deficit candidates in each species: Ey15 `TXF000002`, `TXF000154` and
`TXF001517`; Macadamia `TXF000496`, `TXF000563` and `TXF000695`. Notably,
`TXF000563` was not an old-to-new assembly-proxy positive; it enters because
the scientific question is the read--newer-assembly deficit rather than whether
an older assembly gained sequence. All other eligible families remain required
negative/context rows.

## Independent data

### Ey15-2

`ERR8666067` is a complete PCR-free Illumina HiSeq 3000 paired-end library from
BioSample `SAMEA13018399`. The source paper states that PCR-free DNA was
extracted independently from the same ground tissue sample used for HMW DNA.
This controls accession and tissue pool while changing extraction, library, and
sequencing platform. It is the primary orthogonal abundance evidence. The ENA
BioSample inventory contains this Illumina run and the already used HiFi run,
but no ONT run. Bionano data are structural context, not a quantitative family
abundance measurement.

### Macadamia jansenii

`SRR11191912` is a complete Illumina NovaSeq 6000 paired-end library and
`SRR11191910` is the paper-declared PromethION library for BioSample
`SAMN14217788`. The earlier paper describes one DNA extraction from young leaf
material of clonally propagated tree accession 1005 used for the sequencing
platforms. The update paper describes the HiFi material as the same sample used
previously, but archival identifiers differ; evidence is therefore
paper-same-sample rather than proven same-extraction across the update. The
Illumina library is not PCR-free and retains possible library/GC bias. ONT is a
second orthogonal platform. The three later ENA PromethION pass-read objects
share experiment `SRX7812077` and must not be added to `SRR11191910`, because
they may be alternative partitions of the same reads.

Exact accessions, file bytes, MD5 values, and source limits are archived in
`paper/evidence/orthogonal_abundance_source_audit_v1`.

## Analysis contract

1. Verify complete downloads against ENA bytes and MD5, retain SHA-256, and
   reconcile FASTQ record/base totals with ENA metadata before analysis.
2. Estimate Illumina abundance with frozen family-exclusive diagnostic k-mers
   at k=21 and k=31. Normalize by the mean depth of deterministic empirical controls that occur
   exactly once in both enrolled old and new assemblies and do not occur in any
   catalogue monomer. Disable the HiFi global quality-survival correction for
   the primary short-read estimate; retain any corrected run only as a labelled
   sensitivity. Report control median/MAD/zero fraction as diagnostics; do not
   silently substitute the median for the public estimator's frozen mean-depth
   rule.
   Controls use seed 6101, a 251-bp per-contig candidate stride, the lowest
   200,000 stable SHA-256 priorities, and the first 20,000 candidates verified
   as exactly single-copy in both assemblies.
3. Independently reproduce diagnostic and control k-mer depths with KMC or a
   separately implemented exact counter. This is an implementation check, not a
   biological comparator.
4. For Macadamia, use ONT as a direction-only secondary validation with a
   frozen competitive mapping/occupancy rule, calibrated on the already used
   HiFi reads before inspecting family outcomes. Do not interpret noisy ONT
   exact-k-mer absence as copy-number evidence.
   Tandemize each of the 43 eligible representatives to 250 kb and map
   competitively with minimap2 2.31-r1302 (`map-hifi` or `map-ont`, `-c --eqx
   --secondary=yes -N 5`). Retain `tp:A:P` alignments with at least 500 aligned
   bases and `nmatch/alignment_block >= 0.75`; union query intervals by family
   and exclude a read if primary intervals are assigned to multiple families.
   Raw occupancy is accepted query bp divided by total library bp times genome
   size. Before ONT inspection, freeze a single mapping-efficiency correction as
   the median raw-HiFi-mapping/frozen-HiFi-k-mer abundance ratio across eligible
   families with positive values; divide raw ONT occupancy by that factor.
5. Report Macadamia abundance under predeclared haploid-genome-size
   sensitivities of 616, 653, 738, and 780 Mb. The 616- and 653-Mb values are
   the source paper's stLFR and Illumina GenomeScope estimates; 738 Mb is the
   newer assembly span; 780 Mb is the earlier cytometric estimate used in the
   original TandemX analysis. This sensitivity is mandatory because 780/653 is
   a 1.194-fold normalization difference.
6. Preserve family-level results and an `unresolved` state whenever platform,
   k, control, or genome-size sensitivities disagree.

## Frozen interpretation rules

Binary classification and continuous magnitude are separate outputs.

- `orthogonal_supports_residual_collapse`: the orthogonal read estimate leaves
  `new_assembly_bp / orthogonal_read_estimated_bp < 0.6`, the HiFi and
  orthogonal directions agree, and the decision is stable across k and the
  predeclared applicable genome-size range. Macadamia requires directionally
  concordant ONT support for the strongest wording.
- `orthogonal_supports_quantification_bias`: the orthogonal estimate is
  consistent with the newer assembly (`ratio >= 0.6`) while the frozen HiFi
  result implies a deficit, with HiFi abundance at least 1.5-fold above the
  orthogonal estimate.
- `unresolved`: neither rule passes or the sensitivity analyses disagree.

For each family report the orthogonal read estimate, newer-assembly localized
bases, their ratio, and
`max(orthogonal_read_estimated_bp - new_assembly_bp, 0)` as estimated
under-representation. Do not score the latter as missing-bp prediction error;
no independent physical copy-number truth is available. Aggregate binary
counts and magnitude summaries must remain separate.

## Stop rule

If the validation supports residual collapse, freeze the scientific method and
proceed to final manuscript, figure, Bioconda, Zenodo, and release work. If it
supports quantification bias, repair only the demonstrated quantification
problem, repeat this gate without changing family eligibility or interpretation
thresholds, and then freeze. No other method development is in scope.
