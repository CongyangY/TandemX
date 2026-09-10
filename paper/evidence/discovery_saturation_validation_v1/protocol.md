# TandemX discovery saturation validation protocol

## Objective

This validation asks how repeat-candidate and operational-family discovery
changes with nominal read depth. It is a coverage-saturation experiment for the
unchanged production discovery workflow. It does not tune the detector,
clustering threshold or abundance estimator.

## Frozen design

- Three independently generated 10-Mb haploid source genomes use validation
  seeds 6501, 6502 and 6503.
- Each genome contains the established 55-family factorial design: 54 families
  crossing periods 61/171/421 bp, copy counts 20/80/200, GC 0.3/0.5/0.7 and
  unit divergence 0/0.02, plus one 171-bp, 6,000-copy family.
- Read lengths follow the complete Mo17 HiFi length distribution.
- One fixed HiFi-like IID error model is used: substitution, insertion and
  deletion probability 0.001 each.
- Nominal source coverages are 0.5x, 1x, 2x, 5x, 10x, 20x and 30x.
- Within each source genome, read streams are nested: a lower-depth collection
  is a prefix of the collection at every greater depth. The three genome seeds
  are the independent simulation units; depths within a seed are paired.
- TandemX uses periods 30-1000 bp, minimum repeat span 100 bp, minimum support
  one read, elastic discovery, sequence clustering at identity 0.95, related-
  family audit, the Rust k-mer backend and one thread per run.
- Runs may execute concurrently across seeds because runtime is not a primary
  endpoint. Timings are descriptive and will not be used for performance claims.

The design is stored in
`benchmarks/configs/discovery_saturation_validation_v1.json`. Validation seeds
are distinct from previous development seeds 6301-6303, abundance-validation
seeds 6401-6403 and reserved held-out seeds.

## Pre-specified metrics

For every seed and depth, report:

1. raw candidate count from the discovery summary;
2. unique candidate-monomer count after cyclic and reverse-complement
   canonicalization;
3. operational family count, verified against emitted representatives;
4. planted-family recall overall and within abundance tiers;
5. wall time and direct-child peak RSS as descriptive run diagnostics.

Abundance tiers are defined from planted copy count before results are viewed:

- low: 20 copies (18 families);
- medium: 80 copies (18 families);
- high: at least 200 copies (19 families, including the 6,000-copy family).

Family recovery uses one-to-one maximum matching at cyclic global Levenshtein
similarity at least 0.90, including reverse-complement equivalence. The threshold
matches the established independent family-recovery endpoint and is not the
production clustering threshold.

For every adjacent depth pair within a seed, operational catalogues are matched
with the same one-to-one sequence rule. Let *M* be the number of matched
families, *L* the lower-depth family count and *H* the higher-depth family count.
Report:

- family-level Jaccard = *M* / (*L* + *H* - *M*);
- new families = *H* - *M*;
- lost families = *L* - *M*;
- new-family proportion = (*H* - *M*) / *H*;
- new families per added 1x = (*H* - *M*) / depth increment;
- newly recovered planted families and their number per added 1x.

## Operational saturation rule

An adjacent transition passes when its new-family proportion is below 0.05 and
family-level Jaccard is above 0.95. A depth is declared saturated only when the
transition ending at that depth and the immediately preceding transition both
pass in all three validation seeds. The earliest depth meeting this rule is the
primary saturation depth. Median and range curves are descriptive and cannot
override the all-seed rule.

A finite, reproducible primary saturation depth is eligible for the main text.
Failure to meet the rule, unstable high/medium-family recall, or a result driven
by only one seed remains a Supplementary result. All runs and metrics are
retained in either case.
