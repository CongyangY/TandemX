# Ey15-2 donor-matched collapse validation v1

## Why this experiment is now first priority

Rabanal et al. explicitly compare CLR and HiFi assemblies of the same Ey15-2
sample (1001 Genomes accession 9994; stock CS76399). The methods describe the
Ey15-2 HMW-DNA preparation used for the CLR and HiFi libraries, and ENA assigns
the archived HiFi run `ERR8666125` to BioSample `SAMEA13018399`. This satisfies
the repository's Tier A source rule for a donor-matched retrospective truth
candidate. It corrects the earlier audit, which had considered Col-0 but had not
enrolled the same-sample Ey15-2 CLR/HiFi comparison.

The source-data bundle is Zenodo record `7326462`. The primary contrast will use
the Bionano-scaffolded `9994.CLR_Canu` assembly as the historical reference and
the Bionano-scaffolded `9994.HiFi_Hifiasm` assembly as the newer high-quality
reference proxy. Holding the scaffolding evidence type constant makes this a
cleaner primary comparison than using the authors' final HiFi-Hifiasm plus
CLR-Canu hybrid. The final hybrid is retained as a predeclared sensitivity
assembly.

This is not absolute biological copy truth. The same HiFi evidence used by
TandemX also contributed to the newer assembly, so the two measurements are not
fully independent. Bionano maps, PCR-free Illumina data, author repeat
annotations and assembly-to-assembly alignment are orthogonal audit layers, not
permission to tune the primary analysis.

Primary sources:

- Rabanal et al., *Nucleic Acids Research* 2022,
  <https://doi.org/10.1093/nar/gkac1115>.
- Official supporting-data record, <https://doi.org/10.5281/zenodo.7326462>.

## Frozen inputs and analysis

The machine-readable preregistration is
`benchmarks/configs/ey15_donor_matched_collapse_v1.json`. It was written before
inspection of TandemX old/new localization differences.

1. Discover the catalogue from the already frozen seed6101 6% whole-file hash
   sample: 50,267 reads, 1,119,053,887 bp, SHA-256
   `1f05335279f80e604d0f451a72741d45fb53e31d35e139023e0afbe528b26baf`.
   Use the validated cascade detector, sequence clustering at 0.95 identity and
   the predeclared 30--1,000-bp period range.
2. Quantify the frozen catalogue against the complete 837,586-read,
   18,636,790,429-bp HiFi FASTQ. The primary depth is total read bases divided
   by the independently published 143.12-Mb PCR-free Illumina k-mer genome-size
   estimate. The paper's approximate 107x q20 depth is a sensitivity analysis,
   not a value selected from assembly agreement.
3. Localize the unchanged catalogue independently in the old, new and final
   sensitivity assemblies with the frozen k=21 IID-base model, minimum identity
   0.9.
4. Define the historical binary assembly-transition label as old/new localized
   bp below 0.6 and the corresponding read label as old/read-estimated bp below
   0.6. A family is source-eligible
   for the primary denominator when the newer assembly contains at least 15 kb
   of localized sequence. The 5-, 15- and 50-kb denominator sensitivities are
   fixed in advance.
5. Report every family. Rows below the reference-size gate are
   `not_source_eligible`; a source-eligible family without a positive read
   estimate is a `technical_failure`. Neither category may be silently dropped.
6. Report TP/FN/FP/TN, Wilson intervals, balanced accuracy, MCC, and the
   relationship between HiFi-versus-old deficit and new-versus-old assembly
   gain. This is assembly-proxy concordance, not missing-bp prediction accuracy.
   New/read agreement is reported but is not an eligibility rule, because using
   it to define truth would condition the denominator on the predictor.

## Independent audit and stopping rules

Primary scoring is frozen before author annotations or assembly differences are
used to interpret individual families. Afterwards, author repeat classes,
published centromere/5S totals, PCR-free Illumina evidence, Bionano support and
old-to-new alignments may explain concordant and discordant loci. These audits
cannot change v1 thresholds, discovery parameters or the family catalogue.

A zero-family result, localization failure, reference discordance or resource
failure remains a formal result. No parameter rescue is allowed inside v1. Any
method change requires a separately named development experiment followed by a
new untouched validation input.
