# Retrospective biological collapse validation

This document defines the source-eligibility gate for testing whether read-based
TandemX estimates predict tandem arrays that are absent or shortened in an older
assembly and present in a newer, more complete assembly. A cultivar name or one
BioProject is insufficient donor identity. The read library, historical assembly
and newer assembly must be tracked separately.

## Evidence tiers

| Tier | Required relationship | Permitted interpretation |
| --- | --- | --- |
| A | Same BioSample or explicit paper statement that the same DNA extraction or individual produced the reads and both assemblies | Donor-matched retrospective truth candidate |
| B | Same study and reported accession/cultivar, but different or absent BioSample identifiers | Material-matched sensitivity analysis; biological differences may mimic collapse |
| C | Same cultivar/strain name across studies or repositories | Exploratory locus nomination only |
| D | Species match only | Ineligible for collapse validation |

No current candidate reaches Tier A. The exact ENA query outputs and their
filtered PacBio genomic-WGS rows are hash-frozen in
`paper/evidence/retrospective_collapse_source_audit`.

## Current candidate audit

| Priority | Material and reads | Historical assembly | New assembly | Current tier and blocker | Local state |
| --- | --- | --- | --- | --- | --- |
| 1 | Rice Nipponbare `SRR25241090`, BioSample `SAMN36368305`, 32,966,159,623 bp | IRGSP-1.0 `GCA_001433935.1`, BioSample `SAMD00000397` | AGIS1.0 `GCA_034140825.1`, BioSample `SAMN36344332` | B/C: same reported cultivar/project context, three different BioSamples; exact donor unresolved | New complete assembly and read/QC data present on T7; historical assembly not yet enrolled |
| 2 | Maize Mo17 seven selected CCS WGS runs from `SRR15447414`--`SRR15447421`, BioSamples `SAMN20604742`--`SAMN20604748`, 151,123,941,618 bp total | Mo17ref_V1, candidate accession `GCA_003185045.1` | T2T Mo17 `GCA_022117705.1`, BioSample `SAMN20854702` | B/C: same inbred/study context but read and assembly BioSamples differ; historical accession must be source-receipt verified before download | New T2T assembly present; only `SRR15447419` is currently complete/QC on T7 |
| 3 | Arabidopsis Col-0R `ERR8666127`, BioSample `SAMEA13018400`, 17,746,722,015 bp | TAIR10.1 `GCF_000001735.4` / GenBank `GCA_000001735.2` | Col-0 assembly `GCA_946499705` | B/C: new reads and assembly are from the same study/accession, but identical plant/DNA extraction is not established; TAIR10 is historical Columbia material | Col-CEN v1.2 reference and Col-0R read/QC data present on T7 |
| 4 | Soybean Williams 82 `SRR23004521`, BioSample `SAMN32622558`, 136,471,999,387 bp | `Wm82.gnm4.4PTR` | `Wm82.a5`, Figshare DOI `10.6084/m9.figshare.22688281` | B/C: same study/cultivar for reads/new assembly, but repository assembly has no matched BioSample and the historical stock relationship is unresolved | No full read/assembly acquisition yet |

The maize run interval above is shorthand, not a claim that every intervening
accession is selected. The authoritative seven-run set is the exact filtered TSV
in the evidence directory.

## Execution gate

1. Acquire and checksum only the smallest complete candidate pair first, with
   official assembly reports and source receipts.
2. Confirm contig naming, total non-N sequence, organellar inclusion, ploidy and
   whether each assembly represents haploid, primary or alternate sequence.
3. Run TandemX on reads without looking at old-to-new array differences; freeze
   family catalogue, diagnostic k-mers, copy estimates and confidence labels.
4. Localize the frozen catalogue independently in both assemblies and compute
   read/assembly representation ratios, missing-array calls, boundary changes and
   control-family false calls.
5. Audit candidate loci with raw read alignments and assembly-to-assembly
   alignment. Separate true sequence gain, gap replacement, haplotype/stock
   difference and mapping ambiguity.
6. Report Tier B/C results as retrospective concordance or locus nomination.
   Reserve biological collapse sensitivity/specificity for Tier A truth or
   independent orthogonal validation.

The first bounded pilot should use rice because the complete new assembly and
read data are already available locally and the genome is much smaller than
maize or soybean. It can test workflow feasibility, but its different BioSamples
prevent a donor-matched accuracy claim.
