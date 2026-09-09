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

Ey15-2 now reaches Tier A because the primary paper explicitly compares CLR and
HiFi assemblies of the same sample. This corrects the earlier source audit,
which enrolled the Ey15-2 read metadata but considered the Col-0 historical
comparison instead. The exact ENA query outputs and their filtered PacBio
genomic-WGS rows remain hash-frozen in
`paper/evidence/retrospective_collapse_source_audit`. The analysis is
preregistered in `benchmarks/configs/ey15_donor_matched_collapse_v1.json` and
documented in `docs/ey15_donor_matched_validation.md`.

Macadamia jansenii provides a second species-level Tier A reference proxy. The
HiFi update paper explicitly describes its material as the same sample used for
the earlier CLR comparison. The archival records nevertheless use different
BioSamples, aliases and collection dates, so the allowed claim is paper-level
same-sample matching; the same DNA extraction is not established. The newer
assembly also shares the HiFi evidence used by TandemX and is not absolute copy
truth.

## Current candidate audit

| Priority | Material and reads | Historical assembly | New assembly | Current tier and blocker | Local state |
| --- | --- | --- | --- | --- | --- |
| 1 | Arabidopsis Ey15-2 `ERR8666125`, BioSample `SAMEA13018399`, 18,636,790,429 bp | Bionano-scaffolded `9994.CLR_Canu` | Bionano-scaffolded `9994.HiFi_Hifiasm`; final HiFi-Hifiasm+CLR-Canu hybrid as sensitivity | A: the paper explicitly compares CLR and HiFi assemblies of the same Ey15-2 sample; newer assembly shares HiFi evidence and is a high-quality reference proxy, not absolute independent truth | Complete frozen analysis and compact evidence archive |
| 2 | Macadamia jansenii `SRR13557763` and `SRR13557762`, BioSamples `SAMN17524927/8`, 22,546,488,654 bp | 84x CLR Falcon-Unzip primary contigs after Purge Haplotigs, tree accession 1005 and BioSample `SAMN14217788` | two-cell HiFi IPA primary assembly | A reference proxy: update paper explicitly states same sample; archival identifiers, aliases and collection dates differ, so same extraction is not claimed | Publisher/ENA checks, complete frozen analysis, both depth normalizations and independent verifiers passed; compact evidence archive complete |
| 3 | Rice Nipponbare `SRR25241090`, BioSample `SAMN36368305`, 32,966,159,623 bp | IRGSP-1.0 `GCA_001433935.1`, BioSample `SAMD00000397` | AGIS1.0 `GCA_034140825.1`, BioSample `SAMN36344332` | B/C: same reported cultivar/project context, three different BioSamples; exact donor unresolved | New complete assembly and read/QC data present on T7; historical assembly not yet enrolled |
| 4 | Maize Mo17 seven selected CCS WGS runs from `SRR15447414`--`SRR15447421`, BioSamples `SAMN20604742`--`SAMN20604748`, 151,123,941,618 bp total | Mo17ref_V1, candidate accession `GCA_003185045.1` | T2T Mo17 `GCA_022117705.1`, BioSample `SAMN20854702` | B/C: same inbred/study context but read and assembly BioSamples differ; historical accession must be source-receipt verified before download | New T2T assembly present; only `SRR15447419` is currently complete/QC on T7 |
| 5 | Arabidopsis Col-0R `ERR8666127`, BioSample `SAMEA13018400`, 17,746,722,015 bp | TAIR10.1 `GCF_000001735.4` / GenBank `GCA_000001735.2` | Col-0 assembly `GCA_946499705` | B/C: new reads and assembly are from the same study/accession, but identical plant/DNA extraction is not established; TAIR10 is historical Columbia material | Col-CEN v1.2 reference and Col-0R read/QC data present on T7 |
| 6 | Soybean Williams 82 `SRR23004521`, BioSample `SAMN32622558`, 136,471,999,387 bp | `Wm82.gnm4.4PTR` | `Wm82.a5`, Figshare DOI `10.6084/m9.figshare.22688281` | B/C: same study/cultivar for reads/new assembly, but repository assembly has no matched BioSample and the historical stock relationship is unresolved | No full read/assembly acquisition yet |

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

The first bounded validation used Ey15-2, and Macadamia now supplies a second
species-level paper-same-sample reference proxy. Both newer assemblies share
HiFi evidence with the estimator, and neither is absolute repeat-copy truth.
Rice remains a useful workflow sensitivity, but its different BioSamples
prevent a donor-matched accuracy claim.
