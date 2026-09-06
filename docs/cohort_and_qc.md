# Real and simulated cohort, scale and QC programme

Updated 2026-09-06 following the user's explicit requirement for substantial
real data, simulations, taxonomic breadth and literature-informed QC. This is
the prospective design and acquisition record, **not completed validation**.

Full-file QC passed for ten included libraries across eight species: Mo17,
Col-0N, Col-0R, Ey15-2R, Nipponbare, Morex, Lo7, Chinese Spring, Victoria and
YSD56. Together they contain 14,937,608 reads and 262,731,255,175 bp. These are
file-QC records, not ten independent specimens or eight completed biological
accuracy validations. Seed6101 whole-library sampling is complete for all ten
libraries; the largest nested samples range from 1.109 Gb for Victoria to
42.360 Gb for Lo7. FASTQ and ID products remain on T7,
while compact plans, receipts and distributions are archived with paper evidence.
Mo17 comparison includes 1.129 Gb, with all three methods successful; its seven
TandemX products are identical before/after dictionary-index optimization. All
real-call metrics remain descriptive. The original-IPK MorexV3 pseudomolecule
FASTA now matches its published SHA-256 and passed streaming FASTA QC; raw and
reference donor identity and satellite-array truth remain unresolved.
YSD56 wild-soybean HiFi passed exact download, full-file QC and deterministic
whole-library nested sampling; no accuracy inference has been made.
The prospective enrollment table below retains metadata-level source details;
use the latest receipts and `current_status.md` for completed versus running work.
Published Mo17 region evidence is documented in
[published_mo17_regions.md](published_mo17_regions.md).

The [PanOat primary-source audit](../paper/evidence/PanOat_source_audit/README.md)
now verifies Victoria's HiFi method, exact raw sample and assembly project from
the original supplement. Its Iso-Seq yields are explicitly excluded from genomic
coverage. Public reviewer comments motivate clear sample units, alternative
quantification under homoeolog ambiguity, model diagnostics and biological cases
that directly test the resource. They do not impose a universal sample-size rule.

## Scope and independent units

The development target is **8–10 plant species and at least 20 independent
materials**, covering small and very large genomes, monocots/eudicots, and
diploid/polyploid contexts. Prioritize multiple materials in Arabidopsis, rice,
maize, barley, rye, wheat, oat and soybean. These are project design targets,
not asserted journal rules. Final enrollment and precision must be justified
from between-material/species variability before freezing the publication test.
At least two independently sourced materials per core species should be sought;
key large-genome species require more than one reference accession. Missing
source data must appear as missing coverage of the design, not a completed gate.

A BioSample ID, sequencing run or SMRT cell is **not automatically** a biological
replicate. Resolve cultivar/accession, donor plant or pool, tissue, extraction,
library and run relationships from the original paper/supplement and XML.
The wheat project has many run-labelled BioSamples; the Morex project has
multiple BioSamples for one cultivar. Never inflate specimen or species counts.
Different accessions in one species test within-species generalization; batches
from the same DNA test technical reproducibility. Pooled plants remain pooled.

## Primary-source lessons and corresponding actions

| Primary study | Relevant design/QC lesson | TandemX action |
| --- | --- | --- |
| [SRF, Genome Research](https://genome.cshlp.org/content/early/2023/11/02/gr.278005.123.full.pdf) | Four Arabidopsis HiFi datasets at about 40×, broader species experiments, task-specific short-read TAREAN comparison; k=101 rescues maize CentC missed at k=151 | Match inputs, full-workflow cost and coverage-dependent count cutoffs; retain k sensitivity, donor/pool distinctions and out-of-range HORs |
| [Mo17 T2T, Nature Genetics](https://www.nature.com/articles/s41588-023-01419-6) | Independent technologies, k-mer QV/completeness, local coverage and orthogonal rDNA/FISH evidence complement a complete assembly | Use same-material reads/assembly; triangulate arrays with ONT/Illumina and published assays where available; separate evidence types |
| [Merqury, Genome Biology](https://doi.org/10.1186/s13059-020-02134-9) | Read k-mer spectra support assembly quality/completeness and copy representation checks | Include spectra-cn, reliable-k-mer completeness and QV as background QC; these do not establish exact satellite copy truth |
| [HiTE, Nature Communications](https://www.nature.com/articles/s41467-024-49912-8) | Nine species, curated references, several evaluation definitions and different assembly versions | Borrow evaluation discipline, not TE accuracy endpoints; report full/fragmented/duplicated sequence recovery and assembly-version sensitivity |
| [Barley pangenome, Nature](https://www.nature.com/articles/s41586-024-08187-1) and [PanOat, Nature](https://www.nature.com/articles/s41586-025-09676-7) | Multiple materials and heterogeneous assembly histories; PanOat discloses a chromosome-orientation issue | Resolve exact version/material and known corrections; don't treat a pangenome as interchangeable truth genomes |

BUSCO gene completeness alone does not validate satellite arrays. Likewise,
read mapping rate alone can hide multimap ambiguity, collapsed copies, GC/GA
dropout and organellar contamination. The [Genome Research sequencing-bias
study](https://genome.cshlp.org/highwire_display/entity_view/node/1046439/full)
motivates explicit read initiation, length and quality diagnostics near repeats;
its Drosophila findings are not assumed to hold identically in plants.

[MorexV3 gap analysis](https://onlinelibrary.wiley.com/doi/10.1111/pbi.13816)
provides an especially relevant orthogonal case: HiFi/ONT/Illumina, optical maps,
CENH3 and FISH were compared with pseudomolecules. Its methods resolve the raw
projects as PRJEB40587/40588/31444 and reference PRJEB40589. Repeat quantities
are method-dependent estimates, not exact truth labels. Preserve cpDNA correction,
genome-size uncertainty, chrUn versus pseudomolecule content and differences
between technical read batches. The paper uses HvT01 X16095.1:1–118; its 120-bp
5S gene marker is not a complete intergenic-spacer-containing 5S repeat unit.
Known-motif BLAST+union and TRF comparisons complement de novo tools here.

## Verified source entry points and current state

| Species/material | Verified project/run or source | Current boundary |
| --- | --- | --- |
| Arabidopsis Col-0N / Col-CEN | PRJEB46164, ERR6210723; [official Col-CEN v1.2](https://github.com/schatzlab/Col-CEN/tree/main/v1.2) | Full file QC passed: 933,904 reads / 14,646,601,458 bp; whole-library nested samples completed. Pooled plants, not individual replicates. Reference and known issues verified; nuclear/organelle denominators separate |
| Arabidopsis Col-0R / Ey15-2R | ERR8666127 / ERR8666125 from SRF Table 3 | Complete-file QC passed: 17.747/18.637 Gb. Seed6101 nested sampling completed through 10.644/10.628 Gb. Col-0R is reported as one plant; Ey15-2R is pooled. These are two libraries/materials, not species or interchangeable replicates |
| Rice Nipponbare | PRJNA953663, SRR25241090; [reference audit](nipponbare_reference.md) | Complete 32.966-Gb FASTQ QC and GCA_034140825.1 reference QC passed. Raw/assembly sample dates differ; exact donor match unresolved. chr9 rDNA model sequence excluded from exact copy truth |
| Maize Mo17 | PRJNA751841, seven PacBio WGS runs | First full CCS batch SRR15447419 passed QC (407,670 reads / 5.625 Gb); complete reference and whole-library samples verified. A single batch is not the full published data set |
| Maize B73 | SRR11606869 from SRF | 48.075 Gb bases, 41.984 GB compressed in metadata; matching assembly/material review pending |
| Barley Morex and pangenome | PRJEB40587, PRJEB57567, PRJEB58554; MorexV3 DOI 10.5447/ipk/2021/3 | ERR4659246 complete 19.525-Gb FASTQ QC and nested sampling through 10.745 Gb passed. The original MorexV3 FASTA matches the published SHA-256; eight records contain 4,225,605,719 bases. Same-study context is not identical-donor or exact satellite-copy truth. Pangenome material enrollment remains separate |
| Rye Lo7 | PRJEB91463; [source paper](https://www.nature.com/articles/s41467-026-76753-4) | One of three HiFi Revio technical runs, ERR15194059, passed complete 77.092-Gb FASTQ QC; seed6101 nested sampling completed through 42.360 Gb. Single diploid Lo7 plant is explicit; version-matched reference still needed |
| Bread wheat Chinese Spring / CS-IAAS study | PRJNA1062539; [source paper](https://www.nature.com/articles/s41588-025-02137-x) | SRR28200549 passed complete QC at 1,500,000 reads/24.954 Gb; seed6101 nested sampling completed through 11.215 Gb. The study's 80 PacBio WGS records are technical batches. Exact donor/assembly matching remains unresolved |
| Oat PanOat | PRJEB56828; [source paper](https://www.nature.com/articles/s41586-025-09676-7) | Victoria ERR10422581 passed complete QC at 398,850 reads/7.346 Gb; seed6101 sampling completed through 1.109 Gb. Raw SAMEA111508775 and assembly GCA_947311595 are source-supported. ERR10422482/CN25955 has A. fatua in ENA versus A. occidentalis in paper figure labels; retain the discrepancy and do not substitute it for cultivated oat |
| Soybean ZH13, Jack, Wm82 | [Genome Biology source](https://doi.org/10.1186/s13059-025-03924-9); ZH13 PRJCA015269 | Data source verified in paper; files, material/version matching and accessibility still pending. New relevant comparator: unitFinder |
| Wild soybean YSD56 | PRJNA1095640, SRR28726931; [source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12350797/); assembly GCA_040083835.1 | Exact ENA size/MD5/SHA-256 and full gzip/FASTQ QC passed: 2,617,227 reads/44.193 Gb, no duplicate IDs. Seed6101 nested sampling completed at 11.670 Mb/0.0116×, 110.438 Mb/0.1095×, 1.099 Gb/1.0897× and 11.050 Gb/10.9570× using the exact 1,008,523,555-bp assembly denominator. Cultivated and wild soybean remain distinct species/material strata; file QC and nominal total-base coverage are not donor, nuclear-depth or accuracy validation |

GB means 10^9 file bytes; Gb means 10^9 sequence bases. These totals are ENA
metadata observations, not downloaded or QC-passed amounts. Metadata snapshots
are at the T7 data root `data/manifests/cohort_screen_20260906`.

Selected run/sample/experiment XML is archived with the cohort evidence. It
confirms Mo17 cultivar and CCS-labelled experiment, Nipponbare/AGIS-1.0, and a
single diploid Lo7 plant with HiFi library metadata. It does not turn technical
batch identifiers into independent biological specimens. Col-CEN v1.2 nuclear
Chr1-5 total 131,559,676 bp; ChrM/ChrC add 521,402 bp. Both partitions and the
published issues file must be retained during downstream QC.

## Data quantity and scale experiments

1. Full-file validation precedes analysis sampling: size, MD5, SHA-256, gzip
   trailer, record/sequence-quality agreement and observed read/base totals.
   Never substitute an archive-order prefix for a random whole-library sample.
2. Use seeded hash sampling across each complete included library, independent
   of repeat status. Keep exact read IDs/hashes so all tools see identical inputs.
   Report the included batches and unavailable batches. Whole-study uniform
   sampling cannot be claimed after choosing only some sequencing batches.
3. Coverage ladder: nominal 0.1/0.25/0.5/1/2/5/10/20/40× where available, plus
   the full included dataset. Record achieved bases and nuclear genome-size
   denominator; organellar bases do not establish nuclear coverage. Don't create
   higher coverage by duplicating reads. Three sampling seeds test sampling
   stability, not additional specimens.
4. Separate equal-base throughput from equal-coverage detection/quantification.
   Measure 10 Mb, 100 Mb, 1 Gb, 10 Gb and larger inputs, including at least a
   completed large-genome run in each claimed genome-size domain. Whole-genome
   assembly scans must include the full chromosome files, not only easy loci.
5. Record wall/CPU, parent-and-child aggregate peak memory, scratch **peak**,
   threads, I/O, stage costs, failures and resume behavior. Single-thread and
   2/4/8-thread tests use available hardware; unsupported sizes/timeouts are
   reported. Download/QC/compilation must not overlap final timing runs.

## QC and truth hierarchy

File QC implementation: `fetch_ena_complete.py` resumes only validated byte
ranges, enforces the expected size budget and verifies full MD5 before renaming.
`qc_complete_fastq.py` parses every read, checks gzip completion, verifies exact
ID uniqueness with disk-backed SQLite and emits length/N50, GC/N and quality
distributions. Reported Phred scores are not measured biological accuracy.
Archive-generated unique IDs cannot detect duplicate molecules after renaming;
original ZMW IDs or alignment/sequence-based duplicate evidence need a separate
check. No automatic filtering or quality-score correction is applied here.

Biological QC remains required: material/assembly identity, chromosome/organellar
partition, contaminant evidence, empirical depth and GC/GA biases, read length
and quality by repeat context, assembly gaps/known issues, k-mer spectra and
orthogonal read support. Report nuclear and organellar denominators separately.
Alignment mapping quality alone must not discard all genuine repeat evidence.

Truth tiers: (1) planted genome and observed read-coordinate truth; (2) real
background with independently planted/removed arrays; (3) same-material assembly
plus orthogonal long reads or published assays; (4) curated known-repeat recall;
(5) descriptive de novo candidates. Assembly-derived labels are not independent
of reads used to assemble it. No complete real-family recall without an explicit
curated denominator; no FISH success rate without actual independent assays.
Blind/manual candidate review must use a declared sampling frame and criteria.

## Simulation expansion and inference

The 100-read challenge and 199.1-kb abundance genomes remain diagnostic pilots.
Expand to independent monomer families, factorial length/copy counts, diverse
GC/low complexity, substitutions/indels/homopolymer errors, biological unit
divergence, interruptions/TE insertions, related families/HORs, multiple arrays,
mixtures/ploidy, contaminants and empirical coverage/length/quality distributions.
Add megabase arrays and 10 Mb–1 Gb backgrounds before large-scale claims;
use streaming generation for these larger genomes. Fit empirical simulation
settings on development materials only. Reserve families and entire materials/
species for final tests; AI evaluation additionally needs species-held-out
calibration and ablation. Current confidence labels are not learned/calibrated.

Compute uncertainty over independent specimens/genomes and hierarchical species
groups, not thousands of correlated read calls. Freeze primary endpoints,
superiority/noninferiority margins, tuning budget and multiplicity handling
after pilot variance analysis and before held-out scoring. Report every dataset
and trade-off. Enrollment breadth and millions of reads cannot rescue a biased
truth definition or establish universal superiority.
