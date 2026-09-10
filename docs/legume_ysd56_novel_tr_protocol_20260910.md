# YSD56 legume novel-TR pilot protocol

## Purpose

This first legume run tests whether the frozen TandemX production discovery
profile yields a usable high-abundance tandem-repeat family catalogue from a
donor-matched T2T soybean dataset. It is a cost and candidate-chain pilot, not
a novelty result. The candidate sequences are not compared with known repeats
until de novo discovery is complete.

## Fixed sample and inputs

YSD56 is a wild soybean (*Glycine soja*) accession. NCBI records assembly
`GCA_040083835.1` and HiFi run `SRR28726931` under BioProject `PRJNA1095640`
and BioSample `SAMN40909152`. The assembly record reports a current Complete
Genome with 20 contigs and 1,008,523,555 bp. The source paper reports all 20
centromeres and 40 telomeres. The complete local HiFi FASTQ has already passed
ENA MD5, full gzip/FASTQ, read-count, base-count and duplicate-ID checks.

The discovery pilot uses the existing deterministic seed-6101 nested sample
`sample_003.fastq.gz`: 65,106 reads, 1,099,009,791 bp and nominal 1.08972149x
total-base coverage. Its SHA-256 is fixed in
`benchmarks/configs/legume_ysd56_novel_tr_pilot_v1.json`. The matched assembly
must pass the publisher MD5, gzip integrity and independent FASTA length/count
checks before localization.

## Frozen discovery command

Discovery uses the same production profile as the completed Ey15-2 and
Macadamia family catalogues: cascade discovery, sequence clustering at 0.95,
related-family audit, periods 30--1,000 bp, five period proposals, minimum five
supporting reads, minimum 100-bp repeat span, Rust backend and four threads.
The pilot reads the entirety of the fixed 1.0897x subset; it is not the
complete 44.193-Gb YSD56 HiFi collection, does not re-sample, and does not use
a known-repeat library during discovery.

## Pre-specified gates

The pilot passes the execution gate only if the command exits successfully,
all required outputs validate, source/configuration fingerprints are retained,
and resource usage is recorded. Its candidates may be used to test the
downstream family-ranking and novelty-exclusion workflow, but no candidate can
be called novel from this run.

If the pilot passes and the catalogue is not empty, the next discovery input is
the already generated 10.9570x nested sample. Formal candidates must then be
quantified against the complete 44.193-Gb HiFi file, localized on the matched
T2T assembly, and checked with matched ONT and Illumina evidence. Catalogue
similarity to the 1.0897x pilot is descriptive and cannot be treated as an
independent replicate because the read sets are nested.

## Candidate and novelty states

Family ranking will retain repeat period, supporting-read count, read-derived
abundance, assembly array span/distribution, related-family warnings, sequence
complexity and platform concordance. Candidates will receive one of these
states:

- `known_or_homologous`: sequence/period/location supports an existing named
  repeat family;
- `previously_unreported_candidate`: no qualifying match after documented
  species and cross-legume searches, with read and assembly support;
- `assembly_only_or_read_only`: one evidence class is missing;
- `unresolved_related_or_hor`: family hierarchy or period multiplicity is not
  resolved;
- `metadata_blocked`: provenance or donor matching is insufficient;
- `technical_failure`: an input or command failed validation.

Absence from a search result is not proof of novelty. A manuscript-level novel
TR requires a locked sequence, explicit search databases/versions and
thresholds, cross-platform support, assembly coordinates, exclusion of
telomeric/rDNA/organellar/simple-low-complexity artifacts, and preferably FISH
or another independent biological assay.
