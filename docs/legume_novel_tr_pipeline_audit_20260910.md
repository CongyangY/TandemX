# Legume novel-tandem-repeat pipeline audit — 2026-09-10

## Scope and current status

This began as a local audit and future-run protocol.  It now records the
completed bounded YSD56 pilot and points to its retained execution receipts;
this document does not itself rerun TandemX. A read-discovered operational
repeat family is only a discovery candidate. The word `novel` is reserved for
a later evidence tier that passes the declared sequence-library, literature,
array-structure and independent-support gates; discovery alone does not
establish biological novelty, a validated HOR, physical copy number, or a new
algorithm.

The executable source and options for each result are governed by its retained
run configuration and source snapshot, rather than a fixed repository revision
in this audit. The repository's top status records newer comparator and
discovery-saturation work, but neither changes the real-data interpretation
boundary: assembly representation is a proxy, and a read--assembly deficit is
not physical missing-base truth. This protocol therefore retains a separately
gated assembly-aware formal run.

## Local inventory: legume material on T7

I searched `/Volumes/T7/Codex/TandemX` for Arachis/peanut, Medicago/alfalfa,
and other legume paths as well as FASTA/FASTQ/SRA/BAM/CRAM files.  No local
peanut or alfalfa sequence is enrolled.  However, a complete wild-soybean
HiFi WGS FASTQ and four hash-bound nested samples are already present on T7.
The donor-matched YSD56 T2T assembly is also local at
`/Volumes/T7/Codex/TandemX/data/references/YSD56_T2T_GCA_040083835.1_20260910/GCA_040083835.1_ASM4008383v1_genomic.fna.gz`.
Its provider MD5 record (`c8b216d1f7ba5cc4af20adfb14dc8bdc`) and the NCBI
dataset report are retained beside the FASTA. YSD56 is therefore the
immediate first-priority legume dataset for a separately gated assembly-aware
analysis; availability of this assembly does not itself validate a biological
copy-number claim.

| Audit directory | What is present | What it supports | Why it cannot start a run |
| --- | --- | --- | --- |
| `soybean_WM82_GCA030864155v1_20260906` | NCBI dataset report for *Glycine max* Williams 82 `GCA_030864155.1`: 20-contig, 1,011,793,938-bp assembly; the linked BioProject has no ENA read-run rows in the saved query. | Candidate assembly identity and source provenance only. | No assembly FASTA is stored; no enrolled raw read file exists. |
| `soybean_YSD56_PRJNA1095640_v1_20260906` plus `data/raw/SRR28726931_complete/` and `data/references/YSD56_T2T_GCA_040083835.1_20260910/` | ENA metadata and the complete *Glycine soja* YSD56 leaf HiFi WGS `SRR28726931`: 21,866,640,310 compressed bytes, 2,617,227 reads and 44,193,089,411 bases. ENA MD5, local SHA-256, gzip/FASTQ structure and duplicate-ID QC pass; the raw file is `SRR28726931_subreads.fastq.gz`. Nested sample FASTQs are available at 0.0116x, 0.1095x, 1.0897x and 10.9570x nominal total-base coverage. The donor-matched `GCA_040083835.1` / `ASM4008383v1` assembly FASTA is downloaded with its NCBI dataset report and provider MD5 record; the report identifies 20 contigs and 1,008,523,555 bp. | Immediate discovery and, after the next evidence gate, assembly-aware pilot input. The HiFi, ONT and Illumina WGS metadata share leaf BioSample `SAMN40909152`; the complete HiFi receipt, QC and sampling manifest bind the local reads. | The assembly is a representation reference, not independent physical copy-number truth. A frozen real-data normalization and analysis contract is still required before quantify/locate/compare. |
| `soybean_Jack_PRJNA701655_v1_20260906` | Three NovaSeq RNA-seq runs from soybean GFP-Jack seed-specific material. | Demonstrates an excluded source. | Transcriptomic RNA-seq is not a genomic read input for this protocol. |
| `soybean_ZH13_GSA_v1_20260906` | Two failed remote metadata fetch receipts. | Records that metadata retrieval failed. | No verified run, assembly, read object, or input hash exists. |

YSD56 source enrollment, raw-read integrity and matching-assembly acquisition
are complete. Its file QC and the matching donor record do not prove species
purity, ploidy, nuclear coverage, empirical read accuracy, repeat truth, or
that assembly representation equals physical copy number. Williams 82 and
YSD56 are different taxa/cultivars, so `GCA_030864155.1` must not be paired with
YSD56 reads for a read--assembly deficit result.

## Verified current workflow behavior

The runnable environment is the dedicated `tandemx-dev` Python 3.11.15
environment.  `tandemx --help` exposes `run`, `discover`, `quantify`, `locate`,
and `compare`; the `run` wrapper writes a root configuration, automatic-default
record, per-step logs, fingerprints, pipeline summaries, report products and an
output manifest.  It streams supplied FASTA/FASTQ inputs and rejects missing
paths before starting a run.

The actual `run` wrapper constructs the following steps:

| Step | Inputs and command behavior | Core outputs | Important boundary |
| --- | --- | --- | --- |
| discover | `tandemx discover --reads <reads...> --outdir discover --kmer-backend <resolved> --discovery-method <selected> --clustering-method <selected> --cluster-identity <selected> --family-audit <selected> --min-period <selected> --max-period <selected> --top-periods <selected> --threads <selected>`; optional read limits are passed through. | Candidate table/FASTA, `monomers.fa`, `families.tsv`, membership (when sequence clustering applies), similarities, hierarchy and audit summary. | The current general `run` defaults are `legacy`/`auto`, periods 2--2000, five peaks and automatic Rust/Python choice.  Do not describe those defaults as the frozen formal comparator profile. |
| quantify | Uses the discovery `monomers.fa`; invokes `tandemx quantify --reads ... --catalog ... --genome-size ... --outdir quantify --kmer-backend ...`, with optional haploid depth, controls, read-error and input limits. | `copy_number.tsv`, run config and log. | `--genome-size` is required.  If omitted from `run`, quantification is explicitly skipped unless an assembly is supplied, in which case its total length is recorded as a **provisional** denominator. |
| locate | Uses the enrolled assembly, discovery catalogue and copy-number table. | Repeat-density track, `arrays.bed`, compatibility comparison. | Localization is assembly representation, not independent copy truth. |
| compare | Uses `copy_number.tsv` and family-labelled `arrays.bed`. | `compare/assembly_vs_read_cn.tsv`. | The CLI defaults retain the frozen comparison cutoffs (0.6 possible-collapse and 1.5 possible-overexpansion); continuous deficit and binary status must remain distinct. |
| validate/report | `tandemx validate --project <outdir>` plus summaries, manifest and offline report. | Validation status, source-linked report and provenance. | A valid empty catalogue is an explicit outcome; stale downstream products must not be retained as positive evidence. |

The formal SRF protocol used a distinct, explicit controlled-input profile:
cascade discovery, sequence clustering, periods 30--1000 bp, minimum span 100
bp, support 1, Rust backend and one discovery thread; it quantified with k=21
and explicit sampled-read depth 1.  The latter normalizes a synthetic sampled
collection and must **not** be copied into a real-legume abundance analysis.
For a real run, record a genuinely justified genome-size/depth or independent
control normalization rather than inventing one from an unrelated assembly.

Existing large real-input evidence gives only host-specific orientation for
resource planning.  In the Macadamia donor-matched run, discovery took 959.24 s
with 1,447.00 MiB peak RSS and full-read quantification 11,292.61 s with 222.50
MiB RSS on roughly 22.55 Gb of reads.  Input type, catalogue size, filesystem,
backend and host contention differ, so these observations are neither a bound
nor a promised extrapolation.

## Enrollment gate before any future legume pilot

Create a fresh T7 output directory after the applicable conditions below are
met. YSD56 already meets the raw-read and matching-assembly enrollment
conditions, and its bounded discovery pilot has completed validation. The
remaining gate is an immutable real-data abundance/assembly-analysis contract,
not assembly acquisition.

1. Register one target species/cultivar and a primary scientific question:
   discovery only, or an assembly-aware read--assembly audit.  The latter
   requires a matching assembly and a documented material relationship to the
   raw reads.
2. Record authoritative accession, BioProject/BioSample/run IDs, library
   strategy, platform, tissue, extraction relationship, assembly accession and
   expected file bytes/checksums.  Exclude RNA-seq, Hi-C, Iso-Seq and
   unverified/partial transfers from the primary genomic read input.
3. For a new input, download only after the source record is frozen.  Verify provider MD5/bytes,
   local SHA-256, compression/record structure, read/base totals, duplicate IDs
   across input files and assembly sequence statistics.  A mismatch is
   `input_integrity_failure`, not a biological negative.
4. Freeze a source manifest containing input paths, hashes, command-line
   executable/version, Git commit, Python/environment identity, run config and
   the intended T7 output path.  Hash all source FASTA/FASTQ files, not merely
   a metadata table.
5. Predeclare whether abundance is out of scope for the pilot.  Do not allow an
   assembly from another accession/species to become a provisional normalization
   denominator merely because it is locally available.

## Two-stage legume novel-TR protocol

### Stage A: completed bounded YSD56 discovery pilot

Purpose: establish that the enrolled HiFi-like genomic input is usable and
inspect output states, without claiming abundance, assembly deficit, biological
novelty, or species-wide saturation.

* The existing nested YSD56 FASTQs, rather than a newly made prefix, were used.
  Their
  deterministic `blake2b_128_read_bernoulli_v1` selection (seed 6101), source
  SHA-256 and membership ID tables are already receipted.  `sample_001`
  (11,669,565 bp; 0.0116x) and `sample_002` (110,437,949 bp; 0.1095x) are
  input/output-state checks only. `sample_003` (1,099,009,791 bp; 1.0897x) was
  the bounded discovery pilot; it is nested within `sample_004`
  (11,050,418,089 bp; 10.9570x), which is the pre-existing escalation sample.
  These are independent read-inclusion samples, not independent specimens or
  fixed-base prefixes.
* The completed command used cascade discovery, sequence clustering, periods
  30--1000 bp, minimum span 100 bp, minimum support 5, Rust k-mer backend and
  four threads. It is a production-profile pilot, not a controlled-profile
  comparator.
* Discovery wrote 77,559 candidate reads and 1,474 operational families.
  `tandemx validate` then validated seven output files containing 241,099
  records. The profile receipt reports successful discovery in 435.15 s with
  1,963.88 MiB peak aggregate live-process-tree RSS, followed by successful
  validation in 1.82 s. These are execution observations, not claims about
  biological novelty, abundance or complete family recovery.
* The retained outputs include `candidate_reads.tsv`, `monomer_membership.tsv`,
  `families.tsv`, `monomers.fa`, `family_hierarchy.tsv`, family-audit output,
  logs, source snapshot and validation receipt. The 2,700 hierarchy edges
  remain candidate period-multiple or unresolved relationships; they do not
  establish HOR order or complete-array reconstruction.

The observed 1.099-Gb pilot resource receipt replaces generic planning for that
specific sample only. It must not be extrapolated directly to the 44.193-Gb
complete-read set: catalogue complexity, storage and host contention can
dominate. Any future run needs a declared wall-time/disk/resource stop rule;
an operator stop must be retained as an execution fate rather than silently
truncating output.

### Stage B: one registered formal run

The discovery and matching-assembly availability conditions are now met. The
next evidence gate is to freeze, before inspecting any abundance result, a
source manifest and Stage-B contract that specifies the real-data normalization
path, denominator/control rationale, family follow-up rule, resource stops and
the exact output directory. It must retain the assembly length as a recorded
representation proxy rather than treating it as physical copy-number truth.
Only after that contract is recorded may a new immutable outdir execute
discovery, quantification, localization on the matching assembly, comparison,
validation and report generation.

For the already local YSD56 complete HiFi amount (44.19 Gb), naive scaling of the cited
Macadamia observations gives about 31 minutes for discovery and about 6.2 hours
for quantification.  These are planning arithmetic, not performance estimates:
catalogue complexity, read lengths, storage and host state can dominate.  Stage
B therefore needs a reserved T7 scratch/output budget, serial execution unless
the resource plan says otherwise, a stated wall-time stop rule, and retained
native logs/partial outputs.  It must not begin merely because an assembly from
another soybean accession happens to be available.

## Family filtering and interpretation contract

1. Discovery receives raw enrolled genomic reads only; no known repeat library,
   hand-selected monomer, assembly sequence or truth label is used as a
   discovery input.
2. Retain every emitted family, membership state, warning and hierarchy edge in
   the primary tables.  Do not discard low-support, ambiguous or unresolved
   families after inspecting the result.
3. A separate follow-up subset may be defined only by a preregistered,
   machine-computable rule using primary outputs.  Its selection table must
   retain excluded IDs and reasons.  The audit intentionally does not invent a
   new support, abundance, identity or deficit threshold.
4. Optional post hoc annotation may describe similarity to an external known
   repeat library, with library version and hit criteria recorded.  It cannot
   rename a discovery family into validated homology or establish biological
   novelty.
5. `family_hierarchy.tsv` edges remain candidate period-multiple or unresolved
   relationships.  They do not validate HOR order, ancestry or complete-array
   reconstruction.
6. With a valid matching assembly, report estimated read abundance, assembly
   representation and their deficit/ratio with existing confidence/warning
   fields.  Do not call the difference physical missing bases or use a newer
   assembly as absolute truth.

## Required execution fates

Every planned stage needs a machine-readable state and an explanatory log.

| State | Meaning and handling |
| --- | --- |
| `not_enrolled` | Only metadata exists, or no matching read/assembly pair has passed source review.  No analysis may start. |
| `input_integrity_failure` | Provider/local checksum, bytes, decompression, FASTQ/FASTA structure, identifier or assembly-statistics check fails.  Quarantine the file; do not reuse a partial product. |
| `skipped_missing_normalization` | Discovery may be valid, but real-data quantification is not run without a recorded denominator/control path. |
| `skipped_missing_matching_assembly` | Discovery/possibly quantified catalogue exists; locate/compare remain unavailable because the assembly is absent or not material-matched. |
| `no_discovered_families` | Valid, schema-validated empty discovery result.  Downstream steps are explicitly skipped and no stale positive files are reused. |
| `technical_failure` | Nonzero exit, malformed output, validation/fingerprint mismatch, disk/resource exhaustion or an unrecorded interrupted process.  Preserve logs and partial files; do not score it as zero biology. |
| `operator_stop_resource_bound` | A predeclared wall-time/disk/memory safeguard stopped the run.  Preserve the command and observed resource state; it is neither a timeout nor a biological negative. |
| `completed_validated` | Inputs, outputs, fingerprints and required schema validations pass.  This state alone permits downstream reporting. |

## Decision boundary

YSD56 remains the first executed legume dataset. Its 1.0897x and nested
10.9570x discoveries both validated, and the retained 785-bp candidate recurred
exactly at the deeper depth. TRF and TideHunter independently confirmed 100%
coverage of its 124,029-bp assembly locus at dominant period 785 bp. This is a
stable, previously unreported-in-the-source-article candidate; it is not yet a
globally novel family because the deeper reads are nested, three named soybean
satellite sequences remain unresolved and no independent platform or full-
database exclusion has completed.

V14167 peanut and Medicago A17/R108 are the next registered lineages, but no
large payload was downloaded in this checkpoint. V14167's current small-pilot
plan is blocked by the repository contract that requires complete FASTQ
transfer and EOF validation before deterministic sampling. A17/R108 have an
eight-label repeat exclusion inventory but no directly reusable monomer
sequence records in the audited sources. The next scientific gate is therefore
independent support and remaining known-sequence exclusion for the YSD56
candidate, followed by a separately resourced complete-file acquisition for
one non-*Glycine* lineage. This is evidence collection, not new method
development.
