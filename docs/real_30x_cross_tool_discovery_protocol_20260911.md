# Real 30x cross-tool tandem-repeat discovery protocol

Status: frozen before inspection of any new 30x tool output on 2026-09-11.

## Question

This analysis asks whether TandemX recovers repeat families from real HiFi reads
that are not recovered by other applicable de novo read-level tools under a
common input and normalization protocol.  It does not ask whether another tool
can detect the same tandem array after its locus is supplied.

The phrase `TandemX-only` is therefore reserved for
`TandemX_only_de_novo_under_fixed_30x_workflows`.  It is not a claim that the
biological repeat is detectable only by TandemX.

## Cohort and eligibility

The primary cohort is the six peanut assemblies reported as T2T: V14167,
K30076, S245, HN873, HN51 and S83.  YSD56 and a donor-matched high-completeness
cultivated soybean are the first small-genome controls.  Additional small
genomes enter only after their assembly, raw HiFi accession, material identity
and public release state are closed in an enrollment receipt.

*Brachypodium distachyon* Bd21 was withdrawn from the cohort by the user on
2026-09-11. Its only applicable public HiFi endpoint is a
517,292,608,950-byte BAM without a repository checksum, while the assembly is
near-complete rather than T2T. No download, BAM conversion, 30x sampling or
cross-tool analysis will be attempted. The state
`user_withdrawn_large_public_bam_low_operational_value` is a scope decision and
must not be interpreted as a failed biological or software comparison.

## Common input

1. The discovery input is PacBio HiFi FASTQ from the enrolled material.
2. Nominal input depth is 30x, using the verified assembly span as the
   denominator.  The receipt reports exact sampled bases and observed depth;
   30x is not treated as an exact biological coverage guarantee.
3. A source file must pass expected byte size, repository checksum when
   supplied, gzip/FASTQ EOF parsing, unique read-ID checks and total-base
   reconciliation before sampling.
4. Sampling is deterministic and source-order preserving.  Every tool receives
   the byte-identical normalized FASTA generated from the same FASTQ sample.
   The 30x source is split into three disjoint, source-order-preserving
   approximately 10x partitions using the verified assembly span.  The runner
   records each observed partition depth and refuses a formal run unless all
   three partitions are within 1x of the 10x target.
5. NCBI, ENA and NGDC transfers use direct no-proxy access by default. After
   the user's 2026-09-11 update, a bounded proxy A/B may consume at most 0.5 GB
   before a long-transfer decision. A long transfer may switch only when the
   proxy route is verifiable, provides at least a stable twofold speedup and
   the estimated project proxy traffic fits the current 50-GB project cap. The
   acquisition receipt records the URL, route, transferred bytes, checksum and
   downloader. The user's approximate 350-GB remaining monthly allowance is
   not treated as permission to consume the full balance.
6. Failed, incomplete or zero-byte transfers never become analysis input.

## Tools and endpoints

The common read-level discovery endpoint includes the frozen TandemX
production discovery method, TideHunter 1.5.5 and TRF 4.10.0-rc.2.  Commands,
binary/source hashes, thread counts, chunk manifests, wall time, CPU time, peak
process-tree RSS, exit state and output hashes are retained.

For TandemX, each partition uses the frozen YSD56 discovery settings:
`cascade`, sequence clustering at 0.95, `family_audit=related`, periods
30--1000 bp, `top_periods=5`, `min_support_reads=5`, 100-bp minimum span and
the Rust k-mer backend.  This comparison does not substitute a permissive
one-read-support pilot command.

TRASH/TRASH2 and TideCluster are assembly-oriented endpoints and are used only
where their native design is applicable.  SRF may be evaluated on a separately
declared donor-matched short-read endpoint.  A tool without the relevant input
or endpoint is `N/A`; a crash or timeout is `technical_unresolved`, never zero
recovery.

## Family normalization

Native arrays are filtered to the shared 30--1000-bp period and at least
100-bp-array scope.  Their reported monomers are canonicalized across reverse
complement and circular rotation.  Native IDs are not compared directly.

Two families are considered recovered by both workflows when the complete
shorter representative has at least 0.90 glocal identity and a direct
`shorter_length/longer_length >= 0.90`, or at least 0.80 identity with an integer period-multiple
relationship whose relative length error is at most 0.05.  The permissive
period-multiple rule is deliberately conservative against declaring a
TandemX-only family.

## Candidate states

A TandemX family can enter the de novo-only shortlist only when:

1. the TandemX run completed and the family passes the frozen production
   confidence/support rules;
2. every required read-level comparator completed and normalized successfully;
3. no comparator representative meets the common-family rule above;
4. the family is recovered from at least two disjoint approximately 10x
   partitions of the same 30x library, or receives independent platform/sample
   support;
5. an assembly is available and the family localizes to a tandem array with
   unambiguous coordinates; and
6. rDNA, organelle, known satellite/minisatellite and available public sequence
   screens are reported with explicit database scope.

If targeted TRF or TideHunter confirms the localized array, the candidate may
remain `TandemX_only_de_novo_under_fixed_30x_workflows`, but it is explicitly
`not_tool-exclusive_after_targeted_locus_analysis`.

Public-database local matches are reported by aligned query interval, identity
and coverage.  Absence of a near-full-length match is not proof that the
sequence, family or biological function is novel.

## Stopping and reporting rules

- No exclusivity conclusion is made when any required comparator fails,
  times out or lacks a valid normalized catalogue.
- Operational representatives are not reported as biological family counts.
- A new combination of known sequence blocks is labelled a composite or
  previously-unreported unit candidate until independent evidence supports a
  stronger description.
- New real-data results remain outside the manuscript until the user decides
  whether a sixth Result should be added.
- Bioconda, Zenodo and a formal tagged release remain paused.
