# Macadamia native-read controlled-collapse development, 2026-09-18

## Why the source screen resumed

The 2026-09-17 Macadamia screen stopped before inspecting two complete HiFi
streams because T7 had previously shown an EIO/unmount event. The 2026-09-18
T7 read, gzip, sequential write/sync/readback and mount gate passed. USB2
remains a throughput limit, not a source-eligibility condition. The historical
negative screen remains a correct description of its earlier declared scope;
this is a new bounded development attempt, not a post hoc replacement of its
result.

The complete Macadamia assembly was 737,617,956 bytes, and the two original
HiFi FASTQ.gz streams were 6,716,963,362 and 11,433,466,175 bytes. All three
full SHA-256 values matched their prior receipts; both FASTQs decompressed to
completion (16,906,457,156 and 28,240,289,032 bytes). See
`docs/evidence/t7_stability_20260918/mac_source_{manifest,readback}.json`.
The two FASTQ streams contain 1,642,394 SRA records in total. SRA spot IDs
do not expose original PacBio ZMWs; distinct record IDs cannot be upgraded to
verified distinct original molecules. The paper-level same-sample relation
does not verify exact individual, DNA extraction or haplotype pairing.

## Selection, mapping and exact edits

The previously frozen assembly-only selection protocol was applied before
inspecting original reads. Seven of 28,368 earlier TandemX `locate` rows met
the length filter and two met the period filter. The selected interval was
`ctg.000105F:[118439,121594)` (0-based, 3,155 bp), family `TXF000219`,
operational period 144 bp; 3-kb natural flanks gave a 9,155-bp context.
The period-shift identity was 2,796/3,011 = 0.928595. The scored candidate
and exact context are in
`benchmarks/controlled_collapse/macadamia_assembly_candidate_v1_20260918/`.
Although this selection preceded **this** read mapping, the prior `locate`
catalogue and family ID came from earlier TandemX analysis of the source; it
is not wholly independent of all source reads or a held-out locus.

The context mapping protocol was frozen before mapping all 1,642,394 original
SRA records with pinned minimap2 2.31-r1302 (`map-hifi`, CIGAR, reported
secondaries up to 20). Of 86 primary geometric array spanners, seven
distinct records passed the >=99% whole-alignment identity and >=1-kb
natural-flank gates. The record-count gate (>=3) passed; a ZMW/molecule gate
could not be evaluated. Exactly those seven original records (five from
`SRR13557763`, two from `SRR13557762`) were re-extracted with full source
rehash and gzip/FASTQ checks. The selected combined FASTQ SHA-256 is
`da4c778e76c99c8489d0dfad429509cd1a64ec4a5c88051b92b5799fac7a9b70`.
Full 737-MB assembly mapping was separately frozen: all seven selected
records had one reported primary spanning the nominated locus at MAPQ 60,
identity 0.99173–0.99964 and no reported secondary under these settings.
This is bounded source anchoring, not proof of globally unique biological
origin. See `benchmarks/controlled_collapse/macadamia_native_mapping_v1_20260918/`.

The exact-edit generator made nine assembly contexts: one intact control,
terminal 25/50/75/100% deletions, internal 25/50/75%, and left-boundary 25%.
The minimum nonzero injected deletion was 788 bp. The independent verifier
reconstructed 9/9 edited sequences, confirmed unchanged original-read hashes,
and a clean replay produced byte-identical generated files. The only ground
truth is the **injected assembly-base difference**. Actual native repeat
copy number, monomer boundaries and biological collapse status are unknown.
See `benchmarks/controlled_collapse/macadamia_bp_provisional_v1/`.

## Frozen read-span result

The 5% length rule was frozen in `score_protocol.json` before scoring.
All seven original records had unique 1,024-bp natural-flank trims whose
boundaries agreed with primary PAF CIGAR projections within 10 bp. Their
spans were 3,155, 3,156, 3,156, 3,145, 3,157, 3,157 and 3,153 bp; median
3,156 bp. The scorer measured each edited FASTA directly before checking
its ledger. The fixed threshold was 158 bp at the original size. It marked
8/8 large injected deletions `DISCORDANT` and 1/1 intact control
`SUPPORTED`, with zero abstentions. Nine outcomes are correlated edits of
one chosen locus using the same seven records. This is a technical detection
check at these large edit magnitudes, not a biological sensitivity, FPR,
physical copy-number or cross-genome generalization estimate. Output:
`benchmarks/controlled_collapse/macadamia_bp_provisional_v1/span_score/`.

## Frozen M2 path result

The existing 144-bp read-derived `TXM000219` template is a supplied prior;
it is not independent of the source reads. Its sequence was checked against
the archived discovery monomers FASTA, and the unchanged M2 alignment
prototype and nine-case protocol were frozen before evaluation. Every one
of the seven original read intervals returned `AMBIGUOUS` with
`no_full_monomer_tiling`. Eight edited assembly intervals returned the same
reason; the fully deleted empty interval was formally `RESOLVED` but its
read comparison still abstained. Thus technical decision coverage is **0/9**
and all nine cases are `ABSTAIN`. These are refusals, not false negatives or
evidence of 0% biological accuracy. No monomer, parameter, trimming rule or
case denominator was changed after the result. Output:
`benchmarks/controlled_collapse/macadamia_bp_provisional_v1/m2_score/`.
The existing M2 synthetic Gate B NO-GO remains unchanged.

## Current scientific and hardware boundary

This Macadamia development set adds a second plant taxon for exact injected-bp
read-span feasibility but does not establish independent biological copy
truth, donor-held-out accuracy, or M1 absolute-abundance validity. The existing
production `quantify_primary` configuration, full-catalogue monomers and
`copy_number.tsv` were rehashed after the source-read requalification. The
`TXF000219` row is internally arithmetically consistent: k=21, 115 diagnostic
k-mers, 766.454 median depth / 28.9058 assumed haploid depth = 26.5156
estimated copies, or 3,818.2494 estimated family-wide bp. The row has
`medium` confidence, a broad 4.2282–62.5697 copy interval, zero empirical
single-copy controls and warnings about depth normalization, base-error
survival and unverified genome-background uniqueness. The 3,155-bp selected
array and seven spanners describe one locus; they cannot validate the
family-wide physical copy estimate. Exact hashes, bytes and the unchanged
row are in `benchmarks/controlled_collapse/macadamia_bp_provisional_v1/m1_existing_quantify_audit.json`.
M1 Gate A remains NO-GO; no new abundance
estimator was introduced. The mapping times
and memory recorded in source receipts are diagnostic run metadata only; no
formal TandemX runtime or scalability inference is drawn under USB2.

After full source readback, mapping and extraction, `/dev/disk5s2` remained
mounted at `/Volumes/T7` with about 1.9 TiB free. A one-hour macOS log query
for `I/O error`, `unexpected unmount` and `disk5s2` returned no matching
failure event. This scoped check cannot guarantee future storage stability.
