# Development evidence, not a finished manuscript

This directory intentionally versions compact executed benchmark tables and
editable figures. Large inputs, normalized predictions and logs are at the
configurable data root (`/Volumes/T7/Codex/TandemX` on the originating machine).

`cascade_native_screen_development/` combines two non-overlapping runs covering
all 16 development scenarios. The hash-checked archive shows no regression in
the common accuracy/error fields versus elastic, lower runtime in 16/16 cases
(geometric-mean ratio 0.562930), lower RSS in 15/16 (ratio 0.837814), and lower
boundary MAE in the 421-bp and 729-bp clean scenarios. Seed 1101 is development
evidence; it is not independent promotion evidence.

`cascade_native_screen_heldout_v1/` preserves the once-only seeds 3101--3103
matrix and its failed promotion decision. Of 432 planned executions, TandemX
and TideHunter completed 144/144 each and TRF completed 135/144; all nine TRF
failures were frozen 180-s timeouts on the low-complexity negative control.
TandemX passed the accuracy, determinism, negative-control, related-family and
RSS gates, but its paired wall-time geometric-mean ratio to TideHunter was
2.457943 against a 2.0 limit. Comparator completion was also mandatory. The
archive retains all failure commands, receipts and logs. `figures_v2` is the
accepted six-panel figure with editable SVG text; `figures_v1` is retained as a
draft because its text was converted to paths. Seeds 3101--3103 must not be
rerun or used for tuning.

`cascade_gap_free_development_v1/` preserves the fresh seed-1201 speed
development that followed the failed promotion. Its 29 hash-checked entries
include three identical 96-run matrices, a 1,600-read observable-feature audit
and native profiles. The baseline TandemX/TideHunter runtime ratio was 2.357698.
An unguarded candidate reached 1.933332 but worsened the 0.1%-indel boundary MAE
to 14.221 bp and was rejected. The selected guarded candidate reached 1.989148,
retained array/family endpoints, minimum base-union F1 0.997597, maximum
boundary MAE 2.921 bp and zero negative-control calls. This is development
evidence; seed 2201 was still unobserved when the rule was selected and frozen.

`cascade_gap_free_validation_v1/` preserves the one-time seed-2201 validation
run performed only after commit `054b935` and both hosted workflows passed.
All 96 commands and all 14 frozen gates passed. TandemX's minimum positive array
recall/precision were 1.0, minimum base-union F1 was 0.997445, maximum boundary
MAE was 2.35 bp and negative-control call rate was zero. Its
TandemX/TideHunter runtime and direct-child peak-RSS geometric-mean ratios were
1.977877 and 0.396864. Some condition-level runtime ratios exceeded 4, so the
result does not imply per-scenario or real-data dominance. The archive has 13
verified manifest entries. `figures_v2` is the accepted six-panel rendering
with editable SVG text; v1 is retained after clarifying its generic NA label.
The source receipt's revision warning records two untracked local native
extensions included and hashed in the execution snapshot; scoped tracked source
matched commit `054b935`.

`retrospective_collapse_source_audit/` freezes the exact ENA metadata queried
for four candidate old-to-new assembly comparisons and a 12-row PacBio genomic-
WGS selection. The metadata supports accession, platform and reported material
matching, but not identical DNA extraction, plant, stock or BioSample between
the read data and both assemblies. It is a source-eligibility audit, not a
completed biological collapse benchmark.

`tme204_donor_matched_source_audit_v1/` retains a negative eligibility result
for cassava TME204.  The study explicitly used CLR and HiFi reads from the same
DNA sample, but the 86-file GigaDB inventory and linked Mendeley record expose
the final HiFi haplotypes and supporting annotations, not the published
CLR-Falcon/Falcon-Unzip FASTA.  The candidate is therefore
`metadata_blocked_no_public_old_assembly`; no 56.14-GB read download,
reassembly, benchmark execution or accuracy claim was made.

`tomato_heinz1706_source_audit_v1/` retains a second negative source screen.
SL5.0, its exact HiFi run and the official 831,451,202-byte SL-T2T FASTA are
public, but the SL-T2T paper explicitly allows mixed-seed heterozygosity and
ONT-versus-HiFi sample differences. The pair is therefore
`not_source_eligible_donor_mismatch_risk`: same cultivar is not donor matching.
No approximately 28.77-GB HiFi download or benchmark execution was started.

`quantify_calibration_development_v1/` archives 108/108 successful public
`quantify` executions across three independently simulated genomes, 27 read
conditions and four methods. Empirical controls reduced aggregate MARE from
0.401898 to 0.365728 and improved each genome mean, while regressing at nominal
1x. The planted-error oracle reached 0.359640 but is unavailable in blind data.
Controls plus oracle were estimate-identical to controls and added runtime. The
post-hoc mean-control-depth >=2 candidate reached MARE 0.356117 and improved all
three development seeds, but remains unpromoted. The archive has 19 verified
manifest entries; `figures_v1` is the visually accepted six-panel rendering.

`quantify_calibration_fast_fasta_replay_v1/` verifies the committed FASTA
survival shortcut against the complete development matrix. All 5,940 metric
rows, 108 copy-number files, three control panels and non-resource summary fields
agree exactly. Full-matrix time changed 610.507→349.217 s and peak RSS
108.234→106.844 MiB. The 112 product hashes and 15-file compact archive are
verified. Resource changes are one same-machine replay, not publication timing.

`tidecluster_simulated_smoke/` records the first complete assembly-level
TideCluster 1.21.2 smoke test, its strict normalization, exact image ID,
container recipe, dependency versions, logs and internal GNU-time measurements.
All three planted families and arrays were recovered, but the test is only a
199.1-kb controlled assembly. Clustering used 7,669,232 kB maximum RSS inside
the emulated container. The archive must not be described as plant-scale
accuracy, portability or resource superiority evidence.

`tidecluster_morex_reference_scaling_v1/` contains hash-checked nested 10- and
100-Mb MorexV3 reference-window runs with the pinned TideCluster 1.21.2 image.
The compact archive retains all three provenance GFFs, cluster consensuses,
normalized arrays, finalization receipts and GNU-time stage resources. The
100-Mb output required family-membership-aware handling of 41 clipped and 131
merged final intervals; 1,200 intervals remained exact. Its six-panel SVG has
68 editable text nodes and no raster nodes, and its PDF was independently
rendered and visually checked. These are descriptive reference calls without
independent family/array truth or whole-chromosome context.

`tidecluster_1gb_resource_gate_v1/` preserves the frozen MorexV3 1-Gb
preflight refusal as a formal resource fate. The successful 100-Mb clustering
stage used 7,957,438,464 bytes, 96.836% of the same 8,217,432,064-byte Docker
limit, exceeding the predeclared 85% threshold. The 1-Gb command was not
started, so runtime and accuracy are unavailable rather than zero. This is a
same-host, same-container infeasibility result, not a universal TideCluster
memory limit.

`abundance_classifier_depth_gated_validation_v1/` preserves the frozen v3
development-to-held-out chain. Commit `62892a6` and dual-platform hosted CI
preceded one-time use of seeds 5801-5803. The rule passed all predeclared gates:
sensitivity increased 0.594650 to 0.653635 with unchanged FPR 0.016461 and
precision increased 0.981880 to 0.983488. This is known-catalogue IID simulation,
not biological validation or a universal superiority result. Figure 7 v2 is the
accepted six-panel rendering; v1 is retained as an inspected draft. Evidence
commit `3b21825` passed hosted Ubuntu/macOS Python, Rust and wheel runs
`34047548628`/`34047562996` on the work branch and main.

`discovery_packed_trace_batch_v1/` records two completed exact-output Mo17
discovery replays after packing native traceback directions and batching all
periods from one read through one native call. At 11.681/111.506 Mb, elapsed
time fell 28.75%/25.72% and peak RSS fell 23.52%/13.27%; six/seven stable
products were byte-identical. The compact archive retains the dirty-source
warning, exact source digest and file hashes. Each row is one historical
baseline plus one replay, so it is engineering evidence rather than a repeated
publication timing claim or universal external-tool superiority result. Source
commit `a73398d` passed hosted Ubuntu/macOS Python, Rust and wheel runs
`34048998182`/`34049010825` on the work branch and main.

`challenge_v1_baseline/` records the discovery algorithm at Git revision
`08e100d` using the initial challenge harness hashes in `environment.json`.
The complete data and command logs are in `results/challenge_v1_development_20260906`
at the data root. 144 attempts produced 135 successes and 9 explicit TandemX
errors on negative controls. Retain failures as failures after the CLI fix.
Timing is exploratory because other development activity overlapped.

`negative_control_fix/` records nine successful post-fix TandemX runs with zero
calls on the same three negative datasets. It tests zero-result handling only.

## Elastic development checkpoint

- `elastic_development_initial/`: first elastic attempt, 16/16 successful
  executions, but three AT-rich negative reads were called and one divergent
  array consensus had an incorrect length. These are retained failures of the
  method, not discarded inconvenient datasets.
- `elastic_development_corrected/`: same inputs after composition correction
  and support-based consensus-template selection. All 13 positive scenarios had
  array recall/precision and strict family recovery 1.0; all 300 negative reads
  yielded zero calls. This is one development seed with one run per scenario.
- `elastic_validation/`: seed 2101, 48/48 runs completed with matching output
  digests across three repetitions. Array recall/precision were 1.0 throughout,
  but related-family recovery was **2/3**, because two related monomers were
  merged. Other family-recovery scenarios scored 1.0; negative controls had no
  calls. These results must accompany, not be hidden by, the development figure.
- `ultra_pilot/`: ULTRA 1.2.2, pinned source/build provenance, two ten-read
  feasibility scenarios. This is a different input size from the main challenge.
  It is not a head-to-head resource comparison. Strict equal-length family
  recovery can fail on a one-base consensus-length error; a complementary
  independent indel-aware homology endpoint is required before paper claims.
- `ultra_tuned_pilot/`: the same ten-read inputs with upstream `--tune
  --tune_indel`. Both scenarios reached array recall/precision 1.0. Total times
  were 166.58 and 176.81 s, including 18 parameter settings and shuffled-input
  checks; these are not final annotation-only timings. The tuned competitor's
  recovery is retained alongside its less accurate default pilot. Strict
  equal-length family recovery was 0.0 for both; inspect the actual consensuses
  and add a gapped homology metric before interpreting that as missing families.

`indel_detector_gap_audit_v1/` reuses the consumed seed-2201 outputs without a
new run or threshold change. Both TandemX and TideHunter reached array recall
1.0 at 0.1%, 1% and 4% total indels; TandemX precision remained 1.0 while its
median elapsed time was 3.34--3.80 times TideHunter's. The audit therefore
defers a seed-and-chain replacement because the observed deficit is runtime,
not frozen synthetic read-local accuracy. It does not settle real-read truth or
generalize beyond the simulated error distribution.

These early runs preserve full source/executable hashes but not the complete
intermediate dirty source. They are diagnostics, not frozen release benchmarks.
The updated runner snapshots source/build inputs and loads the TandemX snapshot;
the final publication experiments must use committed sources and fresh runs.

The corrected four-panel figure uses identical input-file hashes for legacy,
TRF, TideHunter and elastic. Panels a–c use the definitions below; panel d adds
the recovered elastic interval on the same original failure example. Baseline
rows have three timing repetitions and elastic rows have one. Resource timings
are not shown because concurrent development makes them exploratory. The
validation failure above is not part of this single-seed development heatmap.

## Four-panel caption

**Development challenge benchmark.** a, One-to-one array recall on 13 positive
scenarios, requiring the same read, IoU >=0.5 and period error <=max(2 bp, rounded
2% of truth period). b, Raw array precision; unmatched overlapping or harmonic
calls count against precision. c, Family recovery at equal monomer length and
>=90% exhaustive circular ungapped identity on either strand. A single distinct
phase/strand-canonical consensus cannot recover two truth families. TandemX
uses its final catalog; per-array finders use emitted consensuses. d, Planted
and predicted intervals on the first deterministic indel-case mismatch,
illustrating incomplete boundaries despite repeat detection. All panels use
development seed 1101. Each positive scenario has 100 reads and 70 positive
reads; the two-array scenario has 140 arrays. Three timing repetitions reuse the
same inputs and are not independent biological replicates. Negative controls
are retained in tables but excluded from the heatmaps because positive recall
has no denominator. These results motivate improvement, not final superiority.

Generated by `benchmarks/scripts/plot_challenge_benchmark.py`. SVG check: zero
raster images, 182 editable text elements. Full PNG visually inspected. Source
tables and figure-input SHA-256 digests accompany the figure.

## Sequence-clustering validation and independent endpoint audit

`sequence_validation/` archives the fd656b6 validation run: 16 scenarios, seed
2101, three identical-output repetitions, 48 successes. All 13 positive scenarios
have array recall/precision and cyclic planted-monomer recall 1.0; each of three
100-read negative datasets has zero calls. Related-monomer recovery increased
from 2/3 to 3/3. The three new observed representatives have support of 32, 26
and 12 reads and match the three planted sequences exactly. This reused seed
informed development and is not untouched validation evidence.

`independent_sequence_audit/` re-scores first-repetition archived outputs using
edlib 1.3.9.post1 and per-read unions of intervals. It preserves original failed
rows. An equal-length-only score wrongly suggests absent clean ULTRA monomers;
the added cyclic edit score recovers them. TRF's duplicated calls can lower raw
call precision while base coverage precision stays high. Both endpoints remain
visible; these findings must accompany comparisons, not be hidden in supplements.

The four-panel `sequence_validation/figures/sequence_audit.pdf` and editable SVG
show all four ULTRA pilot conditions, all 13 positive TRF development scenarios,
a 90–100% clustering-resolution sweep of a fixed candidate set, and independent
similarities of old/new catalogues to planted units. Panels C/D use identical
seed-2101 input hashes. Replay uses rounded exported candidate scores and verifies
that 95% reproduces the actual catalogue. There is no claim that 95% defines
biological family ancestry. `panel_source.tsv` and `figure_provenance.json` link
every panel to inputs. The final SVG has 95 editable text elements and no images.

`comparator_builds/` records SRF, KMC, minimap2 and k8 provenance, including KMC's
small libc++ compatibility patch, failed build attempts, and an independent exact
k-mer count validation. These build records alone are not benchmark evidence.

`unitfinder_source_scope_v1/` fixes the official source commit, publisher
supplement identity and the evidence boundary for unitFinder.  The 60 Table S1
soybean intervals are published workflow outputs that can test reproducibility,
not independent biological truth.  The audit keeps the failed PMC workbook
download and makes CENH3 ChIP-seq the required route to an accuracy endpoint;
it reports no successful unitFinder execution yet.

## SRF workflow evidence

`srf_pilot/` retains the initial four successful native workflows and the exact
original runner (restored and hash verified). `srf_initial_development/` retains
all 32 expanded attempts, including 11 native assertions on empty KMC dumps.
`srf_guarded_development/` records the full guarded rerun: 17 ordinary successes,
four native no-catalogue outcomes, 11 observed-empty-count skips and no process
failures. Empty-count skips create no fake native FASTA. Both summaries remain
available; do not retroactively relabel original crashes as successful tool runs.

At ci20, ten of thirteen positive scenarios recover all three planted monomers;
three high-mutation scenarios have no eligible k=151 counts. All negative controls
have zero in-scope calls. These are conditional development results with two
count settings, no technical repetitions and no tuned k-mer-length sweep. A
higher-order unit is not automatically a missing constituent monomer. SRF's
native abundance estimates are based on mapped sequence, not validated biological
copy number. See `docs/srf_workflow.md` for all workflow and build caveats.

## Native seed processing: paired development performance

`native_seed_paired/` retains 96 successful executions of frozen old/new sources
on 16 identical seed-2101 datasets, three repetitions each and shuffled pair
order. All six discovery data files are byte-identical. Old code uses Python
seed processing; new code is 3f2482f. The per-dataset speedup range is 1.200–4.447
(median 1.371); positive scenarios range 1.200–1.657. RSS changes range -13.09%
to +2.03%, with six increases. This is not a claim of an all-metric gain or an
external competitor comparison. Reused synthetic inputs and technical timing
repetitions do not establish generalization or biological uncertainty.

Four-panel caption: A, median wall time for both implementations in all 16
scenarios. B, ratio of baseline to native median wall time. C, percent change
in median peak RSS, including increases (orange). D, median user CPU time.
Each input has 100 reads and three repeated executions per implementation;
positive scenarios have 70 positive reads (two-array case: 140 arrays).
Scientific parameters and source snapshots are fixed; six data-file hash sets
match in every pair. Raw CPU system time remains in the source table. All
figures preserve editable text; the PNG was visually inspected.

## Abundance baseline: failures are retained

`abundance_baseline/` records all conditional CN/localization/comparison scores
from 8218beb: 177 successful commands, three independent 199.1-kb genomes, nine
coverage/error conditions and five retained-copy assemblies. Known monomers are
supplied, so this does not measure discovery. Only 6/81 diagnostic-k-mer spread
intervals contain genomic truth; those intervals are not sampling confidence
intervals. At 20× with 1% substitutions, mean signed CN error is -18.31% and
estimator-minus-sampling-oracle is -20.08%. All raw conditions, including false
calls/misses and 1× sampling fluctuations, are retained. See
`docs/abundance_benchmark.md`; this baseline is not scientifically accepted.

`abundance_heldout/` records the first and only execution of the predeclared
5101–5103 conditional seeds on source f16596b. All 177 commands completed. Across
243 positive and 162 control family conditions, TP/FN/FP/TN were 208/35/8/154;
sensitivity was 0.855967, false-positive rate 0.049383 and precision 0.962963.
All positive exact-copy localization rows had base recall 1, but 20×/1%-error
50%-retention sensitivity was only 2/9. The compact archive revalidates source,
configuration, matrix dimensions and all receipts; these seeds are now consumed.

`abundance_heldout_v2/` records the fresh 5201–5203 baseline that was configured,
committed and pushed before its first execution. All 177 commands completed. The
original k=21 rule had TP/FN/FP/TN=198/45/14/148. The paired
`abundance_multik_collapse_heldout/` archive applies the frozen development rule
without held-out fitting and changes these counts to 208/35/12/150: sensitivity
0.814815→0.855967, false-positive rate 0.086420→0.074074 and precision
0.933962→0.945455. It preserves 15 unavailable multi-k-only rows, seed-5202
false positives and all exact calibration/baseline/source hashes. This remains
small known-catalogue exact-copy simulation evidence, not biological validation.

`read_cluster_replay/` adds the non-default read-sampling reference model replay
from 6f29db5: 81 point estimates, 33 available intervals, 24 covering truth and
48 unavailable. Mutation-dependent CN underestimation remains. This is not
calibrated real-data inference or a new default quantifier.

`Mo17_input_qc/` contains complete CCS-file QC, seeded nested sampling receipts,
raw histogram source data, a four-panel editable figure and verified metadata
for the 10-chromosome Mo17 T2T reference. Exact study/cultivar context does not
establish identical DNA donors. One included technical batch and its nested
samples are not independent biological replicates or the whole published dataset.

`Mo17_real_pilot/` preserves three completed three-tool comparisons, including
the original slow audit and exact six-file native parity. The larger input is
111.506 Mb, not a completed whole-library or plant-genome benchmark.
`Col0N_input_qc/` contains the complete 14.647-Gb input QC. `Mo17_published_regions/`
contains the original paper supplement and 96 extracted regions, with coordinate
compatibility checks and the explicit mixed-region/base-truth distinction.

`Nipponbare_input_qc/`, `Morex_input_qc/`, `Lo7_input_qc/`,
`ChineseSpring_input_qc/`, `Victoria_input_qc/`, `Col0R_input_qc/` and
`Ey15R_input_qc/` extend complete-file QC and seed6101 nested sampling.
`YSD56_input_qc/` adds exact acquisition, complete-file QC and a completed
seed6101 sampling ladder. Together with Mo17 and Col-0N, the included cohort is ten
libraries from eight reported species (14,937,608 reads; 262,731,255,175 bp).
These are validated files and technical sampling strata, not independent
biological replicates or accuracy truth. Large FASTQ/ID products remain at the
recorded T7 data root.

`MorexV3_reference_qc/` archives the original IPK source page, pinned plan,
completion receipt and a manifest that rechecks the external 4.30-GB FASTA.
Eight unique FASTA records contain 4,225,605,719 bases and the file matches the
published SHA-256. The FASTA remains outside Git. This is file integrity and
sequence-content QC, not identical-donor or satellite copy-number truth.

`Morex_115Mb_index_interface_v2/` contains the compact clean-source replay of
the Python word bridge and sequence-native Rust interface on 7,094 fixed Morex
candidates. Complete 4,380-family payloads are byte-identical. The one-run time
and RSS reductions are favourable engineering diagnostics; fixed order and
native build-artifact provenance prevent a publication timing claim.

`factorial_multik_replay/` and `factorial_joint_multik/` retain the complete
three-genome, 27-condition conditional copy-estimation and joint-read interval
experiments, including missing fits, intervals without truth and conditions in
which the multi-k estimate is worse. Version 2 of the joint-read figure changes
layout only and preserves the version-1 source rows.

`multispecies_real_diagnostics/` contains compact one-thread TandemX/TRF/
TideHunter execution receipts from whole-library random samples. It preserves
calls and called-base extent as descriptive outputs; no real accuracy truth or
final isolated resource ranking is implied. `reference_mapping_diagnostics/`
contains the separate Col-CEN and Nipponbare alignment audits with reference-
content, donor-matching and organellar-denominator limits.

`multispecies_input_qc/figures_v2` contains the current inspected six-panel
cross-cohort input figure, exact panel source and output/input hashes for all ten
included libraries. It visualizes complete-file data volume, read length, GC and
reported quality; it does not measure empirical read accuracy or biological
replication. The nine-library v1 render remains as history.

`known_repeat_sources/` records exact source accessions and curated historical
repeat queries. `known_query_recovery/` scores those selected sequences against
actual tool consensus outputs, including the Morex TandemX family-representation
miss at the 0.90 threshold. It is source-query recovery rather than donor-matched
genome-wide family recall or prediction precision.
