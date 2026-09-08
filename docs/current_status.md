# TandemX current status and handoff

Updated 2026-09-08. Read completely after `AGENTS.md`, then verify Git/tests.
The user has authorized autonomous development and GitHub updates toward mature
software and a full evidence-backed paper (multi-panel figures and supplement).
The goal is **not complete**. Acceptance gates: `docs/release_program.md`.

## Active 2026-09-08 continuation checkpoint

This section is the restart point for a new Codex window. Verify the state below
before continuing and do not stage the untracked `.codex/` directory.

- **A Tier A Ey15-2 donor-matched collapse candidate is now preregistered.**
  Rabanal et al. explicitly compare CLR and HiFi assemblies of the same Ey15-2
  sample (9994/CS76399), correcting the earlier audit that considered Col-0 but
  missed this same-sample contrast. The primary pair holds the Bionano
  scaffolding context constant: `9994.CLR_Canu` versus `9994.HiFi_Hifiasm`;
  the final HiFi-Hifiasm+CLR-Canu assembly is sensitivity-only. The source and
  complete ERR8666125 read library are already fixed; the 1.605-GB official
  Zenodo bundle is still being acquired and must pass its published byte count
  and MD5 before use. The frozen v1 rules are in
  `benchmarks/configs/ey15_donor_matched_collapse_v1.json`; do not change them
  after inspecting old/new localization results.
- **The donor-matched evaluator is frozen before result inspection.** Primary
  source eligibility depends only on at least 15 kb localized in the newer
  assembly, not on new/read agreement. It preserves all family fates, separates
  `not_source_eligible` from `technical_failure`, reports predeclared 5/15/50-kb
  denominator sensitivities and computes confusion, Wilson intervals, balanced
  accuracy, MCC and missing-bp agreement. Five focused tests and the full suite
  pass: **566 tests in 64.12 s**. This is a design/software checkpoint, not a
  completed biological result. The newer assembly shares HiFi evidence with
  the predictor and remains a donor-matched high-quality reference proxy rather
  than absolute independent truth.

- **TideCluster comparator is active and must not be restarted blindly.** The
  pinned `tandemx/tidecluster:1.21.2` image has immutable ID
  `sha256:62691b116427394984a8bace59a6f4ba7f3ea834373772734fbcfc567c3a27a7`
  and contains the required TideHunter 1.4.3, KiteHOR 0.13.2, MMseqs2 commit
  747c64 and BLAST 2.16.0+. Deterministic nested MorexV3 samples are at
  `/Volumes/T7/Codex/TandemX/results/MorexV3_reference_windows_s8101_v2_20260907`:
  10/100/1,000 windows and 10/100/1,000 Mb, with SHA-256 values recorded in its
  receipt. The v1 sample has identical sequence records but unsafe `=`-bearing
  identifiers and is retained as a failed-interface case.
- **The 10-Mb real-reference TideCluster gate is complete.** Do not rerun
  `/Volumes/T7/Codex/TandemX/results/MorexV3_tidecluster_docker_10mb_s8101_v4_20260907`.
  Its two external stages succeeded and the corrected three-GFF provenance join
  was applied without rerunning them. It reports 87 arrays, 28 operational
  families, 239,939 union bp (0.0239939 of the sampled bases), and positives in
  10/10 windows. TideHunter used 25.03 s/561,556 kB internal GNU-time maximum
  RSS; clustering used 45.92 s/7,809,052 kB. Copy number is retained for 85
  exact TideHunter intervals and explicitly unavailable for two merged/resolved
  intervals. All 20 files in `finalization_receipt.json` independently passed
  byte-size and SHA-256 rechecks. These are descriptive calls on real reference
  windows; no independent array/family accuracy truth exists.
- **The 100-Mb TideCluster resource gate is complete.** Do not rerun the two
  external stages in
  `/Volumes/T7/Codex/TandemX/results/MorexV3_tidecluster_docker_100mb_s8101_v1_20260907`.
  The first normalizer stopped after both stages because final intervals can be
  clipped or merged relative to intermediate GFF rows. A family-consensus map,
  overlap-coverage check and deterministic representative selection resolved
  all 1,372 final intervals without rerunning TideCluster: 1,200 exact, 41
  clipped and 131 merged. The run reports 130 families, 2,356,743 union bp
  (0.02356743) and positives in 99/100 windows. TideHunter used 331.52 s/
  5,152,444 kB; clustering used 60.99 s/7,770,936 kB. Copy number is unavailable
  for 202 merged/resolved intervals. All 20 finalization files independently
  passed hash/size checks. The 1-Gb sample remains unrun because the 100-Mb
  memory gate does not yet justify that resource risk.
- **TideCluster scaling evidence is compact and visualized.** The 10/100-Mb
  archive is at `paper/evidence/tidecluster_morex_reference_scaling_v1`; all 32
  archive-manifest entries pass independent hash/size verification. Its accepted
  six-panel SVG/PDF/PNG has 68 editable SVG text nodes and no raster nodes. Both
  the direct PNG and independently rendered PDF were visually inspected. Calls
  are descriptive real-reference output, not accuracy or whole-genome evidence.
- **`tandemx cohort` was hardened concurrently.** Commit `50de982` rejects
  unknown, duplicate and incompletely quantified catalogue
  families; records per-sample input counts and SHA-256 values; propagates the
  least local confidence; and adds wide lower/upper abundance-endpoint matrices.
  The current continuation adds an automatic four-panel cohort SVG/PDF, stable
  abundance-ranked `--top-families` selection, machine-readable plot source and
  a hash/node-count figure receipt. Its toy workflow completed two end-to-end
  samples and validated 11 files/30 records at
  `/tmp/tandemx-toy-cohort-visual-v2`; an independent PDF render passed visual
  inspection, and the SVG has 59 editable text nodes and no raster image nodes.
  The complete source suite passes 548 tests in 66.69 s; compileall, 72 tracked
  local Markdown links and diff checks pass.
- **Current verification after TideCluster evidence.** The full Python suite
  passes 547 tests in 62.84 s; compileall and diff checks pass. The archive has
  32/32 independently rehashed entries. Supplementary Figure S5 has six panels,
  68 editable SVG text nodes and no raster image elements; its direct PNG and a
  Poppler render of the PDF were both visually inspected with no overlap or
  clipping observed.
- **Published cohort-visualization checkpoint.** Commit
  `9723aad7e031dabd5c2c6a29f3273a2bd3e6d51d` is on both
  `codex/publish-current-progress` and `main`. Hosted runs `34105542000` and
  `34105542475` passed Ubuntu/macOS Python tests, Rust checks and wheel builds.
  Commit `50de982` is the cohort hardening/TideCluster runner checkpoint;
  `c8593e0` archives the completed Morex 10/100-Mb evidence. The earlier
  commit `9eda196` froze
  the untouched depth-gated quantify validation source, seeds, data-generation
  hashes, decision rule and nine scientific gates before any validation data
  were generated; `9291fd6` publishes the completed evidence and opt-in CLI.
- **External TideHunter calls now have a first-class downstream route.** The
  current continuation adds `tandemx import tidehunter`: strict 11-field `-f 2`
  parsing, disk-backed validation against the original reads, full/native ID
  reconciliation, 1-based-inclusive to 0-based-half-open conversion, native
  field audit, fixed-representative family clustering and standard catalogue
  outputs for quantify/locate/compare/probe/cohort. The first real toy attempt
  at `/tmp/tandemx-tidehunter-import-v1` correctly failed because TideHunter
  retained semicolon header metadata while TandemX normalized it; the corrected
  mapping then completed on fresh `/tmp/tandemx-tidehunter-import-v2` using
  TideHunter 1.5.5: 37 imported candidates, one family, seven validated files/
  150 records and a successful downstream quantify run. External provenance is
  explicit and is not counted as TandemX native detector evidence. Nine import
  output hashes independently rechecked, and the complete source suite passes
  556 tests in 64.55 s; compileall, shell syntax, 72 tracked local links and diff
  checks pass.
- **The TideHunter importer checkpoint is published and portable.** Commit
  `20aca270a4c51ccedfb59a917c4633582a3ef8cf` is on both
  `codex/publish-current-progress` and `main`; hosted runs `34107440725` and
  `34107440959` completed successfully. An independently built CPython 3.11
  macOS ARM64 wheel includes `tandemx/importers/tidehunter.py`; a temporary
  target installation outside the source tree imported the module and exposed
  `tandemx import tidehunter --help`.
- **Candidate monomer/HOR architecture is published and machine-readable.**
  Commit `9f9196498e0e3b14ccf86e6daa68335ff19aedf2` is on both branches;
  hosted runs `34109053190` and `34109049735` completed successfully. Every
  native discovery and TideHunter import writes `family_hierarchy.tsv` while
  streaming the existing pair audit. Each
  `possible_higher_order_or_partial` pair is directed shorter-to-longer and
  classified as `putative_period_multiple` only when the representative-length
  ratio is within 0.05 of an integer of at least two; other related edges remain
  unresolved. All alternatives are retained, so the output is an evidence graph
  rather than a forced tree or validated HOR call. An actual 171/342/684 bp
  import at `/tmp/tandemx-hierarchy-import-v1` produced three families and all
  three pairwise candidate edges; eight files/36 records validated and 10 import
  hashes independently rechecked. The current source suite passes 561 tests in
  65.59 s. A release wheel installed outside the source tree, repeated the
  import with the Rust backend, validated the same eight files/36 records and
  has SHA-256
  `c6c44bfc43fdf037848d3aae2ae9034f5e2f8d3ce19c6784680c78054b4bc504`.
- **Frozen cascade held-out is consumed; never rerun seeds 3101--3103.** The
  once-only run is complete at
  `/Volumes/T7/Codex/TandemX/results/cascade_native_screen_heldout_v1_20260907`:
  432/432 planned rows, 423 successes and nine explicit TRF timeouts on the
  low-complexity negative control. TandemX and TideHunter completed 144/144;
  TRF completed 135/144. Failed processes remain missing measurements and must
  never be converted to zero accuracy.
- **Cascade promotion failed 2/12 predeclared gates.** TandemX itself had no
  process failures, deterministic output fraction 1.0, minimum positive array
  recall/precision 0.985714, zero negative-read calls, related-family cyclic
  recall 1.0 and direct-child RSS ratio 0.455130 versus TideHunter. Promotion
  failed because all comparator runs were required to complete (nine TRF
  timeouts) and TandemX/TideHunter wall-time geometric-mean ratio was 2.457943,
  above the frozen 2.0 limit. The cascade therefore remains non-default.
- **Held-out evidence is archived and visualized.** A tested archiver retained
  the complete matrix, gate receipt and every failure command/receipt/log at
  `paper/evidence/cascade_native_screen_heldout_v1`; an independent check found
  40/40 manifest hashes valid. `figures_v2` is the accepted six-panel SVG/PDF/PNG
  with 188 editable SVG text nodes and no raster image elements. `figures_v1`
  remains an inspected draft because its SVG text was converted to paths.
- **Quantify development ablation is complete.** The committed source snapshot
  ran once at `/Volumes/T7/Codex/TandemX/results/quantify_calibration_development_v1_20260907`:
  108/108 successful public-command executions, 5,940 family rows and 36 strata.
  Across 1,485 family conditions per method, mean absolute relative error was
  0.401898 for total-bases depth, 0.359640 for oracle error survival, and 0.365728
  for empirical controls. Controls improved all three genome means but regressed
  at nominal 1x. Controls plus oracle produced exactly the same 1,485 estimates
  as controls alone while increasing median runtime from 1.771 to 4.934 s.
  Diagnostic-spread truth inclusion was far below 0.95 and remains explicitly
  not a sampling CI. The hash-checked compact archive and accepted editable
  six-panel figure are in `paper/evidence/quantify_calibration_development_v1`.
- **FASTA speed replay has exact scientific parity.** After commit `6863155` and
  both hosted workflows passed, the same 108-command development matrix was
  replayed once into
  `/Volumes/T7/Codex/TandemX/results/quantify_calibration_development_v2_fast_fasta_20260907`.
  All 5,940 metric rows, 108 copy-number products, three control panels and every
  non-resource summary field are byte/field-identical to v1. Driver time changed
  from 610.507 to 349.217 s (-42.80%, 1.748x) and peak RSS from 108.234 to
  106.844 MiB (-1.28%). Median time changed 4.658 to 1.944 s (-58.26%) for oracle
  error and 4.934 to 2.178 s (-55.85%) for controls plus oracle. This is one
  same-machine replay, so it is an engineering check rather than publication
  timing. The 112 product hashes, resources and replay receipts are compacted at
  `paper/evidence/quantify_calibration_fast_fasta_replay_v1`.
- **Frozen quantify validation passed; seeds 6401--6403 are consumed.** Commit
  `9eda196` and both hosted workflows passed before the three 10-Mb validation
  genomes were generated. The input audit verified 30 manifests and all 93
  declared payload files (2,380,726,368 bytes). The once-only result is at
  `/Volumes/T7/Codex/TandemX/results/quantify_depth_gated_validation_v1_20260907`:
  54/54 successful public commands, 2,970 raw rows, 1,485 candidate rows and
  9/9 passed predeclared gates. Candidate MARE was 0.363222 versus 0.408767 for
  total-bases normalization; each seed improved by 0.037855--0.055248. Paired
  outcomes were 661 improved, 495 equal and 329 worse. Nine 1× conditions used
  the baseline branch and 18 5×/20× conditions used controls. Ungated controls
  had slightly lower aggregate MARE, 0.359443; retain this limitation and never
  tune on 6401--6403. Reserved 7401--7403 remain refused.
- **Validation evidence and public option are published.** The compact archive
  at `paper/evidence/quantify_depth_gated_validation_v1` has 67
  hash-checked entries, all 54 execution-artifact rows, compact dataset manifests,
  the input audit and independently recomputed decisions/resources. `figures_v3`
  is the accepted six-panel SVG/PDF/PNG with 94 editable SVG texts, no raster
  elements, and 54 uniquely paired runtime/RSS points. v1 is rejected for legend
  overlap; v2 fixed the layout but lacked unique panel-source resource keys. The
  opt-in `quantify --single-copy-min-depth 2` implementation records a total-
  bases fallback below the threshold without changing the controls default.
- **Guarded cascade speed development is complete.**
  Fresh development seed 1201 was used for three identical 96-run TandemX/
  TideHunter matrices. The baseline at
  `/Volumes/T7/Codex/TandemX/results/cascade_speed_profile_development_v1_20260907`
  had a paired runtime geometric-mean ratio of 2.357698. An unguarded 30%-span/
  95%-identity revision at `cascade_speed_profile_development_v2_20260907`
  reduced the ratio to 1.933332 but worsened 0.1%-indel boundary MAE from 1.364
  to 14.221 bp, so it is rejected and retained. The guarded candidate at
  `cascade_speed_profile_development_v3_20260907` requires at least 95% valid
  shifted columns and at most 2% unit-span residual in addition to the existing
  composition gate. It completed 96/96 runs, reached ratio 1.989148, retained
  array/family results, had minimum positive base-union F1 0.997597 and maximum
  positive boundary MAE 2.921 bp, and made no negative-control calls.
- **Cascade selection evidence is reusable.** The observable-feature audit at
  `fast_path_audit_v2` contains 1,600 read rows and accepted 403 development
  reads; truth-side scoring found zero negative, wrong-single-array or multi-
  array acceptances. Three cProfiles identify native banded alignment as the
  main baseline cost. The compact archive at
  `paper/evidence/cascade_gap_free_development_v1` contains the baseline,
  rejected intermediate, guarded candidate, audit, branch counts, profiles and
  decision; an independent check found 29/29 manifest entries hash-valid.
  `benchmarks/configs/cascade_gap_free_validation_v1.yaml` froze seed 2201 and
  14 gates.
- **Frozen cascade validation passed; seed 2201 is consumed.** Commit `054b935`
  was pushed to the working branch and `main`; hosted runs `34095006477` and
  `34095171697` passed Ubuntu/macOS Python, Rust and wheel jobs before the
  validation directory was created. The once-only run at
  `/Volumes/T7/Codex/TandemX/results/cascade_gap_free_validation_v1_20260907`
  completed 96/96 commands and passed 14/14 gates. TandemX had minimum positive
  array recall/precision 1.0, minimum base-union F1 0.997445, maximum positive
  boundary MAE 2.35 bp and zero negative-control calls. Its paired runtime
  geometric-mean ratio to TideHunter was 1.977877 and direct-child peak-RSS ratio
  0.396864. Some condition-level runtime ratios exceeded 4, so retain the
  distribution-level and synthetic-data limits. Reserved seeds 3201--3203 remain
  untouched. The compact archive has 13/13 valid manifest hashes at
  `paper/evidence/cascade_gap_free_validation_v1`; `figures_v2` is the accepted
  six-panel editable SVG/PDF/PNG, while v1 is retained after its generic NA label
  was clarified. The run's revision warning is caused by two untracked local
  native-extension files included in the source snapshot; both are hash-recorded.
  The scoped tracked source was clean at `054b935` before execution.
- **Real-data source audit.** Exact ENA metadata for PRJEB50694, PRJNA751841,
  PRJNA953663 and PRJNA919186 is archived at
  `paper/evidence/retrospective_collapse_source_audit`. Five manifest files and
  all 12 selected PacBio genomic-WGS run rows were independently verified. The
  metadata supports accession, platform and cultivar/material matching only; it
  does not establish identical DNA extraction, individual plant, stock or
  BioSample between historical and newer assemblies. All five manifest entries
  were independently rehashed. Donor-matched collapse truth remains unavailable.
- **Verification state and exact next actions.** Pre-cascade commit `aabaa29`
  passed all 522 Python tests, compileall, diff check, CLI help, 72 local
  Markdown links and both hosted Ubuntu/macOS workflows. At that checkpoint
  only `.codex/` remained untracked locally.
  Rust format and release Clippy with warnings denied pass. `cargo test --release`
  compiles but cannot launch locally because this conda build has only
  `libpython3.11.a`, not the requested dynamic dylib; hosted CI must execute the
  16 Rust tests. The release wheel builds and its isolated install/import plus
  complete toy workflow pass; wheel SHA-256 is
  `7c80951f09c3c52779dc5b81ea1492bd76a3323a11d4f6942a065af5ac45b350`.
  Final Figure 10 PDF was rendered with Poppler and visually passed. Commit
  `054b935` passes all 531 Python tests in 62.54 s, compileall, diff check, Rust
  formatting/release Clippy and both hosted workflows. The post-validation
  archiver/figure/manuscript source passes 49 focused and all 533 Python tests
  in 60.11 s, compileall, diff check, 69 local Markdown links, archive hashes
  and visual/PDF checks. It is published as `9b09bfd`; working-branch/main runs
  `34097112749`/`34097146588` passed Ubuntu/macOS Python, Rust and wheel jobs.
  Only `.codex/` is untracked. Never tune against consumed 2201, 3101--3103 or
  6401--6403 values. Continue the real matched-donor, plant-scale-comparator and
  release-portability gaps.

## Sequence-clustering and scoring checkpoint (current)

### Exact-output elastic alignment optimization (published)

- The native elastic alignment now packs four traceback directions per byte and
  processes all selected periods for one read through one native call, reusing a
  single uppercase sequence buffer. Scores, bands, tie order, traceback and hit
  order are unchanged; focused native/Python parity checks passed.
- Complete Mo17 discovery replays were run once and are finished. On 11,680,888
  bp, runtime changed 17.348 to 12.361 s (-28.75%) and peak RSS 81.703 to
  62.484 MiB (-23.52%), with all six historical core products byte-identical.
  On 111,505,681 bp, runtime changed 134.684 to 100.039 s (-25.72%) and peak RSS
  178.438 to 154.766 MiB (-13.27%), with all seven products byte-identical.
- `paper/evidence/discovery_packed_trace_batch_v1` contains hash-checked compact
  receipts. Each comparison is one historical baseline and one replay, so it is
  an engineering result rather than a publication timing distribution. Existing
  diagnostics still show TideHunter faster; no general superiority is claimed.
- The replay helper now accepts a legacy baseline without the newer auxiliary
  `family_audit_summary.json`, while still requiring and comparing all six core
  products. Missing core output remains an error. A compact archiver validates
  execution, input/source receipts, summary counts and product hashes.
- Commit `a73398d` is on the working branch and main. Full local Python validation
  passed: 472 tests in 47.68 s; Rust formatting and Clippy with warnings denied
  also passed. Local `cargo test --release` compiled but could not launch because
  this conda environment lacks dynamic `libpython3.11.dylib`. Hosted runs
  `34048998182`/`34049010825` then passed Ubuntu/macOS Python tests, executable
  Rust checks and wheel builds, closing that environment-specific test gap.
- Follow-up workspace reuse was tested and rejected rather than folded into the
  reported optimization. Commit `9e10b9c` passed hosted runs
  `34049804480`/`34049815693`, and all replay products remained byte-identical,
  but the 111.506-Mb run changed 100.039 to 96.913 s while peak RSS increased
  154.766 to 164.063 MiB. A narrower scratch-only revision changed the 11.681-Mb
  run 12.361 to 12.479 s while RSS decreased 62.484 to 60.922 MiB. Both are
  metric tradeoffs, so commit `e760471` restored the `a73398d` alignment core.
  Its work-branch/main hosted runs `34050266939`/`34050278620` passed
  Ubuntu/macOS Python, Rust and wheel jobs.
  Failed-attempt results remain on T7 under
  `Mo17_11Mb_alignment_workspace_replay_v2_20260907`,
  `Mo17_111Mb_alignment_workspace_replay_v2_20260907` and
  `Mo17_11Mb_alignment_scratch_replay_v3_20260907`.

### Published article draft and resolved figure export checkpoint

- b85a266 is pushed to main and working branch. Both hosted CI runs
  34033024599/34033018828 passed; earlier 1231743 runs
  34031712723/34031704058 and f16596b runs
  34029886116/34029873822 also passed. The committed checkpoint has 426 local
  Python tests; unchanged native source last passed 15 tests, Clippy `-D warnings`
  and formatting.
- The temporary automatic-review usage error is resolved. Subsequent Git reviews
  succeeded; a read-only source diff proved the figure update only changes the
  legend/whitespace, and retry through the original approval channel succeeded.
  `factorial_joint_multik/figures_v2` is exported, visually inspected and archived
  locally: 101 editable SVG texts, no raster nodes, source TSV byte-identical to
  version 1. Version 1 remains retained. Figure/doc updates are published in 4aafb00.
- The initial manuscript and explicit readiness audit are published in 4aafb00.
  This is a development draft, not a completed paper or software release.
- Col-0R, Ey15-2R and YSD56 full QC and compact archives passed. Col-0R/Ey15-2R
  add libraries/materials, while YSD56 adds wild soybean as an eighth reported
  species. The cohort table now contains ten libraries/eight species, 14,937,608
  reads and 262,731,255,175 bp. All ten seed6101 whole-file sampling ladders are
  complete. YSD56 contains 11.670/110.438/1,099.010/11,050.418-Mb nested samples
  (0.0116×/0.1095×/1.0897×/10.9570× nominal total-base coverage against the
  1,008,523,555-bp assembly denominator). Col-0R contains 11.708/115.753/1,156.486/
  10,643.552-Mb nested samples; Ey15-2R contains 11.698/113.336/1,119.054/
  10,627.736-Mb samples. Compact plans, receipts and distribution tables are
  archived locally; sample FASTQ and ID files remain on T7.
- Chinese Spring 12.698/126.731-Mb and Lo7 11.640-Mb three-tool comparisons are
  complete and locally archived. At 126.731 Mb, TandemX used 139.529 s/223.828 MiB,
  TRF 332.287/311.938 and TideHunter 48.230/413.500. Nipponbare 11.418-Mb reference
  QC is also archived: 626/626 primary-mapped reads, 611 primary MAPQ 20+, primary
  span fraction 0.999158. This GCA reference lacks organelles, so zero organelle
  counts are a reference-content limit, not evidence of no organellar reads.
- A new cross-tool evaluator records actual consensus provenance and exact cyclic
  edit recovery. At 90%, Mo17 CentC and Nipponbare Rice358 queries are recovered
  by all three tools. Morex HvT01 is recovered by TRF/TideHunter and TandemX's
  candidate stage, while the final TandemX representative is 0.8983. The
  0.89/0.90/0.95 sweep identifies a threshold/family-representation boundary,
  not donor-matched false-negative truth. Code, tests and evidence are local.
- The frozen-438d078 Morex 1.169928-Gb three-tool run completed and normalized
  all outputs. TandemX used 1909.525 s/913.891 MiB, TRF 2409.145/286.172 and
  TideHunter 566.196/580.375. TandemX was 20.74% faster than TRF but took
  3.37-fold TideHunter's wall time and used 3.19-fold TRF's peak RSS. The compact
  14-file archive matches T7 source hashes. This is an unfavourable concurrent
  diagnostic, not an isolated ranking. The original IPK MorexV3 pseudomolecule
  FASTA is now complete: the 4,296,032,540-byte file matches the published
  SHA-256, and streaming QC found eight records and 4,225,605,719 sequence bases.
  Compact provenance is archived; same-study context does not establish an
  identical read/reference donor or satellite copy truth.
- Two six-panel figure families passed source/hash and visual QA. Cross-cohort
  input-QC version 2 reconciles 10 libraries/8 species/14,937,608 reads/
  262,731,255,175 bp and has 110 editable SVG text nodes with no raster node.
  Its PDF was re-rendered and visually inspected. Real diagnostic version 2
  contains 10 nested inputs/5 materials, 78 editable texts and no raster node;
  version 1 is retained as a rejected legend-overlap layout with byte-identical
  panel source. The manuscript now includes Figures S3-S4 and Tables S6-S7.
- Tested compact archivers now cover complete FASTQ QC, whole-file sampling and
  successful three-tool real diagnostics. The full-discovery parity replayer can
  override only the thread budget and still requires all seven products to match.
- Soybean source curation distinguishes inaccessible ZH13 GSA records, a Wm82
  assembly project with no ENA read rows, and Jack project RNA-seq from usable
  genomic HiFi. Wild soybean YSD56 `SRR28726931` provides 44,193,089,411 HiFi
  bases from the same leaf BioSample as ONT support. Exact ENA size/MD5/SHA-256,
  gzip/FASTQ structure, expected counts and zero duplicate IDs passed. Seed6101
  nested sampling against the 1,008,523,555-bp assembly denominator completed;
  all four ID hashes, row counts and nested-set relationships were rechecked.
- f16596b publishes sequence-native exact clustering indexing and compact integer/
  32-bit-array family-audit postings. The Morex
  115-Mb v1 fixed-order interface diagnostic retained exact complete-output
  parity but ran with acquisition load. A post-acquisition v2 at committed source
  `a10c309` again produced identical 4,380-family payloads: sequence-native versus
  word-bridge clustering was 12.198 versus 13.611 s and peak child RSS was
  58.047 versus 61.203 MiB. The 10.38% time and 5.16% RSS reductions are one
  fixed-order diagnostic, not publication timing. Its generic revision warning
  identifies only the two ignored local native build artifacts and retains their
  exact hashes; tracked source is represented by the commit and source digest.
  The original YSD56
  transfer failed at 19,152,111,724 bytes with an SSL EOF, resumed from that
  partial, and then passed exact ENA byte/MD5 and complete-file QC checks; the
  failure/resume history remains in the receipt. The IPK server ignores
  MorexV3 Range requests. The successful rerun proved the retained 2.126-GB
  prefix against the new full response before appending the suffix, then required
  the official SHA-256 before canonical rename and FASTA QC.
- The predeclared conditional-abundance held-out seeds 5101–5103 were executed
  once with the unchanged configuration and f16596b source: 177/177 commands,
  81 copy-number, 45 localization and 405 comparison family rows. Across 243
  positive and 162 control comparisons, TP/FN/FP/TN=208/35/8/154 (sensitivity
  0.855967, false-positive rate 0.049383, precision 0.962963). All 36 positive
  exact-copy localization rows had base recall 1.0, but 20×/1%-error 50%-retention
  sensitivity was only 2/9 and 1× full assemblies produced one false call per
  error tier. The compact verified archive is `paper/evidence/abundance_heldout`;
  these seeds are consumed and cannot be used to tune the next model.
- A development-only multi-k/depth decision calibration reused the completed
  4101–4103 inputs without rerunning the original k=21 or localization commands.
  The transparent hybrid changed TP/FN/FP/TN from 195/48/7/155 to
  203/40/3/159, sensitivity 0.802469→0.835391, false-positive rate
  0.043210→0.018519 and precision 0.965347→0.985437. All three leave-one-genome-
  out folds selected the low-depth threshold 0.5; multi-k alone was unavailable
  for 35 rows. This is development evidence only. Config `abundance_v2.json`
  froze the rule and calibration hashes with untouched seeds 5201–5203.
- After 1231743 and both CI runs passed, the 5201–5203 matrix was run once. All
  177 baseline commands and 27 multi-k conditions completed. The frozen rule
  changed TP/FN/FP/TN from 198/45/14/148 to 208/35/12/150, sensitivity
  0.814815→0.855967, false-positive rate 0.086420→0.074074 and precision
  0.933962→0.945455. Each seed improved TP without increasing FP, but seed 5202
  retained 12 FP; the multi-k-only arm had 15 unavailable rows. At
  20×/1%-error/50% retention, recall improved 3/9→9/9; 1× complete-control FP
  decreased 5→3/27. At 1×/0 or 0.1% error/50% retention, recall worsened
  6/9→4/9 in each stratum. The validated compact archives preserve all adverse rows,
  frozen calibration files and baseline receipts. These seeds are consumed.
- The expanded source checkout passes all 462 Python tests. The added tests cover
  compact multi-k evidence integrity, exact pairing in the domain-shift figure,
  and clean Morex recovery while retaining a divergent partial. Rust source is
  unchanged. Four additional tests validate compact reference archiving and
  reject altered FASTA, metadata and contig totals. A further regression test
  covers the legacy helper-path receipt that initially blocked v2 archiving.
  New checks validate localization-only compact archives and exact pairing of
  failed development, selected development, fresh held-out and classifier rows.
- A backward-compatible challenge extension now represents independent
  founder-to-unit substitutions and one or more interrupted same-family arrays.
  Truth scoring aggregates non-overlapping intervals and copies by family. The
  predeclared configuration crosses 1/3/5% unit divergence with one/three array
  segments on unused held-out seeds 5401–5403, retaining the frozen 4101–4103
  model hashes and all original coverage/error/retention levels. Commit c15dad7
  passed hosted Ubuntu/macOS CI before seeds 5401–5403 were consumed once. The
  baseline completed all 1,062 commands and produced 486 copy-number, 270
  localization and 2,430 comparison rows. Frozen multi-k replay completed 162
  read conditions and 7,290 method rows. Sensitivity rose 0.843621→0.866255,
  but false-positive rate rose 0.643004→0.682099 and precision fell
  0.663073→0.655763. Full-assembly localization mean recall was 0.000793 at 3%
  divergence and zero at 5%. Compact paired evidence and an inspected six-panel
  figure retain the failure; seeds 5401–5403 are consumed and cannot tune the fix.
- The first `iid_base` localization-only development run consumed 5301–5303
  without repeating read simulation or quantification. All 90 commands completed.
  Positive-assembly precision was 0.999435 and all 54 absent-family rows had zero
  predicted bases, but full-assembly recall was 0.903608 and failed the declared
  0.95 gate. At 5% divergence, long arrays split into as many as 119 predictions.
  Development v2 retained the 0.90 proxy threshold and allowed exact anchors to
  bridge at most one monomer length; 500-bp planted interruptions remained outside
  that bound. Its 90/90 localization-only commands passed the development gates:
  full-assembly mean base recall 0.975728, positive-assembly mean precision
  0.999441 and 0/54 absent-family false positives. Maximum predicted fragments
  decreased from 119 to 4, while 5%-divergent three-segment arrays remained the
  weakest full-assembly stratum (mean recall 0.941017). Config
  `abundance_localizer_heldout_v1.json` freezes the development hashes, unchanged
  collapse model and then-untouched seeds 5501–5503. Commit 4e662db was pushed to
  both branches; hosted runs 34039414153/34039420693 passed Ubuntu and macOS
  before the seeds were consumed once. All 1,062 commands completed. Held-out
  full-assembly recall was 0.981033, positive-assembly precision 0.999565 and
  absent-family false-positive rate 0/54; the weakest 5%-divergent three-segment
  stratum was 0.949274. The frozen classifier changed TP/FN/FP/TN from
  875/583/12/960 to 1339/119/83/889: sensitivity 0.600137→0.918381, FPR
  0.012346→0.085391 and precision 0.986471→0.941632. This validates the bounded
  localizer under the tested model while retaining the classifier trade-off.
  Four compact archives and the inspected six-panel Figure 5 retain all rows;
  seeds 5501–5503 are consumed.
- The next classifier experiment is predeclared before data generation. Config
  `abundance_classifier_development_v1.json` reserves 5601–5603 for development
  and leaves 5701–5703 untouched. It limits selection to 20 transparent
  log-space single/multi-k blends and thresholds, requires no worse FPR or
  precision than single k=21, and adds full-development and leave-one-seed-out
  acceptance gates. The full 1,062-command baseline completed with all commands
  successful and the frozen config SHA-256 `3210fab6555a...c8c8e54`.
- The first 5601–5603 multi-k execution exposed a workflow failure after all 162
  estimates had been computed: the shared evaluator attempted the older hybrid
  calibration before writing paired inputs, and no legacy threshold met its FP
  constraint. The partial result is retained on T7. The recovery patch makes the
  predeclared blend-grid mode explicit and checkpoints raw paired rows before any
  optional calibration; candidate definitions and acceptance gates are unchanged.
- The recovered paired table contains 4,860 rows and 45 explicit multi-k
  unavailable observations. Classifier development v1 completed and failed its
  scientific gate: the pooled choice improved sensitivity 0.598765→0.714678,
  FPR 0.010288→0.008230 and precision 0.988675→0.992381, but selection differed
  among all three leave-one-seed-out folds. Cross-validated FPR rose to 0.021605
  and precision fell to 0.980374.
- Robust selection v2 is frozen as a post-v1 development refinement in
  `abundance_classifier_robust_selection_v2.json`. It requires FPR and precision
  preservation within every development seed and maximizes the worst seed-level
  sensitivity gain. It selected alpha 0.5/threshold 0.5 and passed all six
  development gates: worst-seed sensitivity delta 0.076132, maximum seed FPR
  delta 0 and minimum seed precision delta 0.001558; full deltas were +0.088477,
  -0.003086 and +0.004387, respectively.
- `abundance_classifier_heldout_v1.json` freezes the v2 classifier and selected
  localizer hashes. Both guards re-hashed and recomputed development gates before
  output creation; commit 300e48d and hosted runs 34043678448/34043693408 passed
  before seeds 5701–5703 were consumed once. All 1,062 baseline commands and 162
  multi-k conditions completed. The frozen blend changed TP/FN/FP/TN from
  916/542/29/943 to 1044/414/45/927: sensitivity 0.628258→0.716049, FPR
  0.029835→0.046296 and precision 0.969312→0.958678. Overall and worst-seed FPR/
  precision gates failed; seed 5703 had FPR delta +0.040123 and precision delta
  -0.025985. All 16 additional false positives occurred at nominal 1×. The same
  baseline had full-assembly localization recall 0.973858, positive-assembly
  precision 0.999631 and 0/54 absent-family false-positive rows. A compact
  evidence chain and inspected six-panel Figure 6 retain the result. Seeds
  5701–5703 are now consumed and cannot be retested as independent validation.
- Post-hoc stratification showed that the added false positives were confined to
  nominal 1×. Config `abundance_classifier_depth_gated_development_v3.json`
  therefore declares a transparent development hypothesis before any new seed
  is touched: use the unchanged single-k21/0.6 decision below estimated depth 2,
  and the alpha-0.5/threshold-0.5 blend otherwise. On the six now-consumed
  development genomes, this derived rule keeps all 39 baseline false positives,
  changes TP/FN from 1789/1127 to 1966/950, and raises precision
  0.978665→0.980549. This is post-failure development, not validation. At that
  point, seeds 5801–5803 remained untouched and were reserved for a future
  frozen test.
- `abundance_classifier_depth_gated_heldout_v1.json` now freezes that exact v3
  rule, all five development-artifact hashes, the unchanged localizer and the
  original held-out gates before generating any new data. Commit 62892a6 was
  pushed to both branches and hosted runs 34046154395/34046160309 passed Ubuntu
  and macOS before one-time use of 5801–5803. All 1,062 baseline commands and
  162 multi-k conditions completed. The frozen rule changed TP/FN/FP/TN from
  867/591/16/956 to 953/505/16/956: sensitivity 0.594650→0.653635, FPR remained
  0.016461 and precision 0.981880→0.983488. Every seed passed the predeclared
  sensitivity, FPR and precision guardrails; the weakest sensitivity gain was
  0.030864. The 810 estimated-depth-below-2 rows were unchanged, gains were
  0.072016 at nominal 5× and 0.104938 at 20×, and no standard-depth fallback
  occurred. Full-assembly localization recall was 0.980059, positive-assembly
  precision 0.999617 and absent-family false-positive rows 0/54. The compact
  hash-checked archive and inspected six-panel Figure 7 v2 are in
  `paper/evidence/abundance_classifier_depth_gated_validation_v1`. This passes
  the frozen known-catalogue IID simulation gate, not biological validation or
  an external-tool superiority claim. The complete Python suite passes 470
  tests after archive and figure QA. Evidence commit 3b21825 passed hosted
  Ubuntu/macOS Python, Rust and wheel runs 34047548628/34047562996 on the work
  branch and main. Seeds 5801–5803 are consumed.

### Completed interval calibration and reference-concordance QC (local, not pushed)

- Source d2aa4d6 was pushed to both branches. Hosted CI 34021933081/34021929197
  passed on main and the work branch. New reference-QC code
  and tests remain uncommitted. Full suite389 passed in86.59 s; Rust unchanged
  since the13-test/clippy/fmt checkpoint. Plot layout change after tests is
  documentation/visualization only; its second render did not execute.
- Frozen d2aa4d6 joint-read replay completed all27 conditions/1,485 family
  conditions in287.624 s/182.109 MiB. Points agree with the previous multi-k
  model. Intervals available/covering truth:1×9/9;5×66/60;20×424/401, each with
  495 total conditions. Overall499 available,920 insufficient support,66 no
  finite fit. Conditional20× coverage94.58%, availability85.66%; no unconditional
  95% calibration claim. Complete source tables and six-panel diagnostic SVG/
  PDF/PNG are archived locally in `factorial_joint_multik`.99 editable SVG texts,
  no raster nodes. Panel A legend overlap remains a layout issue in version1.
- Reference-mapping QC uses strict disk-backed PAF/CIGAR validation, preserves
  unmapped denominators and missing MAPQ, excludes deletions from target coverage,
  and reports primary organelle/other-reference query spans and their overlap.
  Col-0N controls completed at11.766 Mb and118.497 Mb. Larger sample maps7,554/
  7,557 reads;26,290,920 query bp (22.1871%) have organellar primary spans,
  with4,936 bp overlapping other-reference spans. This raises a concrete total-
  library normalization concern, not a proven22.2% copy-number bias/correction.
  Archives include exact development helper snapshots under `reference_mapping_diagnostics`.
- Wheat whole-file sampling completed through11.215 Gb; receipts archived in
  `ChineseSpring_input_qc/sampling`. Victoria11.361-Mb three-tool pilot completed
  (TX16.827 s/122.97 MiB, TRF44.997/165.78, TideHunter11.097/257.20).
  Victoria110.203-Mb run completed:TX164.738 s/250.45 MiB, TRF362.689/176.72,
  TideHunter81.982/621.45. Both Victoria sizes are archived locally. TX uses
  more RSS than TRF in this larger oat sample; do not imply a universal memory gain.
  Lo7 sampling, Morex reference transfer and Col-0R/Ey15-2R full QC remain active.
- Automatic approval review rejected the requested T7 figure-v2 render because
  the review service hit its usage limit (tool reported next reset Sep13,
  2026 11:23 AM; timezone not stated). It did not execute. Do not bypass that
  rejection through an alternative route. Unaffected read-only checks and local
  evidence/code writes continue. Git commit/push and new T7 writes requiring
  review await restored review availability; user authorization already exists.
  This is not scientific completion and does not justify abandoning the goal.

- An evidence-backed article draft now exists at `paper/manuscript.md`, with
  abstract/background/results/discussion/methods, existing multipanel figures,
  complete source-table links and a source-derived seven-species Table1. All
  seven input rows and cited local evidence paths were checked. This draft is
  explicitly incomplete; `paper/submission_readiness.md` retains every major
  scientific/release requirement. These files remain local and uncommitted.

### Joint-read uncertainty, complete TRASH evaluation and seven-species file QC

- Previous9fbd7a0 is pushed to main and the working branch; hosted
  CI34021111743/34021111471 passed. New joint-read multi-k research API preserves
  the point model but estimates sampling variance from full read-level cross-k
  moments. Sparse Rust counts match a naive word oracle through k31; Python/
  Rust collectors and batching agree with the existing point estimator. An
  identical-read regression prevents spurious tiny intervals from cancellation.
  Native13 tests, formatting and clippy-Dwarnings passed; only tandemx-dev rebuilt.
  Full Python suite after the joint and native-consensus changes:381 passed
  in82.51 s. No public quantification default changed.
  Actual27-condition coverage calibration is prepared but not yet run here.
- Committed96fd5dc Mo17 full1.129-Gb replay completed:1,241.876 s/697.109 MiB;
  all seven complete products are byte-identical to the prior indexed run.
  Time decreased23.61% from1,625.652 s, peak RSS increased0.81%. This includes
  native clustering and alignment changes and remains slower than prior
  TRF/TideHunter runs. Archive: `Mo17_alignment_workspace/full_1129Mb`.
- TRASH1 default10-Mb s6301 run completed in2,299.44 s/513.227 MiB. All11,257
  native units match the reference with audited−1bp offset. Unit-base precision
  0.99937075/recall0.99270646. Primary family recovery43/55; including native
  secondary consensus48/55. Period+region matching26/55 is distinct from base
  recall. Fractional native period171.5 is preserved, not rounded or omitted;
  the first failed integer-parser receipt and later evaluations remain archived
  under `TRASH_factorial_s6301`. No cross-input/platform speed ranking.
- Morex115.272-Mb three-tool diagnostic completed:TX159.600 s/173.03 MiB,
  TRF274.546/179.00, TideHunter60.409/489.81. Added compact source/native-hash
  evidence under `multispecies_real_diagnostics`; real accuracy remains unscored.
- Victoria complete file QC passed:398,850 reads/7,346,159,178 bp,
  median18,137/N5018,295, no N or duplicate archive IDs. This raises completed
  file QC to7 species, not7 biological accuracy validations. Whole-library
  seed6101 samples completed:11.361/110.203/1109.305 Mb. Raw/QC/sampling evidence
  archived in `Victoria_input_qc`. Col-0R and Ey15-2R complete downloads passed
  source MD5; full-file QC now runs. Lo7/wheat sampling and Morex reference
  transfer remain active; exact status must be checked on disk.

### TRASH2 truth evaluation and six-species file QC

- Source96fd5dc pushed to both branches; CI34020558074/34020558009 passed. Full suite
  after native-membership and cgroup-CPU additions:365 passed in81.70 s.
  Native core is unchanged since96fd5dc. Mo17 1.129-Gb exact replay is active.
- TRASH2 default10-Mb s6301 run completed:55/55 cyclic family recovery;
  all11,436 monomer sequences match reference at standard1-based coordinates.
  Native approximate windows match49/55 arrays; using explicit native unit-to-array
  IDs and outer unit bounds matches55/55, with3.1545-bp boundary MAE. Unit-base
  precision0.9998987 and recall0.9999437. Both evaluations are archived under
  `TRASH2_factorial_s6301`. Do not call coarse-window differences missed families.
- A repeated author control validated whole-cgroup CPU accounting (29.761 s
  versus3.85 s GNU time user+system). The repeated default output is not byte
  identical:7 region geometries match, units360 versus363, coverageJaccard0.996439.
  All outputs remain archived. R-worker stochasticity is a hypothesis, not an
  isolated cause; no claim that the new resource wrapper preserves exact native
  products. TRASH1 10-Mb run is still computing with bounded8-GiB allocation.
- Chinese Spring full QC passed:1,500,000 reads/24.954 Gb, median15,985/N5016,579,
  no N or duplicate archive IDs. This raises completed file QC to6 species,
  not6 biological accuracy validations. Wheat seed6101 whole-file sampling runs.
- Rice sampling completed through10.542 Gb; source receipts archived. New three-
  tool real comparisons completed on Morex12.542 Mb and Nipponbare11.418/113.627 Mb.
  Rice113.627 Mb:TX104.784 s/212.75 MiB, TRF140.957 s/333.86 MiB,
  TideHunter46.330 s/450.64 MiB. Morex115.272-Mb comparison still running.
  No real accuracy truth or isolated timing ranking is inferred. Compact evidence
  is under `multispecies_real_diagnostics`.
- Victoria oat full download passed source checks; complete7.346-Gb QC runs.
  Col-0R/Ey15-2R downloads, Lo7 sampling and Morex reference transfer remain active.

### Exact native alignment workspace and assembly-evaluation checkpoint

- Source99c0c66 pushed to both branches; CI34019743265/34019743319 passed.
  New in-place score/peak rows preserve every Python-reference traceback in
  30 focused tests. Full Python suite363 passed in94.10 s;11 Rust tests,
  clippy-Dwarnings and formatting passed. Only tandemx-dev was rebuilt.
- Actual111.506-Mb Mo17 profiles preserve all seven complete products. Original
  native alignment123.803 s/full151.262 s/174.75 MiB; reused rows104.134 s/
  139.015 s/173.36 MiB; final in-place93.374 s/124.860 s/184.22 MiB.
  Exact uncommitted-development source snapshots/hashes archived in
  `paper/evidence/Mo17_alignment_workspace`; both speed improvement and higher
  observed RSS are retained. Concurrent profiles are not final timing rankings.
- TRASH/TRASH2 are executing default de novo annotation on the same10-Mb
  factorial genome(s6301). The tested evaluator separates region matches,
  cyclic founder recovery and actual monomer-base union; TRASH1 window-grid/R
  extraction and native-peak/consensus-length alternatives are explicit. No
  benchmark completion or accuracy is yet claimed; inspect run receipts.
- Morex whole-file seed6101 sampling completed:572/5294/53715/493503 reads,
  12,542,086/115,272,341/1,169,928,427/10,744,915,847 bp. Receipts and distribution
  tables archived under `Morex_input_qc/sampling`. Rice/Lo7 sampling remains active.
- Nipponbare original GCA_034140825.1 reference acquisition and full FASTA QC
  completed; compact evidence archived under `Nipponbare_input_qc/reference`.
  The generic exact-version/official-MD5 downloader is tested. GCF differs by
  added MT/Pltd/B1; raw and assembly BioSamples differ in age/date; exact donor
  match is unresolved. Paper's chr9 rDNA model sequence is excluded from exact
  copy truth. See `docs/nipponbare_reference.md`. Other transfers/QCs still run.

### Five completed full-file QCs and tested comparator execution

- Source4ae8ca8 pushed to both branches; CI34016167100/34016167099 passed.
  Native index ablation completed on85,663 actual candidates,28,586 families,
  exact entire-payload parity. Python218.015 s/250.22 MiB versus native164.312 s/
  259.63 MiB:24.63% shorter clustering but3.76% higher child peak. Archived in
  `paper/evidence/Mo17_native_index`; no all-metric superiority claim.
- Full111.506-Mb live-pipeline replay with source4ae8ca8 also has all seven
  products byte-identical. cProfile shows native banded self-alignment123.8 s
  of151.1 s total; this is the next hotspot. Profiling adds overhead. Old-source
  indexed1.129-Gb comparison completed all3 tools, TX1625.65 s/691.50 MiB,
  TRF1073.00 s/321.12 MiB, TideHunter337.43 s/425.66 MiB. All seven TX outputs
  equal the previous1.129-Gb run. Both result sets retain unfavorable rankings.
- Full FASTQ QC now passed for5 species/materials: Mo17, Col-0N, Nipponbare
  (1,785,885 reads/32.966 Gb), Morex(896,701/19.525 Gb), Lo7(4,576,975/77.092 Gb).
  Their compact receipts/histograms are archived. These are file-integrity and
  distribution checks, not complete biological QC. Rice/Morex/Lo7 seed6101
  whole-library size ladders are running. Wheat full download passed and QC is
  active. Victoria oat, Col-0R and Ey15-2R full downloads have started.
- Original PanOat supplement confirms Victoria HiFi +Hi-C, exact raw sample
  SAMEA111508775 and assemblyPRJEB56706/GCA_947311595. Selected source cells and
  reviewer/QC lessons archived; Iso-Seq counts are excluded from genomic depth.
- TRASH1/2 offline author-example execution completed in the pinned ARM64
  container. The initial Docker startup timeouts are retained; their delayed
  empty containers were explicitly removed. Native monomer-sequence audits
  identify TRASH1's systematic-1bp offset(355/355) versus TRASH2 standard1-based
  agreement(360/360). These are installation/coordinate controls on an author
  human example, not plant benchmarks. No silent tool-output correction.
- New container runner, small-control coordinate audit and full-pipeline replay
  are tested. Full Python suite336 passed in75.65 s; native code unchanged since
  the11-test Rust/clippy/fmt checkpoint. Ongoing profiles/runs freeze source.

### Native candidate gate and known-source curation

- Source9d6ace8 CI34015330952/34015330671 passed. Its completed85,663-candidate
  storage replay has exact28,586-family/member output parity: dictionary172.420 s
  /361.50 MiB versus packed199.764 s/247.86 MiB. Both31.44% lower RSS and15.86%
  longer clustering time are archived in `paper/evidence/Mo17_compact_clustering`.
- New native exact index keeps the Python gate as reference; no thresholds or
  alignment rules change.11 Rust tests, clippy and20 focused Python tests passed.
  The isolated runner now offers an identical-alignment native-index ablation;
  actual85,663-candidate replay follows this source checkpoint. Full suite330
  passed in77.00 s; editable native rebuild stayed inside tandemx-dev.
- Four original ENA record pairs are archived under known_repeat_sources. The
  curator validated exact versions/MD5/bases/materials and extracted three queries;
  639-bp CentO cloneAF058902 is excluded. Cross-cultivar and native unit-boundary
  limits are explicit. Actual bank/receipt archived in known_repeat_sources/curated.
  No actual known-repeat matching result yet.
- Isolated ARM64 R4.4.3 TRASH/TRASH2 image built successfully; offline preflight
  and native-output validation follow. No actual comparator success claimed yet.
  Rice/Lo7/Morex full FASTQ QC, ChineseSpring download, Morex reference acquisition
  and indexed1.129-Gb Mo17 comparison remain active; inspect receipts.

### Compact clustering storage checkpoint

- Representative postings now store exact ID/multiplicity pairs in contiguous
  uint64 arrays and release the index before output construction. Alignment,
  candidate thresholds, ordering and membership rules are unchanged. Out-of-range
  fields fail, never wrap. Existing exhaustive/gated parity tests still pass.
- New isolated replay checks every other core/native source hash and runs both
  implementations in fresh children from frozen snapshots. It retains entire
  output hashes, command/resource receipts and serialized-score precision limits.
  Its toy source-isolation/parity/failure test passed. Full suite **322 passed in
  69.82 s**. Actual85,663-candidate memory/parity replay follows this source commit;
  no measured compact-index improvement is claimed yet.
- RiceSRR25241090 complete transfer passed source size/MD5 and SHA-256
  7bc90f777c995268d81f6f39d2baf9ed7c6fb123bcad0964930b7cfaeaa6b16d.
  Full32.966-Gb FASTQ QC is running. Rye/barley/wheat transfers and Morex reference
  acquisition remain active; inspect receipts before claiming completion.

### Disk-backed real-input evaluator

- Source8ff7ab6 pushed to both branches; CI34014817881 and34014817815 passed.
  Completed Mo17111.506-Mb evaluator replay preserves all three normalized TSVs
  byte-for-byte and all metrics (4.117 s /43.33 MiB). Col-0N5.280-Gb /336,613-read
  preparation-only stress check completed (161.180 s /54.59 MiB); this is not a
  discovery benchmark. Compact receipts: `paper/evidence/disk_evaluation_replay`.
- Old-source Mo171.129-Gb comparison is now complete: TandemX1930.701 s /
  710.59 MiB; TRF1019.346 s /323.52 MiB; TideHunter326.372 s /390.55 MiB.
  All unfavorable values retained in `Mo17_real_pilot/Mo17_1129Mb_related_audit`.
  Indexed-clustering1.129-Gb rerun started with8ff7ab6 and disk evaluator;
  not yet complete. Concurrent jobs mean all resource values are diagnostics.
- New Morex original-IPK reference downloader pins published SHA-256, archives
  source metadata and requires complete streaming FASTA QC. Its full4-GB-class
  transfer started; no reference-QC completion yet. Full suite **320 passed in
  68.06 s**; original/native core unchanged since8ff7ab6.

- Default real comparator evaluation now streams FASTQ/native rows into SQLite;
  read IDs and intervals are not retained in Python dictionaries/lists. The
  100,000-read cap remains only for the optional memory reference backend.
  Scope, native call order and duplicate/union semantics are unchanged.
- Tests compare all four adapters to an independent base mask, verify exact
  normalized TSV parity, accept 100,001 generated reads and reject duplicate IDs,
  changed/truncated input, late malformed rows and invalid coordinates. Failed
  outputs stay `.partial`; no false completion. **312 passed in 68.12 s**.
- New replay command retains exact input/source/native hashes and can check
  previously completed normalizations or only large-input preparation. Actual
  real-data parity and a >100,000-read stress replay follow this source commit.
  This removes an evaluator limit; it is not proof of whole-genome tool scaling.
- ce33eeaf9cce8d6a8b0fac5aa803af552a620e07 hosted CI 34014359486 and
  34014359473 both passed. Mo17 1.129-Gb old-source TandemX completed in
  1,930.701 s / 710.59 MiB, with85,663 candidates and28,586 families. TRF is
  still running; preserve all receipts before claiming the comparison complete.

### Exact multiset index for sequence-clustering candidates

- e57ae542f1d3ac50ada057a8d6155821f358fb5c pushed to both branches;
  hosted CI 34014078933 and 34014078840 passed. Serialized Mo17 replay completed:
  8,401 candidates, 4,219 families, all family/member values equal; stage13.132 s
  to5.294 s. Full independent 111.506-Mb rerun also completed all three tools:
  all seven deterministic TandemX data products byte-identical to the previous
  related-audit run. Evidence archived in `paper/evidence/Mo17_clustering_index`.
- Full rerun TandemX134.684 s/178.44 MiB versus previous142.919 s/155.45 MiB:
  observed memory increased despite shorter time. TRF88.068 s/258.03 MiB and
  TideHunter33.352 s/328.53 MiB. Concurrent jobs mean these are diagnostics;
  no final speed/memory ranking. The 1.129-Gb old-source run is still clustering.
- Wheat/oat run→experiment→sample XML chain archived in cohort_screen/xml_extended.
  SRR28200549 is explicit Chinese Spring HiFi; one 19.359-GB complete batch
  download started. Victoria oat is identified, but detailed protocol/reference
  matching is pending; CN25955 has a taxonomy-label discrepancy and is not
  silently treated as cultivated oat. Four full-file transfers (rice, rye, barley,
  wheat) remain active; only Mo17 and Col-0N have full raw QC completion.

- Canonical circular-word multiplicities are now accumulated through the existing
  representative index before oriented q-gram/alignment comparisons. The gate
  implies rejection by both old oriented tests; rounding/order/assignments and
  low-identity/short-query exhaustive paths are unchanged. Index multiplicities
  may change memory use; no memory improvement is claimed without measurement.
- Randomized/edge cases compare complete old/new family and membership objects
  at five identity thresholds; pruned pairs are checked against the original
  comparison. Full suite **296 passed in 66.88 s**. Native Rust is unchanged.
- Serialized replay uses rounded exported scores and is narrower than the now
  completed full 111.506-Mb live-pipeline parity check described above.

### Fixed multi-k replay and expanded discovery scoring

- Source 57ee82d8c78eb9ced1be135068c8c99567ca723e is pushed to both branches;
  hosted CI 34013737329 passed. All six factorial discovery executions/scorings
  completed for source6301, 5x clean/high-error conditions (~50 Mb each).
  All methods recover55/55 founders and eligible array recall1. High-error
  array precision: TandemX1, TRF.698867, TideHunter.920175; base precision remains
  .999952/.999988/.996805, respectively. Duplicate penalties are not wrong bases.
  TandemX50.706/58.898 s, TRF46.873/81.786 s, TideHunter14.434/17.137 s;
  concurrent-job timing diagnostics only. Compact evidence archived in
  `paper/evidence/factorial_discovery_s6301_5x`; no family-recall advantage here.
- Multi-k replay source tables and receipts are archived at
  `paper/evidence/factorial_multik_replay`. Six-panel PDF/SVG/PNG inspected;
  SVG has116 editable text elements and zero images after vector heatmap export.
  All input/output/script hashes verified; per-panel source and detailed legend
  retained. Initial raster heatmap draft is preserved under `/private/tmp` only.
- Mo17 1.129-Gb sourceb53a193 run scanned81,775 reads in1,269.837 s, producing
  85,663 candidate intervals; clustering is still running. Do not infer total
  runtime/family count/success before its completion receipt. Next scaling target
  is the repeated per-pair q-gram construction/selection in sequence clustering.

- 0aad0ab82c50f93f0c80351b3b4b3c1749c1908e pushed to both branches. Its fixed
  k15/21/27/31 replay completed all 27 inputs, 1,485 fits and 4,455 paired method
  rows. At 20x/2% unit divergence/high errors, mean bias changes -56.586% to
  +1.590%; mean absolute error 56.586% to 13.529%. These are IID-development
  results, not real-data calibration or external superiority.
- Trade-offs retained: 66 fits have zero support at some k (63 at 1x, three at
  5x), and two fitted positive slopes violate the simple loss model. At 20x,
  384/495 paired absolute errors decrease, 111 increase. No sampling CI is
  supplied and the public quantifier remains unchanged. T7 result directory:
  `results/factorial_multik_replay_v1_20260906`.
- New factorial discovery controller and independent truth views keep all planted
  bases, exclude partial-fragment reads only from the declared eligible-array
  endpoint, and use per-array consensuses with two sequence-recovery denominators.
  Exact bounded-edlib threshold decisions match the exhaustive score on tests.
  Source/tool/input hashes and failures are preserved. Full suite **294 passed
  in 65.11 s**; actual three-tool factorial pilot follows this commit.

### Completed factorial inputs and conditional scoring controller

- Source 6330d5625c102787c7a894822d1b3ba9599cc92d pushed to both branches;
  hosted CI 34013149414 and 34013149482 passed. All 27 factorial quantification
  commands now finished, producing 1,485 rows. Only 77 native k-mer spreads
  contain truth; at 20x/2% biological divergence/high read errors mean bias is
  -56.586%. All baseline rows/commands/configs are archived, not scientifically
  accepted. Three complete simulated input metadata bundles are also archived.
- Experimental `multik.py` and fixed k15/21/27/31 paired replay are implemented
  separately from public quantify. Independent naive counting, Rust/Python parity,
  multiplicity/ambiguity/invalid-input and small actual replay checks passed.
  Mean k21 and original median k21 are retained as ablations. No sampling interval
  or AI novelty is claimed. Actual large replay follows the tested source commit.
  Full suite after these changes: **291 passed in 65.20 s**.
- Col-0N completed pilot receipts are archived at `paper/evidence/Col0N_real_pilot`.
  Morex gap-analysis primary methods now establish the relevant HiFi/ONT/Illumina,
  optical-map, CENH3 and FISH evidence path. CentIER is added as a task-matched
  assembly centromere-region comparator; not yet installed/run.

- All three 10-Mb development genomes (6301/6302/6303) finished generation,
  each with 55 planted families and nine read conditions. Observed read bases:
  780,054,525 / 780,064,624 / 780,060,539 (total 2,340,179,688). Receipts under
  T7 `data/simulated/factorial_scale_sSEED_v1` are complete. These remain three
  controlled IID-background genomes, not three species or 165 plants.
- New `benchmarks.abundance.run_stream_quantify` uses only observed reads,
  founder catalogue, fixed genome size and k=21 in the public estimator; source
  occupancy/error metadata remain evaluation-only. It checks hashes, preserves
  failures and compares per-family bias and native k-mer spread to truth.
  Actual scoring is next. Full suite **288 passed in 64.69 s**.
- Col-0N 11.766-Mb related-audit pilot completed all tools: TandemX 17.520 s /
  89.47 MiB, TRF 12.981 s / 125.34 MiB, TideHunter 4.247 s / 254.58 MiB.
  Counts/coverage remain descriptive; these concurrent-job resource values are
  not final timing rankings. The Mo17 1.129-Gb pilot is still running.
- Complete-file downloads for Nipponbare, Lo7 and Morex continue; only Mo17 and
  Col-0N have completed full raw-file QC. Morex ERR4659246 is one technical CCS
  batch from PRJEB40587, not an independent five-plant cohort.

### Factorial scale generator and completed related-audit verification

- Full suite: **286 passed in 62.58 s**. New generator writes indexed FASTA
  with bounded memory and samples empirical lengths with independent error
  events/source coordinates. Config `factorial_scale_v1.json`: 10-Mb genome,
  54 length/copy/GC/divergence cells plus one 1.026-Mb array, nine coverage/error
  conditions. Development 6301/6302/6303; 7301/7302/7303 refused/reserved. Tiny
  independent base-mask, mutation-event, random-access and controller tests passed.
  Actual 10-Mb generation and scoring are next; no large-simulation result yet.
- Mo17 111.506-Mb exact-related audit completed with 9e0b87f source. TandemX
  142.919 s / 155.45 MiB; TRF 83.620 s / 210.84 MiB; TideHunter 32.653 s /
  269.34 MiB. All five main files and all 4,595 non-distinct rows match full mode
  byte-for-byte. 10,471 pairs scored out of 8,897,871 possible, audit about .69 s.
  Receipt and full-table streaming parity check are archived in the real pilot.
  These remain one-run diagnostics, not isolated performance rankings.
- Col-0N complete-library sampling finished: 752 / 7,557 / 74,600 / 336,613
  reads, 11,765,868 / 118,496,530 / 1,170,032,165 / 5,280,163,009 bp.
  All full-source hashes and exact IDs retained. `Col0N_input_qc/figures/input_qc.pdf`
  has four panels; PNG inspected, SVG 77 editable text elements / zero images,
  all input/output provenance hashes reverified. These are nested technical
  samples of pooled Col-0N, not independent plants or corrected nuclear depth.

### Complete-file QC, real scaling and exact-related audit (newest)

- 762618dc9728973da8219af27c09b658c521bf15 pushed to both branches;
  GitHub CI 34011654703 passed. Native-audit rerun on identical Mo17 840 reads
  completed in 17.348 s / 81.70 MiB versus baseline 258.829 s / 168.58 MiB;
  all six data files byte-identical. Repeated external calls remained faster.
- Next 8,084 reads / 111,505,681 bp completed all three tools: TandemX
  247.539 s / 165.47 MiB; TRF 87.493 s / 212.52 MiB; TideHunter 34.115 s /
  324.36 MiB. Scanning 127.534 s; 4,219 families produce 8,897,871 pair rows.
  Compact evidence is at `paper/evidence/Mo17_real_pilot`. All these single-run
  resource values are diagnostic during acquisition, not publication ranking.
- New opt-in `--family-audit related` on discover/run uses an exact inverted
  k-mer index and a monotone upper bound on existing rules. It emits all
  non-distinct pairs, preserving catalogue/warnings/collapse; full remains
  default. Audit receipt counts scored/pruned/omitted pairs and is hashed by
  the completion receipt. Dense catalogues can remain quadratic. Unit/CLI checks
  verify related rows against full results, both modes' main files, empty sketches,
  ambiguous sequences and pipeline propagation. **281 pytest passed in 65.14 s**.
  Native code unchanged since the previous 8-test Rust/clippy/fmt checkpoint.
  Real rerun still needed after this source commit.
- ERR6210723 complete QC passed: 933,904 reads / 14,646,601,458 bp, N50 15,663,
  median 15,514, GC .3676189867, N=0, duplicate archive IDs=0. Reported mean
  error .0015386005 is not empirical accuracy. Archive `paper/evidence/Col0N_input_qc`.
  Whole-library fractions .0008/.008/.08/.36 seed6101 are being generated;
  inspect T7 `data/subsets/ERR6210723_seed6101_v1/sampling_receipt.json` for completion.
  Last fraction targets approximately 40 nominal nuclear genome equivalents.
- Original Mo17 supplement downloaded and independently read. `curate_mo17_regions.py`
  extracts 64 satellite / 2 rDNA / 10 CENH3-centromere / 20 telomere regions,
  preserving raw coordinates and composition. All ten chromosome ends and
  genome size match the NCBI reference numerically. Half-open origin is inferred,
  not explicitly labelled; Table13's mixed size arithmetic is excluded. Two
  readers agreed on 2,033 selected cells. `paper/evidence/Mo17_published_regions`
  retains original XLSX and sources; these mixed regions are not base-level truth.
- Rice and Lo7 complete downloads continue under T7 raw data; inspect receipts.
  There are still only two completed raw-file QC libraries, not a validated
  eight-species cohort. No complete manuscript or production release exists.

### Real-input family-audit bottleneck (latest checkpoint)

- Commit 3d3201a0a168f5a36eea3f06b0a473576c42c962 is pushed to both branches;
  GitHub CI 34011141035 passed Linux/macOS. Mo17 real pilot completed all three
  tools on the same 840 reads / 11,680,888 bp. T7 result directory:
  `results/Mo17_complete_random_11Mb_pilot_v1_20260906`. Wall seconds / RSS MiB:
  TandemX 258.829 / 168.58; TRF 11.859 / 194.72; TideHunter 3.508 / 147.80.
  Calls/union coverage are descriptive, not accuracy. Concurrent QC occurred.
- Scanning took 14.514 s; 599 representatives then incurred 179,101 pair
  comparisons. New native ungapped audit, cached representative k-mer sets and
  streamed pair table preserve Python values/order and optional collapse.
  Quadratic output/time and candidate/warning state remain scaling limits.
- Checks: **273 pytest passed in 65.43 s**, 8 Rust tests passed, clippy passed.
  Exact parity tests cover randomized scores, ties, table bytes and warnings.
  Native extension installed in `tandemx-dev` with explicit `--no-user`.
  Rerun the identical real input after this source commit; no measured speedup
  from the new audit is claimed yet.
- ERR6210723 full FASTQ QC is ongoing. Rice SRR25241090 and rye ERR15194059
  complete transfers started, budgets 30/31 GB respectively; inspect receipts
  before claiming completion. Three Lo7 batches still represent one plant.

- New `distance.py` and `clustering.py`, with a Rust bounded global edit-distance
  kernel, implement explicit operational monomer clusters. Elastic `auto` now
  selects sequence clustering at 95%; legacy discovery/default is unchanged.
  Fixed observed representatives, no transitive merges, unique-read support,
  candidate FASTA and membership/ambiguity/support-filter audit are documented.
- `edlib==1.3.9.post1` installed inside `tandemx-dev`, pinned in the benchmark
  extra/environment. It is an independent benchmark evaluator, not detector code.
  Cyclic edit similarity exhausts query rotations/strands; independent DP and
  edlib comparisons test both the scorer and new cluster acceptance/rejection.
- Added period-independent union coverage and duplicated-base metrics. Old raw
  call and strict length endpoints remain visible. Re-scoring archive:
  `results/independent_rescore_v1_20260906` under T7. It predates the subsequent
  config-aware rescorer update; all these source configurations used IoU 0.5.
  ULTRA clean strict recovery 0 becomes cyclic recovery 1; tuned indel4 cyclic
  recovery remains 0. TRF indel4 raw array precision .5556 coexists with base
  union precision .999714. Do not interpret duplicate penalties as wrong bases.
- macOS CI for a5aaa50 failed two multi-array tests because they requested two
  threads above the runner's one-thread cap; Ubuntu passed. The tests now obey
  `discover_thread_limit()` while retaining repeated-run comparison. Hosted CI run 34007602407 for fd656b6 passed on both Linux and macOS,
  including source installation, pytest, Rust checks and wheel build.
- Initial local checks: 230 pytest passed in 59.50s, 7 Rust tests passed, clippy
  and fmt passed. Focused final clustering/CLI/rescoring checks: 16 passed.
  Final full suite: **231 passed in 60.12s**. Toy simulate/discover/validate
  succeeded with 29 candidate sequences, membership records and two families.
- `benchmarks/configs/sequence_clustering_v1.yaml` is ready for explicit 95%
  operational clustering, with three repetitions. Full rerun at `results/sequence_clustering_validation_v1_20260906`: 48
  successes, identical repeated outputs, all 13 positive scenarios at array
  recall/precision and cyclic monomer recall 1.0, zero calls on the three
  100-read negative datasets. Related-monomer recovery is corrected to 3/3. Held-out
  3101/3102/3103 were unused at that checkpoint and were later consumed by the
  frozen cascade audit described at the top. Candidate evidence supports a comparison
  at matched per-array output granularity instead of catalog versus raw calls.

## Broader comparator and figure checkpoint

- Source fd656b6a9171710a25ef48aa5b8b038a3d0537aa was pushed to both main and the
  working branch. Its CI is verified above. Current added scripts/evidence/docs
  form the next checkpoint; inspect Git for its exact final revision.
- Four-panel sequence audit finalized at
  `paper/evidence/sequence_validation/figures/sequence_audit.pdf` (also SVG).
  All values have source rows/hashes; final SVG has 95 text elements, zero image
  elements. PNG inspected; PDF opened in the app (queued). Resolution sensitivity
  gives two clusters at 90% and three at 91–100% for the same 70 exported candidates.
  Default 95% replay equals the actual catalogue; no biological family claim.
- SRF, KMC, minimap2 and k8 are available under T7 `tools/src`. Pinned commits,
  native hashes, full build logs and KMC compatibility patch are archived in
  `paper/evidence/comparator_builds`. KMC k=17/151 counts agree with an independent
  Counter on bounded N/RC/multiline examples. KMC's actual minimum -m is 2, not 1.
  See `docs/srf_workflow.md` for exact provenance and the mixed incremental build
  caveat; clean rebuild is needed for final publication benchmarks.
- `results/srf_workflow_pilot_v1_20260906`: four successful seven-stage runs,
  clean171/indel0.1%, ci100/ci20, seed1101, 100 reads each. Array and monomer recall
  1; base union recall about .951/.947. Workflow wall time 1.92–1.96s, maximum
  sequential native-stage RSS 48.1–51.4MiB, one repetition only. These are
  developmental observations, not paired superiority or biological CN evidence.
- `results/srf_workflow_development_v1_20260906`: all 16 scenarios x ci100/ci20,
  17 successful catalogues, four successful no-catalogue outcomes, 11 native SRF
  empty-count assertions. Original failures remain NA and are archived.
- `results/srf_workflow_empty_guard_v1_20260906`: full rerun with a disclosed
  workflow guard, 17 ok, four no_catalogue, 11 no_eligible_kmers, zero process
  failures. Empty dumps skip native SRF and create no fake FASTA. At ci20 all
  positives except divergent_units, indel4% and substitution5% recover all units;
  those three have zero eligible 151-mers at the count threshold. All negative
  controls have zero in-scope calls. Higher-order motifs outside 30–1000 bp are
  preserved; no HOR decomposition or equally budgeted tuning is claimed.
- New scripts: plot_sequence_audit.py, prepare_kmc_libcxx.py, verify_kmc.py and
  run_srf_pilot.py. SRF parser/empty-guard unit checks were added. Full pytest
  passed: **233 tests in 62.93s**. Core Rust files are unchanged since fd656b6.
  The final figure/SRF evidence checkpoint still needs its hosted CI after push.

## Native seed-histogram optimization

- Commit cf2ebd5ce5d485f31f4ceba3190bf1d68bfcb1a5 (SRF/figure evidence) was
  pushed to both branches; GitHub run 34008566906 completed successfully.
- A cProfile run on 100 low-indel reads is at T7
  `results/profiling_20260906/indel_01pct.prof`. It measured 2.069s including
  profiler overhead: Python extraction/histogram 1.018s cumulative, native
  self-alignment .703s. These are diagnostic profile values, not ranking timings.
- Elastic Rust mode now reuses the existing native extraction/histogram functions,
  exposed independently of the legacy refiner in `rust-core/src/spacing.rs`.
  Python is unchanged as the reference. Exact histogram/overflow parity covers
  376 randomized/edge parameter combinations plus explicit invalid-input checks.
  The native vector is bounded by read length plus bin-rounding slack even if
  the requested maximum period is huge. The GIL is released for seed processing.
- Current checks: 235 tests passed in 62.25s; 7 Rust tests, clippy and fmt passed.
  Added paired-performance gate test passed separately after that full run.
- Commit 3f2482f5c501e7d4f2e7e92d17db2bfa76e5cead is pushed to both branches;
  GitHub run 34008863635 passed Linux/macOS source tests, Rust checks and wheels.
- Paired source-snapshot experiment completed at T7
  `results/native_seed_paired_v1_20260906`: **96/96 executions, all six outputs
  byte-identical**, same seed-2101 inputs, shuffled old/new order, three repeats.
  Per-dataset median-wall speedups range 1.200–4.447, median 1.371. Positive
  scenarios range 1.200–1.657; the maximum is the low-complexity control.
  RSS changes range -13.09% to +2.03%, increasing in six of 16 scenarios.
  No concurrent local tests/profilers/builds were run during this experiment.
  This is same-tool optimization, not external or publication superiority.
- Compact raw/summary/provenance and a four-panel editable figure are archived
  in `paper/evidence/native_seed_paired`; all 16 scenarios are displayed. PNG
  visually inspected; SVG contains editable text and no raster images.

## Conditional abundance experiment (new)

- Independent generator/scorer/runner in `benchmarks/abundance`, with bounded
  genomic truth, uniformly sampled circular reads, sampling oracle, and five
  retained-copy assembly versions. Documentation: `docs/abundance_benchmark.md`.
- Development seeds 4101–4103; held-out 5101–5103 reserved. Known monomers are
  supplied to isolate quantification/localization, not discovery recovery.
- Current copy-number intervals are diagnostic-k-mer spreads, not calibrated
  sampling confidence intervals. The experiment measures bias, truth coverage,
  localization and false under-representation calls before modifying the model.
- Initial independent truth/scoring/CLI checks: four passed in 1.30s. Full
  suite now passed: **240 tests in 63.67s**. Core Rust code unchanged since
  3f2482f. Commit 8218beb2496c9d796b18833efc9a103b8e258286 was pushed to both
  branches; GitHub run 34009300960 succeeded on Linux/macOS.
- Full baseline at T7 `results/abundance_baseline_v1_20260906`: **177/177 command
  executions** (15 locate, 27 quantify, 135 compare), 81 CN/45 localization/405
  comparison family rows. Compact source tables at `paper/evidence/abundance_baseline`.
  Intervals contain truth in only **6/81** conditions; they are k-mer spread, not
  sampling CI. At 20x/1% substitution, mean signed CN error is -18.31%; estimator
  minus sampling oracle is -20.08%. Error-free mean estimator-minus-oracle is
  about -0.36% to -0.52% across coverages, but low-coverage sampling fluctuation
  itself is large. All positive assembly base recall equals 1 in this simple
  exact-copy setting. False under-representation calls and misses exist; see
  every conditional row. The baseline is **not scientifically accepted**.

## Real cohort, complete downloads and QC (ongoing)

- User further requires real and simulated evidence with enough species/data.
  `docs/cohort_and_qc.md` records a prospective target of 8-10 species / >=20
  independent materials, scale ladders and primary-source QC lessons. This is
  a project target, not claimed journal policy or achieved sample size.
- ENA metadata for Arabidopsis, rice, maize, barley, rye and wheat is at T7
  `data/manifests/cohort_screen_20260906`. Oat PRJEB56828 returned 285 PacBio WGS
  run records; ERR8666127/ERR8666125 metadata retrieved. Twelve queried tables
  are summarized with hashes in `paper/evidence/cohort_screen`. Never use
  run/BioSample counts as biological replicate counts; wheat
  has run-labelled BioSamples and Morex multiple cells under different IDs.
- Complete ENA downloads started for ERR6210723 (11.074 GB) and SRR15447419
  (5.346 GB). Inspect download_receipt.json and processes before claiming done.
  Mo17 filename includes `_subreads`, but original experiment title explicitly
  says CCS; initial observed Phred range is 3-93, not an accuracy validation.
  Full biological/library QC still required. Previous prefixes remain pilots.
- New complete-file downloader verifies size+MD5 before rename; full FASTQ QC
  checks gzip trailer, valid records, exact archive ID uniqueness on SQLite,
  length/N50, GC/N and reported quality distributions. Ten focused tests passed
  in .04s including truncation, duplicate IDs, invalid resume and checksum errors.
  These file checks do not validate species/material/coverage or molecule IDs.
- Full suite after adding these helpers: **250 passed in 64.58s**. Dataset
  downloads remain in progress; no full raw-read QC or final timing claimed.
- Col-CEN v1.2 plus known issues and README downloaded at T7
  `data/references/Col-CEN_v1.2` from pinned source abb9b614d91c8a0bbd05a199c694fb5eafb5fe30.
  Git blob hashes, full gzip/FASTA verified; 132,081,078 bp including ChrM/ChrC.
  SHA256 of compressed FASTA b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8.
  Nuclear denominator/known-issue masks must be explicit in comparisons.
- Selected ENA run/sample/experiment XML now retrieved after a transient TLS
  retry and archived under `paper/evidence/cohort_screen/xml`. Mo17 cultivar
  and CCS experiment title are explicit; rice is Nipponbare/AGIS-1.0; Lo7 is
  a single diploid plant with HiFi library metadata. More general independence
  of all libraries/specimens remains to be resolved from study supplements.
- Col-CEN nuclear Chr1-5 sum to 131,559,676 bp, separate from ChrM/ChrC.
- Added reviewer-relevant unitFinder (Genome Biology soybean study) to comparator
  matrix. SRF paper's k=101 CentC rescue is now explicit; k/count sweep remains
  required. Do not score assemblies/known-motif tools on unsupported de novo reads.

## Source, Git and storage

### Completed Mo17 input and read-cluster replay

- Commit 6f29db5bc3a77db834239d38e91925f80e6b6018 pushed to both branches;
  GitHub run 34010445317 passed Linux/macOS. Full suite after adding reference
  QC and real-input adapter checks: **270 passed in 62.62s**.
- SRR15447419 full QC passed: 407,670 reads / 5,624,644,958 bases; N50 13,968 bp,
  median 13,455 bp, GC .4592192482, N=0, exact duplicate archive IDs=0. Reported
  mean error probability .0028814411 is not empirical accuracy.
- Seed6101 whole-file hash sampling finished at T7
  `data/subsets/SRR15447419_seed6101_v1`: .002 = 840 reads / 11,680,888 bp;
  .02 = 8,084 / 111,505,681 bp; .2 = 81,775 / 1,128,793,699 bp. All source/ID/
  FASTQ hashes recorded. These are nested technical subsets, not extra plants.
- NCBI Mo17 T2T GCA_022117705.1 downloaded/full FASTA verified at
  `data/references/Mo17_T2T_GCA_022117705.1`: 645,529,558 compressed bytes,
  MD5 98b17477eb92144393bc3682e0c5e4dc, 10 contigs/2,178,604,320 bp, no ambiguity.
  NCBI report maps CM039150.1..CM039159.1 to chr1..chr10. Same-study/cultivar
  context is established; identical donor/extraction and identity to the alternate
  MaizeGDB release are not. CyVerse verification page / MaizeGDB403 prompted use
  of the separately public paper-linked NCBI record.
- Compact receipts/histograms/reference reports in `paper/evidence/Mo17_input_qc`;
  final four-panel `figures_checked/input_qc.pdf`/SVG visually inspected. SVG
  70 text elements, zero images. The initial cramped-caption figure was removed.
- Read-cluster replay completed on all 27 development inputs / 81 family rows:
  33 intervals available, 24 contain truth, 48 missing for inadequate support.
  At 20x/1% substitutions mean CN bias remains -17.97%, estimator-minus-oracle
  -19.74%. Archive `paper/evidence/read_cluster_replay`; no calibration or
  production quantifier acceptance. Held-out abundance seeds unused.
- New `run_real_comparators.py` verifies identical selected FASTA for elastic
  TandemX/TRF/TideHunter, retains raw/normalized intervals and explicit failures.
  It is capped at 100,000 reads and records descriptive observations only;
  inspect T7 runs before claiming an actual comparison completed.
- ERR6210723 full download now completed; full QC is the next step. Check the
  current receipt/process rather than older download-in-progress notes below.

### Read-cluster and sampling implementation checkpoint

- Commit 10bedca453a3e24f5065bf1f6b30c0d6578b9aa2 was pushed to both branches;
  GitHub run 34009894957 completed successfully on Linux/macOS.
- New shared bounded FASTQ parser computes raw SHA256 while parsing, including
  multiple gzip members and trailers. QC no longer needs a second hash scan.
- New `sample_complete_fastq.py` requires full QC, checks exact source hash/totals
  and creates nested seeded read-ID hash samples, deterministic gzip, ID audit,
  length/GC/quality histograms and explicit empty-sample status. Order-independent
  membership, reproducible bytes, missing/changed QC, truncation and denominator
  checks passed. Sampling is within included files, not a whole-study claim.
- `tandemx/quantify/read_moments.py` is a non-default Python research model:
  multiplicity-weighted mean target counts / read k-mer opportunities, with
  independent-read ratio SE. It keeps only family moments. Sparse/zero evidence
  has no claimed calibrated interval; error/library/genome-size bias is unresolved.
- Six read-model tests pass, including independent naive counts, residual formula,
  Bernoulli-model coverage, independent-genome replay and held-out/hash rejection.
  This is not empirical plant calibration. `evaluate_read_moments.py` snapshots
  source for a development replay of all prior abundance inputs; inspect results
  before claiming any improvement. Held-out abundance seeds remain unused.
- XML source files' executable bits were removed; their content is unchanged.
- Full Python suite passed: **260 tests in 66.02s**. Native code unchanged.
- SRR15447419 complete download finished: 5,346,359,125 bytes, source MD5
  7d25a15f3faa9ce4bd2bed9789e394d8, SHA256
  3047df8fc93cd0c87d03b968ad702a95a2edf6ae3474d92b08350827d847f05c.
  Full QC started against 407,670 records / 5,624,644,958 bases; inspect
  `data/raw/SRR15447419_complete/qc_v1/qc.json` before claiming QC passed.

- Source: `/Users/ycy/Codex/Sofw/TandemX`; GitHub `https://github.com/CongyangY/TandemX`.
- Active branch: `codex/publish-current-progress`. Previous checkpoint `f4bb14f`
  was verified on both GitHub main and the working branch. This document ships
  with the next tested checkpoint; inspect `git log` for its exact revision.
  SSH remote access works with sandbox escalation.
- `.codex/` is user-owned: never stage/delete it. Ignore caches/generated runs.
- New data/results: `/Volumes/T7/Codex/TandemX`, mounted with ~1.6 TiB free on
  entry. Source and `tandemx-dev` remain on internal disk.
- 553 legacy output files (120,621,151 bytes) copied and SHA-256 verified at
  `results/legacy_20260715`. Manifest: `provenance/legacy_copy_manifest.json`.
  Originals preserved. Existing comparator binaries remain in `.conda-benchmark`.

## Implemented this session

1. Modular independent challenge simulator/adapters/evaluator/runner in
   `benchmarks/challenge/`: disjoint development/validation/held-out seeds,
   16 scenarios, one-to-one array matching, separate read/array/sequence-family
   metrics, strict errors/timeouts and source/input/executable/resource records.
2. Bounded ENA prefix downloader with validation, metadata and subset checksum.
3. Valid negative discovery now succeeds with `discovery_summary.json`, empty
   catalog and header-only tables. Validation requires the receipt's hashes.
   Empty/malformed input still fails. Pipeline dependencies explicitly skip;
   forced/resumed negative reruns remove stale expected positive outputs.
4. LICENSE, CITATION.cff, declared existing visualization/benchmark dependencies
   and Linux/macOS conda CI. Check hosted CI after push; configuration alone is
   not evidence of a hosted test pass.
5. Four-panel development diagnostic with editable SVG/PDF and source tables,
   selectively versioned in `paper/evidence`. It is not a finished paper.

## Elastic checkpoint (2026-09-06)

- User explicitly requires broader reviewer-relevant comparators and improvement
  across speed, memory, accuracy and other metrics. `docs/comparator_matrix.md`
  records task matching, all metrics and the no-cherry-picking contract.
- New opt-in `--discovery-method elastic` in both `discover` and `run`; legacy
  remains default for now. Separate Python reference and Rust banded local
  alignment/global consensus kernels; multiple arrays per read; composition
  filter; support-based unit template; explicit uncalibrated confidence.
- New modules: `tandemx/discover/{alignment,consensus,elastic}.py`,
  `rust-core/src/elastic.rs`. Array identity includes gap columns. No AI claim.
- At data root, `results/elastic_development_v1_20260906` retains the first
  attempt: 3 AT-rich false-positive reads and one divergent consensus-length
  error. `elastic_development_v2_20260906` corrects both: all 13 positive
  scenarios at array recall/precision and family recovery 1.0; three 100-read
  negative datasets have no calls. Both runs have one repetition/seed 1101.
- `results/elastic_validation_v1_20260906`: 48/48 successful runs, seed 2101,
  three identical-output repetitions per scenario. Array recall/precision 1.0
  throughout; no negative calls. **Related-family recovery is 2/3**: clustering
  merged two distinct related monomers. Other family endpoints are 1.0. This
  is the immediate algorithmic deficit; do not claim complete validation.
- The family failure is at `runs/related_families_s2101/tandemx/rep1` in that
  validation directory. `families.tsv` contains only two 171-bp families with
  38 and 32 supporting reads. Existing `candidate_sequences_compatible` uses
  a shared-sketch threshold; an alignment-based family criterion is needed.
- ULTRA v1.2.2 built locally at `tools/build/ULTRA/ultra` under T7. Upstream
  source `tools/src/ULTRA`, commit `99418b9eb396aaaf59e4b793a481b4c7aa8a8104`.
  Provenance: `provenance/ultra_build_20260906.json`. No system package install.
  Default model and automatic indel-aware tuning pilots both ran on 10 reads,
  clean_171 and indel_4pct. Native streaming buffers set to 5000 bases and one
  queued window, based on observed read length. Default resource estimate at
  period1000 was 4.09GB; measured bounded pilot peaks ~421 MiB.
- `results/ultra_pilot_v1_20260906`: clean array recall 1.0, indel4% 3/7.
  `results/ultra_tuned_pilot_v1_20260906`: both array recall/precision 1.0;
  166.58/176.81 sec including 18 grid settings plus shuffles. Strict equal-length
  family recovery is zero even on clean reads; consensuses differ by ~1 bp.
  **Add an independent gapped/circular homology endpoint before interpreting
  this as missing families.** Do not compare 10-read resource rows to 100-read
  runs. Future comparisons need default, tuning cost and frozen tuned inference.
- Four-panel baseline/elastic development figure and source tables archived in
  `paper/evidence/elastic_development_corrected`; no raster SVG elements,
  225 editable text elements. PNG inspected. Validation failures are archived
  separately and explicitly described in `paper/evidence/README.md`.
- Early elastic runs have source hashes but no complete dirty source snapshot;
  they are diagnostic only. Updated runner snapshots Python/Rust/native/build
  files, verifies copied hashes, and runs TandemX with snapshot PYTHONPATH and
  working directory. Final publication runs require a committed source revision.

## Checks

- Latest elastic/ULTRA/snapshot checkpoint: **218 pytest passed in 60.63 s**;
  6 Rust tests passed, clippy `-D warnings` and fmt check passed. Earlier
  217-test runs preceded the final snapshot import-isolation test.
- The documented elastic toy tutorial executed simulate/discover/validate with
  exit 0 and two emitted families. Final snapshot/cwd benchmark smoke on T7:
  `results/elastic_snapshot_cwd_smoke_20260906`, three successful, deterministic
  validation indel4% runs; read-local array recall/precision and family recovery 1.
  The emitted run config confirms cwd is the archived source snapshot.
- The prior development wheel below predates elastic; no newer release-wheel
  validation or production release is claimed by this checkpoint.

- Previous checkpoint `f4bb14f`: GitHub CI run 34005639588 verified successful
  on both Linux and macOS (source install, pytest, Rust checks and wheel build).
  This does not cover the newer elastic checkpoint until its own CI completes.
- Entry: 170 pytest passed in 60.09 s.
- Expanded suite: 200 passed in 59.01 s; after stale-output cleanup the complete
  suite passed again (200 tests in 58.55 s). Focused pipeline regression: 9 passed.
- Rust: 4 tests passed; clippy `-D warnings` and fmt check passed.
  macOS cargo tests need `DYLD_FALLBACK_LIBRARY_PATH` pointing to
  `/opt/homebrew/Caskroom/mambaforge/base/envs/tandemx-dev/lib`.
- Editable install `pip install --no-build-isolation -e '.[test,benchmark]'`
  succeeded in `tandemx-dev`.
- Maturin requires explicit `--interpreter` (use env `sys.executable`), otherwise
  it selected system Python 3.13 despite conda activation. Do not distribute an
  unvalidated auto-selected wheel. The explicit Python 3.11 macOS ARM64 wheel
  passed temporary-target installation, Python/Rust import-origin checks and all
  seven toy workflow steps. Receipt at data root:
  `releases/development_20260906/verification_cp311/wheel_validation.json`.
  Wheel SHA-256: `e872da439d9b8a8e11c24ab47fdacb5c14340bb5fc6dd85926a4652a60407460`.
- Diagnostic SVG verified: zero raster image elements, 182 editable text nodes.
  Full PNG visually checked; panels linked to source-table SHA-256 hashes.

## Actual development benchmark

Full run: `/Volumes/T7/Codex/TandemX/results/challenge_v1_development_20260906`.
144 attempts (16 scenarios x 3 tools x 3 repetitions), 135 successes and 9
pre-fix TandemX errors on negative controls. Runner correctly returned nonzero;
all rows, logs and failures were retained. Do not relabel old failures after a fix.

| Scenario, seed 1101 | TandemX array recall | TRF | TideHunter |
| --- | ---: | ---: | ---: |
| Clean/substitution-only | 1.0 | 1.0 | 1.0 |
| 0.1% total indels | 0.7000 | 1.0 | 1.0 |
| 1% total indels | 0.02857 | 1.0 | 1.0 |
| 4% total indels | 0 | 1.0 | 1.0 |
| Two arrays per read | 0.5 | 1.0 | 1.0 |

These are interval-aware **array** metrics (IoU >=0.5 plus period tolerance).
TandemX still recovered all three families at 1% indels and one of three at 4%.
Inspect per-read matches before changing period inference: much of the failure
is incomplete boundaries. The one-candidate-per-read model explains 50% recall
for two arrays. Raw overlapping/harmonic calls count against comparator precision.
Timings are developmental because source checks and a separate real-read pilot
overlapped; publication timing must be isolated with frozen code/settings.

Fix run: `results/negative_control_fix_20260906` at the data root. All 9 TandemX
reruns succeeded with zero calls on the same three negative datasets. Recall
with no positive/family denominator remains NA. No held-out seed has been used.

## Real data obtained

`data/raw/ERR6210723_prefix1000` at the data root: 1,000 Arabidopsis reads,
15,666,956 bases; FASTA size 15,710,688 bytes. SHA-256:
`4adca195ae3587c52bf503751fd82fe68ae755fe15c7d10f3c09d2d1a1a2305f`.
Source project PRJEB46164, run ERR6210723, sample SAMEA8961650. The study
identifies HiFi reads; ENA reports Sequel. Full remote FASTQ is 11,074,308,726
bytes; its MD5 and gzip trailer are **not verified** for this prefix extraction.
Prefix sampling does not establish unbiased library abundance.

Baseline: `results/ERR6210723_prefix1000_baseline/discover`: 320 candidates,
213 families, min-support=1, periods 30-1000, minimum span 100, one Rust thread.
Top family: 177 bp, 23 supporting reads. No known-repeat, abundance or collapse
validation has been established. Raw black-rye reads remain unavailable; do not
reuse assembly-alignment BAMs from other projects as raw reads.

## Next critical work

0. The current checkpoint and both hosted CI runs pass. Sequence clustering and
   independent gapped/union endpoints are implemented and tested; do not redo the
   earlier completed fix. Profile the cascade speed gap only on development/new
   seeds, then expand independent distributions, matched per-array consensus
   scoring and task-matched abundance/copy-number evaluation. SRF high-k misses
   on highly mutated inputs require broader sensitivity study, not a claim of
   general inferiority. Seeds 3101--3103 and 6401--6403 are consumed.

1. Keep meaningful source/tests/documentation checkpoints on GitHub; source
   snapshots are now available for subsequent benchmark runs.
2. Elastic now addresses indel boundaries/multiple arrays in development and
   validation, with parity tests; investigate remaining family clustering and
   runtime costs before making it the default.
3. Use new development/validation seeds, freeze every decision before untouched
   evaluation and add independent families/processes for publication inference.
4. Add coverage/error/copy-number and engineered assembly-collapse experiments,
   calibrated uncertainty and independently tested probe specificity.
5. Expand real plant validation with matched assemblies and curated repeats;
   measure resources on bounded inputs before larger datasets.
6. SRF family/abundance development workflows now ran as described above. TRASH
   assembly and applicable TAREAN short-read workflow comparisons remain unrun. Docker CLI is present but daemon not running.
   TRASH needs R dependencies; TRASH 2 also needs mafft/nhmmer. ULTRA, SRF, KMC, k8 and minimap2 comparator dependencies are now available. Keep incompatible task metrics separate.
7. Prior art: SRF already supports accurate-read satellite discovery/abundance.
   Optional AI requires transparent baselines, held-out evaluation, ablation,
   calibration and domain-shift evidence; do not add an AI label for novelty.
8. Complete release/reuse and manuscript gates as evidence permits. No complete
   manuscript, final comparative study or production-scale release exists yet.

The objective remains active; a development wheel or toy run does not meet the
user's mature-software and mature-paper endpoint.
