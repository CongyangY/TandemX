# TandemX current status and handoff

Updated 2026-09-06. Read completely after `AGENTS.md`, then verify Git/tests.
The user has authorized autonomous development and GitHub updates toward mature
software and a full evidence-backed paper (multi-panel figures and supplement).
The goal is **not complete**. Acceptance gates: `docs/release_program.md`.

## Sequence-clustering and scoring checkpoint (current)

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
  3101/3102/3103 remain unused. Candidate evidence supports a future comparison
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
- New `compare_discovery_snapshots.py` will compare frozen old/new sources on the
  same datasets with shuffled order, repeated timings, RSS, user/system CPU and
  byte equality of all six discovery outputs. It is not yet run at this snapshot;
  do not claim a speedup until measured. It must run without concurrent local
  tests, profiler or compilation. Final publication superiority is still unproven.

## Source, Git and storage

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

0. Verify this checkpoint's hosted CI after pushing. Sequence clustering and
   independent gapped/union endpoints are now implemented and tested; do not redo
   the earlier completed fix. Next expand independent distributions, matched
   per-array consensus scoring, performance profiling and task-matched abundance/
   copy-number evaluation. SRF high-k misses on highly mutated inputs require a
   broader parameter sensitivity study, not a claim of general inferiority.
   Held-out 3101/3102/3103 remain unused.

1. Keep meaningful source/tests/documentation checkpoints on GitHub; source
   snapshots are now available for subsequent benchmark runs.
2. Elastic now addresses indel boundaries/multiple arrays in development and
   validation, with parity tests; investigate remaining family clustering and
   runtime costs before making it the default.
3. Development/validation seeds first, then freeze settings and use held-out
   3101/3102/3103. Add more independent families/seeds for publication inference.
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
