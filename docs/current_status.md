# TandemX current status and handoff

Updated 2026-09-06. Read completely after `AGENTS.md`, then verify Git/tests.
The user has authorized autonomous development and GitHub updates toward mature
software and a full evidence-backed paper (multi-panel figures and supplement).
The goal is **not complete**. Acceptance gates: `docs/release_program.md`.

## Source, Git and storage

- Source: `/Users/ycy/Codex/Sofw/TandemX`; GitHub `https://github.com/CongyangY/TandemX`.
- Entry branch: `codex/publish-current-progress`, HEAD `08e100d`. Entry remote
  main verified as `b9ef30a`. SSH remote access works with sandbox escalation.
  Recheck current revisions rather than assuming entry values are current.
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

## Checks

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

1. Commit/push tested code and evidence; inspect hosted CI. Do not stage `.codex/`.
2. Fix indel-sensitive boundaries and multiple arrays per read with an ablatable
   baseline and Python/Rust parity. Candidate consensus can survive when full
   array boundaries fail; panel d and `matches.tsv` provide concrete examples.
3. Development/validation seeds first, then freeze settings and use held-out
   3101/3102/3103. Add more independent families/seeds for publication inference.
4. Add coverage/error/copy-number and engineered assembly-collapse experiments,
   calibrated uncertainty and independently tested probe specificity.
5. Expand real plant validation with matched assemblies and curated repeats;
   measure resources on bounded inputs before larger datasets.
6. SRF family/abundance, TRASH assembly, and applicable TAREAN short-read workflow
   comparisons remain unrun. Docker CLI is present but daemon not running.
   TRASH needs R dependencies; TRASH 2 also needs mafft/nhmmer. No new comparator
   installation has been completed. Keep incompatible task metrics separate.
7. Prior art: SRF already supports accurate-read satellite discovery/abundance.
   Optional AI requires transparent baselines, held-out evaluation, ablation,
   calibration and domain-shift evidence; do not add an AI label for novelty.
8. Complete release/reuse and manuscript gates as evidence permits. No complete
   manuscript, final comparative study or production-scale release exists yet.

The objective remains active; a development wheel or toy run does not meet the
user's mature-software and mature-paper endpoint.
