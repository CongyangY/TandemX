# TandemX Current Status and Codex Handoff

Last updated: 2026-07-15

This file is the durable handoff for a new Codex conversation. Read it after
`AGENTS.md`, then verify the live repository state before making changes.

## Repository and Git state

- GitHub repository: `https://github.com/CongyangY/TandemX`
- Default branch: `main`
- Current publication branch: `codex/publish-current-progress`
- Baseline on GitHub before this handoff: `b9ef30a`
- Local implementation commits included in the publication branch:
  - `ac4b7f8 Optimize streaming algorithms and strengthen validation`
  - `1941c6e Make external benchmarks publication-ready`
- `.codex/` is a user-owned, untracked local directory. Do not stage or delete it.
- Benchmark outputs under `benchmarks/results/` are reproducible local artifacts
  and are intentionally ignored by Git unless a future task explicitly selects
  a small publication artifact for version control.

At the start of a new session, run:

```bash
git status -sb
git log --oneline --decorate -8
git diff --check
```

Do not assume that the branch or commit identifiers above are still current
without checking them.

## Implemented progress

The current MVP is an end-to-end, read-first and assembly-aware tandem-repeat
workflow. It provides simulation, de novo discovery, diagnostic k-mer copy-number
estimation, assembly localization, read-versus-assembly comparison, FISH probe
ranking, visualization, validation, and the combined `tandemx run` workflow.

The two latest implementation commits add the following material changes:

1. Streaming and bounded-memory sequence processing, incremental output, and
   reduced retention of read-level state.
2. Faster candidate-period generation and refinement in discovery.
3. Python/Rust behavior parity checks and stable deterministic ordering.
4. Stronger input/output validation, empty-input handling, reverse-complement
   coverage, reproducibility coverage, and clearer warnings.
5. A fair external comparison harness for TandemX, TRF, and TideHunter with an
   explicit one-thread policy, direct-child peak RSS measurement, repeated runs,
   dataset and prediction SHA-256 hashes, runtime variability, and normalized
   accuracy metrics.
6. A reproducible analysis script and notebook/report package for the external
   benchmark.

The detailed MVP command and output inventory remains in `docs/mvp_status.md`.
Algorithm descriptions and benchmark protocol are in `docs/algorithms.md` and
`docs/benchmark_plan.md`.

## Validation completed

All checks below completed successfully before this handoff:

- Python: 170 pytest tests passed in 60.19 seconds.
- Rust: 4 tests passed.
- Rust `clippy` passed with warnings treated as errors.
- Rust formatting check passed.
- Python byte-code compilation passed.
- External comparison: all 45 runs completed successfully
  (5 datasets x 3 tools x 3 repetitions).
- Every normalized prediction file was deterministic across repetitions.
- Benchmark analysis package validation passed. HTML browser verification was
  structural/semantic only because a usable headless browser was unavailable;
  this is a report-rendering limitation, not a benchmark-computation failure.

Use the `tandemx-dev` conda environment for source tests. The external benchmark
also requires the documented TRF and TideHunter executables. Commands are in
`benchmarks/README.md`.

## External benchmark findings

The controlled comparison used the same existing FASTA inputs, one thread per
tool, three repetitions, and the median runtime. Tested versions were TandemX
0.1.0 release build, TRF 4.10.0-rc.2, and TideHunter 1.5.5.

On the 2.5-10 Mb benchmark inputs:

- TandemX was 1.84-4.19 times faster than TRF.
- TandemX was 1.52-4.78 times faster than TideHunter.
- TandemX used less peak memory than TideHunter on all five datasets.
- Memory relative to TRF was mixed: TRF used less on four small/medium inputs,
  while TandemX used 12.2% less on the 10 Mb input.
- TandemX and TRF both achieved macro recall 1.0, precision 1.0, false-positive
  rate 0, and monomer-length MAE 0 on these synthetic datasets.
- TideHunter achieved macro recall 1.0, precision 0.998405,
  false-positive rate 0.0048, and monomer-length MAE 0.000533.

These results support a scoped publication claim: TandemX has a throughput
advantage on the tested data while integrating read-based abundance, assembly
representation, and probe prioritization in one reproducible workflow. They do
not support a claim of universal accuracy or memory superiority over every tandem
repeat finder.

Local reproducible results are under:

```text
benchmarks/results/external_tool_comparison_20260715/
├── summary.tsv
├── raw_runs.tsv
├── dataset_manifest.tsv
├── environment.json
└── analysis/
    ├── pairwise_comparison.tsv
    ├── macro_accuracy.tsv
    ├── analysis_validation.json
    ├── tandemx_external_tool_analysis.ipynb
    └── report.html
```

## Data boundary and unresolved evidence gap

The repository contains a previous result summary for the black-rye analysis
(164 families and 568,482 supporting reads), but no corresponding raw black-rye
FASTQ/FASTA was found during the latest audit. Therefore the real-data comparison
cannot yet be rerun from raw reads, and those summary counts must not be presented
as an independently reproduced external-tool benchmark.

If the raw data become available, preserve them outside Git, record checksums and
provenance, and run the same fixed-input, fixed-thread, repeated benchmark design.

## Recommended next work

1. Complete GitHub publication of the current branch and keep this file updated
   with the resulting pull-request link.
2. Obtain or locate the raw real-data reads and run TandemX, TRF, TideHunter, and
   where practical TRASH on identical bounded subsets.
3. Add truth-controlled noisy simulations covering substitutions, indels,
   fragmented arrays, nested repeats, reverse complements, and negative controls.
4. Measure scaling beyond 10 Mb and add multi-process/chunked execution only after
   profiling identifies the dominant remaining bottlenecks.
5. Treat RepeatExplorer2/TAREAN as a workflow-level comparison rather than a
   drop-in per-read finder; document differences in inputs and biological claims.

## New-session rule

A new Codex session must not repeat completed optimization work merely because it
lacks chat history. It should first read this file, inspect the two implementation
commits, verify the current tests and artifacts, and continue from the recommended
next work or the user's newest instruction.
