# TandemX Benchmarks

## Expanded evaluation

The comparator contract now includes ULTRA and plans task-matched SRF, TRASH,
TAREAN and specialized annotation/probe comparisons; see
[comparator_matrix.md](../docs/comparator_matrix.md). Adding a method to that
contract does not mean it has been benchmarked.

`elastic_development_v1.yaml` runs the experimental elastic algorithm on the
same development distributions without external dependencies. Both first and
corrected development attempts are retained on T7. Use the actual source digest
in `environment.json`, not the config filename alone, to identify an experiment.

ULTRA is supported by the strict challenge adapter. Build the official source
locally (no system install is required), set its executable under `tools.ultra`,
and keep that source revision and build commands in provenance. Optional
`ultra_options` are `window_size` and `windows` (positive integers) and `tune`
and `tune_indel` (booleans). This permits bounded buffers and a separately
reported tuned run. Model defaults are otherwise retained. The same maximum
period/minimum span apply, and minimum-period filtering is applied to outputs.
ULTRA's documented TSV coordinates are already zero-based half-open; `*` in
consensus becomes N and `.` is missing consensus, never invented sequence.
Tuning occurs inside the timed command and its cost is retained.

The first ULTRA feasibility pilot uses only ten reads per scenario; it must not
be compared directly with 100-read timing or memory rows from the main challenge.
Its bounded 5000-base window/one-window queue are documented settings based on
the observable 5000-base input-read length, not repeat truth. The upstream default
at max period 1000 advertises about 4.09 GB of buffer memory; both settings and
actual measured peak memory must be distinguished.

## TideCluster assembly comparator

The pinned TideCluster 1.21.2 comparator image is defined in
`benchmarks/containers/tidecluster/`. It rebuilds the official TideHunter 1.4.3
commit and replaces MMseqs2 with the official 16-747c6 SSE2 release after exact
checksum/version checks, because both conda binaries otherwise execute illegal
instructions under linux/amd64 emulation on Apple Silicon. Build and use this
container separately from the read-local TideHunter comparator.

Normalize a completed two-stage TideCluster run against explicit assembly truth:

```bash
python -m benchmarks.tidecluster.normalize \
  --tidehunter-gff RUN/tc_tidehunter.gff3 \
  --clustering-gff RUN/tc_clustering.gff3 \
  --assembly assembly.fa --truth truth.tsv --catalogue catalogue.fa \
  --outdir RUN/normalized
```

`archive_tidecluster_smoke.py` then verifies both profiled stages, every input
and normalized-output hash, the planted-truth evidence boundary, internal GNU
time and the exact Docker image ID. The current 199.1-kb smoke recovered 3/3
families and arrays with 0-bp period MAE, 29-bp boundary MAE and 0.998247 base
union precision. Its clustering maximum RSS was 7,669,232 kB. These values prove
that the comparator path works; they are not publication-scale ranks.

The first cascade promotion configuration predeclared all 16 scenarios, held-out
seeds 3101–3103 and three repetitions. That once-only split is complete and must
not be rerun. It failed two gates: nine TRF low-complexity controls timed out and
the TandemX/TideHunter runtime geometric-mean ratio was 2.457943, above 2.0.
All failures and missing measurements remain in the archived evidence.

## Cascade gap-free speed development

Fresh development seed 1201 profiles the same 16 scenarios with TandemX and
TideHunter. The baseline completed 96/96 commands and had a paired runtime
geometric-mean ratio of 2.357698. An unguarded 30%-span/95%-identity fast path
reached 1.933332 but worsened the 0.1%-indel boundary MAE from 1.364 to 14.221 bp;
it is retained as a rejected experiment. The selected candidate adds >=95%
valid shifted columns and <=2% residual from a whole number of units. It also
completed 96/96 commands, reached 1.989148, retained array/family endpoints,
kept minimum positive base-union F1 at 0.997597 and maximum boundary MAE at
2.921 bp, and made no negative-control calls.

`analyze_cascade_fast_path.py` records 1,600 read-level observable screen rows
and scores proposed rules against development truth. The selected rule accepted
403 reads with no truth-scored negative, wrong-single-array or multi-array
acceptance. `archive_cascade_gap_free_development.py` preserves the baseline,
rejected attempt, selected candidate, profiles and audit at
`paper/evidence/cascade_gap_free_development_v1`.

`cascade_gap_free_validation_v1.yaml` freezes validation seed 2201, the selected
rule provenance and 14 accuracy/resource gates. It may be executed exactly once
only after the source/config commit and both hosted workflows pass:

```bash
python -m benchmarks.challenge.run \
  --config benchmarks/configs/cascade_gap_free_validation_v1.yaml \
  --split validation \
  --outdir /Volumes/T7/Codex/TandemX/results/cascade_gap_free_validation_v1_20260907

python -m benchmarks.scripts.evaluate_cascade_heldout \
  --config benchmarks/configs/cascade_gap_free_validation_v1.yaml \
  --run /Volumes/T7/Codex/TandemX/results/cascade_gap_free_validation_v1_20260907 \
  --outdir /Volumes/T7/Codex/TandemX/results/cascade_gap_free_validation_v1_20260907/gate_evaluation
```

The evaluator rejects incomplete/duplicated matrices, config-hash mismatches,
missing metrics and comparator failures. A failed speed, memory or accuracy gate
must be retained as a failed validation experiment. Seed 2201 cannot tune a
replacement model after it is observed; seeds 3201–3203 remain reserved.

## Quantify depth-gated validation

`factorial_scale_quantify_validation_v1.json` predeclares untouched validation
genomes 6401--6403 under the same factorial process as development. The source/
config freeze passed hosted CI before they were generated with
`generate_factorial_scale --split validation`. All three receipts report
`validation_used=true`; a separate audit verified 30 manifests and all 93
declared payload files. Seeds 6401--6403 are consumed, while 7401--7403 remain
refused reserved seeds.

`evaluate_quantify_depth_gated_validation.py` runs the public total-bases and
empirical-control modes, then applies the already fixed mean-control-depth >=2
condition-level switch. Its nine gates, exact one-time command and output fields
are documented in [quantify_calibration.md](../docs/quantify_calibration.md) and
[file_formats.md](../docs/file_formats.md#depth-gated-quantify-validation).
The once-only run completed 54/54 commands and passed all nine gates. The optional
public `quantify --single-copy-min-depth 2` mode implements the validated routing
rule without changing the default controls behavior. No replacement rule may be
tuned against 6401--6403.

## Challenge benchmark and public-data pilot

Development release smoke check (inside `tandemx-dev`):

```bash
maturin build --release --interpreter python --out dist
python benchmarks/scripts/verify_wheel.py \
  --wheel dist/<wheel-for-your-platform>.whl --outdir /tmp/tandemx-wheel-verification
```

The verification installs without dependencies into a temporary target, confirms
Python and Rust import from that target, and runs all seven toy workflow steps.
It writes per-command logs and `wheel_validation.json` containing wheel SHA-256,
interpreter, commands/exit codes, import locations, validated steps and scope.
Platform-specific wheel validation does not establish cross-platform readiness.

Run development experiments from the repository root inside `tandemx-dev`:

```bash
python -m benchmarks.challenge.run \
  --config benchmarks/configs/challenge_v1.yaml \
  --split development \
  --outdir /Volumes/T7/Codex/TandemX/results/challenge_development
```

Choose any empty writable output directory on another machine. Use `--scenarios
clean_171 indel_1pct two_arrays` for a smaller experiment. TRF/TideHunter executable
paths are configurable; the checked-in paths use the pre-existing local comparison
environment. TandemX resolves from the active `tandemx-dev` PATH. No truth file is
passed to any finder. The 16 scenarios include substitutions, insertions/deletions,
divergent units, related families, short/multiple arrays and negative controls.
The zero-result CLI behavior was fixed after the first development baseline;
retain the baseline's nine failed TandemX negative-control runs as failures.

`development`, `validation` and `heldout` use disjoint seeds. Once a seed or family
is used for tuning it is no longer held-out evidence. All three tools receive the
same FASTA and period limits. Normalized output uses 0-based half-open coordinates
and a common minimum span. TRF's overlapping/harmonic predictions remain visible
and unmatched duplicate calls count against raw array precision. This is **raw
call precision**, not a post-merged nonredundant annotation score.

The runner records subprocess commands, input/source/executable hashes, timeout
and exit status, direct-child peak RSS, interval matches and final-catalog family
recovery. Run failures or malformed outputs receive `NA`, never a fabricated zero
recall. Timing repetitions are not biological replicates. Use isolated, frozen
runs for publication runtime comparisons. The first development run was used
alongside other development checks, so its timings are exploratory.

After completion, render the four-panel diagnostic with:

```bash
python -m benchmarks.scripts.plot_challenge_benchmark \
  --run /Volumes/T7/Codex/TandemX/results/challenge_development \
  --outdir /Volumes/T7/Codex/TandemX/results/challenge_development/figures
```

Its PDF/SVG/PNG, panel source tables and checksum receipt are developmental
diagnostics. The figure does not claim a final benchmark or general superiority.

Retrieve a bounded, public Arabidopsis pilot (Python 3.11 as in `tandemx-dev`):

```bash
python benchmarks/scripts/fetch_ena_subset.py \
  --accession ERR6210723 --read-count 1000 \
  --outdir /Volumes/T7/Codex/TandemX/data/raw/ERR6210723_prefix1000
```

The downloader limits compressed bytes, validates complete FASTQ records, writes
FASTA, and saves ENA run/sample metadata plus a subset SHA-256 receipt. It reads
a prefix, not a random library sample; the remote full-file MD5 and gzip trailer
are explicitly **not verified** for this bounded extraction. Do not use the
prefix to claim unbiased whole-genome abundance. See [file formats](../docs/file_formats.md)
and [release programme](../docs/release_program.md).

This directory contains synthetic benchmark configuration and runner scripts for measuring the toy-scale TandemX MVP before any real large-genome analysis.

## External Tool Comparison

The truth-aware comparison harness runs TandemX, TRF, and TideHunter on the
same single-family synthetic FASTA datasets using one thread per tool:

```bash
.conda-benchmark/bin/python benchmarks/scripts/run_external_tool_comparison.py \
  --config benchmarks/configs/external_tool_comparison.yaml \
  --outdir benchmarks/results/external_tool_comparison
```

To rerun tools against existing benchmark inputs without rewriting them, use:

```bash
.conda-benchmark/bin/python benchmarks/scripts/run_external_tool_comparison.py \
  --config benchmarks/configs/external_tool_comparison.yaml \
  --datasets-dir benchmarks/results/external_tool_comparison/datasets \
  --outdir benchmarks/results/external_tool_comparison_20260715
```

It writes `raw_runs.tsv`, `summary.tsv`, `dataset_manifest.tsv`, generated
FASTA/truth files when `--datasets-dir` is omitted, an environment manifest,
and per-run logs. Runtime and direct-process peak resident memory are recorded
with `wait4`; three-run medians, runtime CV, and normalized-prediction digests
are included in the summary. TandemX and TideHunter are explicitly restricted
to one thread, while TRF is single-process. Accuracy is evaluated at read level
with a period tolerance of `max(2 bp, 2% of truth period)`. Simulated errors
are substitutions only, so these are engineering comparisons rather than a
complete model of HiFi or ONT error profiles. TRASH and TAREAN are excluded
from the numeric chart because their primary assembly and short-read graph
tasks are not directly equivalent to per-read detection.

All parsed predictions are filtered to the configured period range before
scoring. This is necessary because TRF exposes a maximum-period argument but
does not expose the same minimum-period restriction as the other two tools.

Validate and derive publication-review tables from a completed comparison with:

```bash
python benchmarks/scripts/analyze_external_tool_comparison.py \
  --summary benchmarks/results/external_tool_comparison_20260715/summary.tsv \
  --dataset-manifest benchmarks/results/external_tool_comparison_20260715/dataset_manifest.tsv \
  --outdir benchmarks/results/external_tool_comparison_20260715/analysis
```

The analysis step writes pairwise speed/memory comparisons, macro accuracy,
the largest-dataset chart source, a bounded dataset table, and a validation
receipt. Its report companion notebook uses only the Python standard library.

The tested macOS ARM64 environment lives at `.conda-benchmark/` and is
described by `benchmarks/environment.external-tools.yml`. TideHunter must be
installed from its official ARM64 release archive because the Bioconda package
is unavailable for this platform. The downloaded v1.5.6 archive currently
contains a binary that reports version 1.5.5; this discrepancy is retained in
`tool_versions.tsv` rather than silently normalized.

The benchmark workflow is:

```text
simulate -> discover -> quantify -> locate -> probe -> validate
```

The default discovery step remains de novo. Simulator truth files are used only after the run to calculate benchmark accuracy summaries; they are never passed as analysis command inputs.

## Files

```text
benchmarks/configs/synthetic_scale.yaml
benchmarks/scripts/inspect_reads.py
benchmarks/scripts/run_real_read_pilot_benchmark.py
benchmarks/scripts/run_synthetic_benchmark.py
benchmarks/results/.gitkeep
benchmarks/simulated/.gitkeep
```

## Synthetic Scales

`synthetic_scale.yaml` defines:

1. `tiny`: 1,000 reads and the only scale used by pytest.
2. `small`: 10,000 reads for manual runtime checks.
3. `pilot`: 50,000 reads for manual subset scaling.
4. `real_pilot_manual`: 100,000 reads, never run by default.

Only `tiny` is intended for pytest. Larger scales remain manual because the Python backend is still single-process.

## Real-read Pilot

Use the real-read runner only for bounded engineering pilots. It runs `discover -> validate`, never reads truth files, and writes `tmpfq_benchmark_summary.tsv` with processed reads/bases, runtime, throughput, candidate rate, recovered family count, validation status, and the exact command. Peak memory is recorded as `NA` until a portable reporter is available.

```bash
python benchmarks/scripts/run_real_read_pilot_benchmark.py \
  --reads reads.fastq.gz \
  --max-reads 1000,5000,10000,25000 \
  --kmer-backend rust \
  --outdir /tmp/tandemx_real_pilot
```

The summary includes the selected backend. Run Python and Rust into separate output directories before computing speedups; do not benchmark a debug-mode Rust build.

## Run

```bash
python benchmarks/scripts/run_synthetic_benchmark.py \
  --config benchmarks/configs/synthetic_scale.yaml \
  --scale tiny \
  --outdir /tmp/tandemx_benchmark_tiny
```

Outputs:

1. `benchmark_summary.tsv`
2. `accuracy_summary.tsv`
3. `<scale>/logs/*.stdout.log`
4. `<scale>/logs/*.stderr.log`
5. `<scale>/simulated`, `<scale>/discover`, `<scale>/quantify`, `<scale>/locate`, `<scale>/probe`

`benchmark_summary.tsv` fields:

1. `benchmark_id`
2. `scale`
3. `seed`
4. `read_count`
5. `read_length`
6. `total_read_bp`
7. `monomer_lengths`
8. `command`
9. `runtime_seconds`
10. `exit_status`
11. `output_validated`
12. `recovered_family_count`
13. `processed_reads`
14. `processed_bases`
15. `candidate_reads`
16. `candidates_per_mb`
17. `reads_per_second`
18. `mb_per_second`
19. `peak_memory_mb`
20. `algorithm_mode`
21. `notes`

`accuracy_summary.tsv` fields:

1. `benchmark_id`
2. `expected_monomer_length`
3. `recovered_closest_length`
4. `length_error_bp`
5. `expected_read_copy_bp`
6. `estimated_read_copy_bp`
7. `copy_number_relative_error`
8. `locate_status`
9. `recovered_sequence_identity`
10. `matching_method`
11. `notes`

## Runtime And Memory

The runner records wall-clock runtime and per-process peak resident memory with the Unix `wait4` resource record. `peak_memory_mb` is reported directly in `benchmark_summary.tsv`. Accuracy pairing is length-aware but sequence-driven, so unrelated repeats with the same monomer length are not treated as the same truth family.

## Interpretation

Synthetic benchmark results are engineering signals for runtime, output validity and toy accuracy. They do not validate TandemX for 7-20 Gb plant genome production analysis.

## Sequence clustering and independent endpoint audit

`sequence_clustering_v1.yaml` fixes elastic discovery, operational sequence
clustering at 95% similarity, and three technical repetitions. Development and
validation seeds are reused explicitly for diagnosing earlier errors; they are
not fresh held-out evidence. The preset does not redefine planted truth labels.

The benchmark now adds independent edlib cyclic global edit similarity and
period-independent union-of-interval coverage. Strict equal-length sequence
recovery and raw call precision remain visible. This distinguishes small length
errors from missing sequence, and redundant calls from wrong genomic bases.
The benchmark extra and environment pin edlib 1.3.9.post1; the detector does not
use it. Example re-scoring of an archived completed suite:

```bash
python -m benchmarks.scripts.rescore_challenge --runs /path/to/completed-run \
  --outdir /path/to/new-rescore-directory
```

The rescorer selects first-repetition archived predictions and preserves failure
states, configurations and input hashes. It does not update original outputs or
supply new tool timings. `cyclic_monomer_recall` uses TandemX's final catalogue
and other tools' per-array consensuses; these are different output granularities.
Candidate FASTA evidence now permits a further matched per-array consensus study.

## SRF catalogue and abundance workflow

The pilot runner is `python -m benchmarks.scripts.run_srf_pilot`; it executes seven
native stages rather than timing only graph construction. All original failures
are retained. A disclosed guard identifies empty successful KMC dumps and skips
SRF's empty-input assertion without creating a fake native FASTA. The task and
platform scope, parameters and measurements are in [docs/srf_workflow.md](../docs/srf_workflow.md).
An SRF output may be a higher-order repeat; out-of-period motifs remain preserved.
Minimal-monomer recovery alone does not establish absence of its constituent units.
