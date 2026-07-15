# TandemX Benchmarks

This directory contains synthetic benchmark configuration and runner scripts for measuring the toy-scale TandemX MVP before any real large-genome analysis.

## External Tool Comparison

The truth-aware comparison harness runs TandemX, TRF, and TideHunter on the
same single-family synthetic FASTA datasets using one thread per tool:

```bash
.conda-benchmark/bin/python benchmarks/scripts/run_external_tool_comparison.py \
  --config benchmarks/configs/external_tool_comparison.yaml \
  --outdir benchmarks/results/external_tool_comparison
```

It writes `raw_runs.tsv`, `summary.tsv`, generated FASTA/truth files, an
environment manifest, and per-run logs. Accuracy is evaluated at read level
with a period tolerance of `max(2 bp, 2% of truth period)`. Simulated errors
are substitutions only, so these are engineering comparisons rather than a
complete model of HiFi or ONT error profiles. TRASH and TAREAN are excluded
from the numeric chart because their primary assembly and short-read graph
tasks are not directly equivalent to per-read detection.

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
