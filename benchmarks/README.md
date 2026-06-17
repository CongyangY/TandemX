# TandemX Benchmarks

This directory contains synthetic benchmark configuration and runner scripts for measuring the toy-scale TandemX MVP before any real large-genome analysis.

The benchmark workflow is:

```text
simulate -> discover -> quantify -> locate -> probe -> validate
```

The default discovery step remains de novo. Simulator truth files are used only after the run to calculate benchmark accuracy summaries; they are never passed as analysis command inputs.

## Files

```text
benchmarks/configs/synthetic_scale.yaml
benchmarks/scripts/run_synthetic_benchmark.py
benchmarks/results/.gitkeep
benchmarks/simulated/.gitkeep
```

## Synthetic Scales

`synthetic_scale.yaml` defines:

1. `tiny`: CI-scale, kept below 1 MB of reads.
2. `small`: manual synthetic run around 10 MB of reads.
3. `medium`: manual synthetic run around 100 MB of reads.
4. `large`: manual stress run around 500 MB of reads.

Only `tiny` is intended for pytest. `small`, `medium` and `large` are manual scale tests and may be slow with the current Python toy algorithms.

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
13. `notes`

`accuracy_summary.tsv` fields:

1. `benchmark_id`
2. `expected_monomer_length`
3. `recovered_closest_length`
4. `length_error_bp`
5. `expected_read_copy_bp`
6. `estimated_read_copy_bp`
7. `copy_number_relative_error`
8. `locate_status`
9. `notes`

## Runtime And Memory

The runner records wall-clock runtime with Python's standard library. Peak memory is not recorded portably yet; use `/usr/bin/time -v` on Linux or the platform `time` command on macOS for manual memory measurement until TandemX has a cross-platform resource reporter.

## Interpretation

Synthetic benchmark results are engineering signals for runtime, output validity and toy accuracy. They do not validate TandemX for 7-20 Gb plant genome production analysis.
