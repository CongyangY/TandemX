# Real-input comparator development pilot

Use the identical FASTA converted from a verified whole-library random FASTQ
sample. Original read IDs and sequence content are retained; only qualities are
omitted so every tool sees the same format. The source checksum, selected read
counts/bases, sampling receipt, FASTA checksum, source snapshot and tool hashes
are recorded. This controller is capped at 100,000 reads because normalization
keeps a read-length dictionary in memory; final larger studies need disk-backed
evaluation. It is not the final performance runner.

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.run_real_comparators \
  --sampling-receipt /path/to/samples/sampling_receipt.json \
  --sample-id sample_001 --outdir /path/to/new/pilot \
  --trf /path/to/trf --tidehunter /path/to/TideHunter --timeout 900
pytest -q tests/unit/test_real_comparator_inputs.py
```

This first workflow runs elastic/sequence-clustering TandemX, TRF and TideHunter
at one thread, one repetition and fixed shuffled order. The common evaluation
range is 30–1000 bp period and at least 100 bp called span. Native outputs are
preserved, including out-of-range calls. This does not test long rDNA units or
arbitrary higher-order structures, nor replace broader SRF/ULTRA/TRASH and other
task-matched comparisons in the cohort programme.

`summary.tsv` records `tool`, process `exit_code`, `timed_out`, wall
`runtime_seconds`, direct-child `peak_rss_mib`, user/system CPU seconds,
`observed_in_scope_calls`, `observed_positive_reads`, `observed_union_bp`,
`observed_union_base_fraction`, `normalization`, and `warning`. Union coverage
merges overlapping calls within each read. A successful zero-call result is
valid; failed execution/normalization gives missing observation fields. No
accuracy or superiority is inferred from call counts or agreement between tools.

Each tool retains its command, execution receipt, native output/logs and
`normalized_arrays.tsv` (`read_id`, 0-based half-open `start/end`, `period`,
native `sequence` when available, `family_id` when available). TandemX candidate
TSV alone has no consensus; its candidate FASTA and catalogue remain in native
outputs. Wall/RSS here are engineering diagnostics: downloads can overlap and
controller memory is excluded. Publication timings require isolated repeated
runs, full aggregate resource accounting, fixed tuning budgets and independent
accuracy evidence.
