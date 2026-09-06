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

The first Mo17 whole-library random sample (840 reads / 11,680,888 bp) completed
all three tools at T7 `results/Mo17_complete_random_11Mb_pilot_v1_20260906`.
Wall times were 258.829 s TandemX, 11.859 s TRF and 3.508 s TideHunter; respective
RSS values were 168.58, 194.72 and 147.80 MiB. This unfavorable TandemX timing is
retained. Read scanning took 14.514 s; the subsequent 599-family exhaustive audit
was the main bottleneck. Concurrent QC makes these diagnostic measurements.
Native/streamed audit optimization must be rerun against the same sample and all
six discovery data files checked for exact equality before reporting its effect.

That rerun completed: all six outputs byte-identical, TandemX 17.348 s / 81.70
MiB, TRF 12.019 s / 165.17 MiB, TideHunter 3.644 s / 155.03 MiB. The 111.506-Mb
sample (8,084 reads) then completed: TandemX 247.539 s / 165.47 MiB, TRF 87.493 s
/ 212.52 MiB, TideHunter 34.115 s / 324.36 MiB. Single executions are engineering
diagnostics; the input subsets are nested. All values and source snapshots are
retained at `paper/evidence/Mo17_real_pilot` and the T7 result directories.

At 111 Mb, 4,219 representatives caused 8,897,871 exhaustive pair rows. For the
optional exact-related audit, add `--family-audit related` to the pilot command.
The runner records this non-default policy in the command and environment.
Compare five primary data files byte-for-byte and compare the related table to
the full table's non-distinct rows. Do not demand full table byte equality after
explicitly changing its row policy, or hide the omitted-distinct-pair count.
