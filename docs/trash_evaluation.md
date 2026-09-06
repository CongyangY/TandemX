# TRASH assembly evaluation

Use the pinned [container recipe](../benchmarks/containers/trash/README.md) and
the same independently generated assembly for both versions. Author-example
installation controls and read-based comparisons are distinct experiments.

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.run_trash_container \
  --fasta /path/to/factorial/genome/genome.fa --tool trash \
  --image YOUR_RECORDED_IMMUTABLE_IMAGE_ID --outdir /path/to/run --timeout 7200
conda run -n tandemx-dev python -m benchmarks.scripts.score_trash_factorial \
  --genome /path/to/factorial/genome --run /path/to/run --outdir /path/to/evaluation
conda run -n tandemx-dev pytest -q tests/unit/test_trash_adapters.py \
  tests/unit/test_score_trash_factorial.py
```

`--tool trash2` runs the successor separately. Default de novo settings receive
no known templates or truth. One CPU and8 GiB container memory are allocated.
Timeouts, OOM, startup and native failures remain missing accuracy observations.

The evaluator reports one-to-one array/period matches, conditional boundary and
period errors, independent cyclic sequence recovery and interval-union coverage.
It separately evaluates native broad regions and actual monomer intervals; a
wide discovery window is not equivalent to claiming every contained base as a
monomer. Both region/coverage endpoints must accompany any accuracy comparison.
All55 planted families, including the1.026-Mb array, remain in the denominator
of the current10-Mb input. IID background and a single genome seed are development
conditions; they do not establish performance on real plant assemblies.

For TRASH1, pinned source `f7d53a0` creates window positions with
`start=(ii-1)*window; end=start+window-1`, then calls one-based R `str_sub` with
those values (`src/fn_repeat_identifier.R`). Both interpretations are reported.
`src/fn_extract_all_repeats.R` adds region.start to local one-based positions;
the author control shows a systematic−1bp monomer offset. This is independently
audited by exact reference-sequence extraction for each new input, not assumed
to apply to every region or used to alter original output. Start-zero windows
are a particular reason to retain the sensitivity analysis.

TRASH2 source `f290a5e` uses one-based inclusive windows and subtracts1 from
the regional start before coordinate adjustment. Its author control agrees with
that interpretation. Native peak spacing and representative-sequence length
remain separate period policies for both versions. No policy is chosen after
looking at which gives a preferred ranking.

See [field definitions](file_formats.md#factorial-assembly-comparator-evaluation).
Synthetic parser fixtures in tests verify scoring; they are never reported as
executed comparator data. Resource comparisons need identical input and platform;
Linux container RSS cannot be directly ranked against a different-input macOS
read run. Concurrent diagnostics are excluded from final isolated timing claims.
