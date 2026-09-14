# YSD56 TandemX--TideHunter permissive family postprocess

Status: `blocked_missing_post_remount_source_outputs`.

The remount removed the intended source directory
`/Volumes/T7/Codex/TandemX/results/soybean_ysd56_30x_comparison_v1_20260911`.
Consequently, this archive intentionally reports `NA` for every support count;
it contains no inferred exact, permissive, unmatched, or TandemX-only family
number. See `receipt.json` for the exact missing files and the SHA-256 of the
replay script.

## Input contract after source restoration

The restored run directory must retain these unmodified outputs from the frozen
real-30x runner:

1. `normalization.sqlite`, containing normalized `tandemx` and `tidehunter`
   representatives, their canonical primitive-period keys, and per-tool
   `chunk_count` values;
2. `execution_summary.tsv`, with every TideHunter row marked `completed` and
   `normalization=ok`; and
3. `run_manifest.json`, retaining the frozen input, tool and rule receipt.

The permissive population is fixed to TandemX families with
`chunk_count=3` that have no exact normalized TideHunter key. Exact support is
recorded separately. TideHunter representatives are not additionally required
to recur across all three partitions by the frozen support screen; their
`chunk_count` is reported for every permissive match. Any absent, failed,
timed-out or non-normalized TideHunter partition remains a
technical unresolved state. It is never interpreted as no TideHunter support.

## One-command replay

Run from the repository root, choosing a new output directory:

```bash
PYTHONDONTWRITEBYTECODE=1 conda run -n tandemx-dev python \
  benchmarks/scripts/postprocess_ysd56_tidehunter_permissive.py \
  --run-dir /Volumes/T7/Codex/TandemX/results/soybean_ysd56_30x_comparison_v1_20260911 \
  --outdir paper/evidence/ysd56_tidehunter_permissive_postprocess_replay_v1
```

The script first samples rotation- and strand-aware 4-mer posting sizes as a
diagnostic. It does not use a seed index to make a final negative call. Before
alignment it exactly estimates the number of all length-eligible pairs. If that
count exceeds the declared alignment budget, every nonexact recurrent family is
written as technical unresolved and no unmatched count is reported. If it is
within budget, a temporary local SQLite index retrieves every representative at
an eligible length and edlib `HW` aligns the complete shorter sequence against a
doubled longer sequence in both orientations. A pair supports the comparator
only at either 0.90 glocal identity with
`shorter/longer >= 0.90`, or 0.80 identity under the frozen integer-period-
multiple rule with relative length error at most 0.05.

An unmatched output row means only that the frozen TideHunter postprocess found
no qualifying TideHunter representative. It is explicitly not a TandemX-only,
novelty, biological-recurrence, or full-comparator conclusion.
