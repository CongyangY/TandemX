# Exact Python target-counting opportunity (2026-09-10)

## Finding

`quantify_toy_copy_number` spends its read scan in target-only counting, but
all three Python scan paths construct a canonical k-mer string for every valid
window before testing target membership:

* `tandemx/quantify/mvp.py:536` limited/single-file scan;
* `tandemx/quantify/mvp.py:686` quality-correction scan;
* `tandemx/quantify/mvp.py:793` per-file worker scan.

The prior iterator, `iter_linear_canonical_kmers`, slices and reverse-complements
one Python string per window. The target counter is already sparse; it does not
need a string for non-target windows. The exact 2-bit canonical code primitive
in `tandemx/utils/kmers.py:81-132` provides the same canonical A/C/G/T key for
k=1..31. This is engineering only, not an estimator, threshold, or biological
change.

## Implemented exact contract

The Python path now maps only exactly matchable canonical target strings to
2-bit codes once per scan and tests rolling codes directly. It retains the old
string iterator whenever `k > 31` or the sequence is non-ASCII. Invalid or
noncanonical targets are omitted from the code map because they cannot match
legacy canonical output. `N` resets a rolling window in both paths. The Rust
counter, quality calculation, filters, batching, progress, and file-level
parallelism are unchanged.

The contract is byte equality of `copy_number.tsv` for identical Python-backend
inputs and config, plus equality of selected-k-mer Counters. Tests cover ACGT/N,
lowercase, non-ASCII fallback, k=31, k>31 fallback, and a full quantify run
against a monkeypatched legacy string oracle.

## Complexity and existing coverage

For fixed k, both paths are O(total read bases) time; the legacy string path
has O(total read bases × k) character work/allocation, while the rolling-code
path is O(total read bases) for k<=31. Both retain O(number of selected target
words) counter state. The change removes O(number of valid windows) transient
Python substring/reverse-complement allocation. `iter_canonical_kmer_codes`
now avoids low-complexity-window bookkeeping when `filter_low_complexity=False`;
that flag was already false for selected counting. Rust already performs a
native target-count scan, so the expected gain is limited to explicit Python
fallback users. No new parallelism, sketch, or sequence heuristic is proposed.

`docs/performance.md` records quantify as 32.68 s in the existing 100,000-read
profile, but does not attribute that time by function. Therefore the small
pilot below is a local engineering signal, not a predicted whole-library gain.

## Small reused measurement

Input (already consumed development material, not a holdout):
`/Volumes/T7/Codex/TandemX/results/abundance_baseline_v1_20260906/reads/s4101/c1_e0/reads.fa`
(196 KiB) with its s4101 catalogue. The reproducible driver is
`benchmarks/scripts/benchmark_python_target_counting.py`; its complete JSON
receipt is `paper/evidence/python_target_counting_benchmark_v1/result.json`. In
three fresh worker processes per variant, a complete Python-backend quantify
call gave byte-identical output SHA-256
`78c5d4afd7cde2def3ed9c3901f2ce7d31c1d806806b064366c0d661d193cf85`.

| Variant | Wall seconds, three runs | Process peak RSS bytes, three runs |
| --- | --- | --- |
| Legacy strings | 0.2626, 0.2684, 0.2631 | 24,428,544; 23,789,568; 23,805,952 |
| Rolling target codes | 0.0843, 0.0782, 0.0789 | 24,182,784; 23,920,640; 24,641,536 |

The input is small. The recorded timer starts immediately before
`quantify_toy_copy_number`, so process startup/import is excluded; the observed
roughly 3.34x median wall-time ratio is not a release performance claim. RSS overlap is
within ordinary process noise. A future performance record should run at least
three isolated processes on a fixed larger, already-consumed input, retain CPU,
Python/Rust build, command, file SHA-256 and `ru_maxrss` units, and compare
output bytes before reporting a speed or memory result.
