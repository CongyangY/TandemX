# Isolated clustering resource replay

The representative index now stores each `(representative ID, circular-word
multiplicity)` in one unsigned 64-bit array item: 32 bits per nonnegative ID
and positive count. Counts/IDs outside these fields fail explicitly; they are
never truncated. Keys, multiplicities, insertion order, thresholds and alignment
rules are unchanged. The index is released before family/member outputs are
built. This reduces object overhead in principle; measured whole-pipeline memory
can still be dominated by other stages and must be reported independently.

```bash
python -m benchmarks.scripts.replay_clustering_isolated \
  --candidate-run /path/to/completed/1Gb_real_comparison \
  --baseline-run /path/to/completed/unpacked_index_comparison \
  --outdir /path/to/new/isolated_replay --timeout 3600
pytest -q tests/unit/test_monomer_clustering.py tests/integration/test_isolated_clustering_replay.py
```

The baseline run supplies its original source snapshot. Every core Python/Rust/
native file other than `clustering.py` must match the current source; changed
snapshots or additional algorithm changes fail preflight. The candidate run
supplies the identical `candidate_reads.tsv` and `candidate_monomers.fa` for both
workers. Candidate IDs/cardinalities are checked. These TSV scores have the
precision of serialized exports, not the original live in-memory candidates.

Each implementation runs in a fresh child with a frozen source/worker snapshot,
one thread and one fixed-order repetition. `clustering_seconds` times only the
clustering call. Child wall/RSS includes interpreter imports, candidate loading,
clustering and streamed JSON output. Candidate catalogues are resident; whole
reads/genomes are not loaded. Both complete family/member JSON payloads remain
on disk, and SHA-256 must agree exactly. Failed/time-limited attempts and missing
results remain explicit. No JSON payloads are merged into the parent's memory.

This stage diagnostic separates source representation changes from input
changes. It does not establish full discovery-pipeline parity, three-repetition
isolated timing, external-tool superiority or new biological inference.

## Native exact candidate gate

The Rust backend now keeps compact `(u32 ID, u32 multiplicity)` postings behind
an injective encoding of 1–9-base ACGTN keys. Query overlap accumulation happens
in Rust, with thresholds cached by representative length and only touched sums
reset. The Python gate remains the reference. Input words must be unique and
their counts must sum to sequence length; malformed appends fail before index
mutation. Order, rounding, N treatment, exhaustive paths, alignment, tie-breaking
and family rules are unchanged. Rebuild an older extension before using Rust;
missing native capabilities fail explicitly.

```bash
python -m benchmarks.scripts.replay_clustering_isolated \
  --candidate-run /path/to/completed/1Gb_real_comparison \
  --native-index-ablation --outdir /path/to/new/native_gate_replay
```

This separate mode freezes one current source/native snapshot for both children.
Only the Python-index worker disables creation of the native index, invoking the
retained Python gate while preserving native alignment. Receipts identify this
controlled child-only ablation. It does not relax the historical source-comparison
preflight: old baselines with different native or other core files still fail.
Full payload hashes must agree. Actual resource benefits require completed runs.

The earlier 85,663-candidate Python storage replay is archived under
`paper/evidence/Mo17_compact_clustering`: 31.44% lower child peak RSS but 15.86%
longer clustering time, exact complete-output parity. That run predates the
native gate; its trade-off must not be attributed to the new implementation.

The completed native ablation is in `paper/evidence/Mo17_native_index`: identical
outputs,218.015→164.312 clustering seconds,250.22→259.63 MiB child peak. Time
decreased while memory increased in this diagnostic; no universal win is claimed.

## Sequence-native index interface

The Rust path accepts each complete canonicalized monomer sequence and computes
the same circular canonical ACGTN word multiplicities inside the native index.
This removes the former per-candidate Python string `Counter` and Python-to-Rust
word-list conversion. The Python word interface remains available as a frozen
reference adapter; candidate thresholds, representative order, alignments and
membership outputs are unchanged.

```bash
python -m benchmarks.scripts.replay_clustering_isolated \
  --candidate-run /path/to/completed/real_comparison \
  --native-interface-ablation \
  --outdir /path/to/new/interface_replay --timeout 3600
```

The two fresh children are labelled `word_bridge` and `sequence_native`.
`validation.json` must report complete execution and exact payload-hash parity.
The fixed-order single repetition is an engineering diagnostic; publication
timing still requires shuffled repeated runs under isolated acquisition load.
After a valid replay, compact evidence can be copied without the duplicate full
payloads using `python -m benchmarks.scripts.archive_clustering_replay_evidence`;
the archiver rechecks payload parity and every frozen source hash before copying.

## Full live-pipeline replay and profiling

```bash
python -m benchmarks.scripts.replay_discovery \
  --previous-run /path/to/completed/real_comparison \
  --outdir /path/to/new/full_replay --threads 4
pytest -q tests/integration/test_discovery_replay.py
```

This command checks the baseline input hash and successful discovery receipt,
freezes current source, preserves original CLI settings except an explicitly
requested thread budget, and compares all seven deterministic products using live
candidate scores. The environment records baseline and replay thread counts.
Changed inputs or any product disagreement fail explicitly. Optional `--profile`
identifies hotspots but adds overhead, so profiled resources cannot be used as
unprofiled speed rankings.
