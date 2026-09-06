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
