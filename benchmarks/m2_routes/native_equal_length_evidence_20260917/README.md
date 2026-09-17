# Frozen C3 equal-length structural development score, 2026-09-17

The six-case challenge was frozen at commit `a1d3ad9` before scoring.
`benchmarks/controlled_collapse/native_equal_length_development_v1/protocol.json`
pins the source, generator, frozen M2 prototype, and prior C3 trim pilot
hashes. The new scorer verifies those hashes, all six edited FASTA hashes and
array hashes, and reconstructs the same three original native-read windows
from the unmodified FASTQ. Natural-flank windows must match the prior trim
receipt and primary PAF CIGAR coordinates exactly. The runner fails closed if
the output directory already exists.

## Reproduce

From the repository root in `tandemx-dev`, choose a new output path:

```bash
conda run -n tandemx-dev python -m benchmarks.m2_routes.native_equal_length_score \
  --output /tmp/tandemx_equal_length_replay_new
conda run -n tandemx-dev pytest -q tests/unit/test_m2_native_equal_length_score.py \
  tests/unit/test_m2_native_read_pilot.py tests/unit/test_m2_alignment_route.py
```

The archived `per_case.jsonl`, `summary.json`, `run_log.txt`, and `receipt.json`
are deterministic. A fresh run reproduced all four SHA-256 values exactly.
The scorer uses the frozen M2 `decompose` defaults, including its 4,096-bp
array cap and 250,000 candidate-alignment budget. The candidate catalogue is
the supplied operational B1 v3 C1/C2/C3 178-bp templates, not a native-read
discovery result. The comparator is the same median read-span baseline with
the fixed 5% or 2-bp threshold, evaluated on the same three eligible reads.
It is a **same-read stop-loss comparator**, not an equal-input accuracy peer:
M2 receives supplied 178-bp operational templates while the baseline receives
only anchored span lengths. All six edited arrays remain 3,560 bp, so a
span-only method has no signal.

## Complete case results

| Operation | Injected truth | M2 technical path | Span baseline | Operational label/orientation path |
| --- | --- | --- | --- | --- |
| intact | negative | SUPPORTED | SUPPORTED | same as native |
| invert_tile_05 | positive | DISCORDANT | SUPPORTED | changed; identifiable |
| invert_block_05_09 | positive | DISCORDANT | SUPPORTED | changed; identifiable |
| swap_adjacent_05_06 | positive | SUPPORTED | SUPPORTED | unchanged; unidentifiable |
| swap_distant_05_15 | positive | SUPPORTED | SUPPORTED | unchanged; unidentifiable |
| replace_tile_05_with_06 | positive | SUPPORTED | SUPPORTED | unchanged; unidentifiable |

Across the full intent-to-diagnose denominator, M2 detects **2/5** structural
positives, misses **3/5**, supports the intact negative (**1/1**), and abstains
**0/6**. The span baseline detects **0/5**, misses **5/5**, supports the intact
negative (**1/1**), and abstains **0/6**. M2 detects **2/2** positives whose
change is represented by its label/orientation path. The three missed
positives all resolve to the same 20-copy `C3+` path as the original reads:
the swap and tile-replacement edits are invisible to this **supplied-template
label/orientation representation**. The 2/2 subset is defined using M2's own
source and edited paths, so it is representation-conditioned and close to a
tautology given that the intact native paths match the source. It is neither
independent sensitivity nor evidence that those edits are biologically or
sequence-level unidentifiable. The primary result is **2/5** across all
positive operations. This diagnostic comparison does not establish equal-prior
superiority over the template-free length comparator. These are explicit
negative method results, not a reason to tune frozen thresholds after scoring.

This is one Col-CEN assembly lineage and the same three HiFi ZMWs from the
earlier deletion pilot. The natural 1,024-bp flanks were selected from this
development source after shorter flanks proved repetitive. Whole-genome
mapping uniqueness, same-plant donor/haplotype pairing, native structural
truth, and biological replication are unverified. `biological_audit_state` is
`NOT_EVALUATED_PAIRING_UNVERIFIED` for every case. These engineered edits
establish only injected sequence truth, not a native assembly error or
biological accuracy claim. This extension has no independent held-out test.
