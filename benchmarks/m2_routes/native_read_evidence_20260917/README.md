# C3 native-read M2 feasibility pilot, 2026-09-17

This is a **research-only development pilot**, not biological accuracy
validation. It uses the existing Col-CEN C3 assembly context, three original
HiFi molecules from one pooled Col-0N study, and the nine exact local C3 edit
cases. The original context is also reported as a duplicate technical reference
row. The C3 `terminal_000` row is the one intact negative in the nine-case
ledger. The other eight rows contain injected deletions. No native FASTQ base,
quality, or source file was changed.

## Reproduction and provenance

Run from the repository root:

```bash
conda run -n tandemx-dev python -m benchmarks.m2_routes.native_read_pilot \
  --output /tmp/tandemx_native_m2_replay_new
conda run -n tandemx-dev pytest -q tests/unit/test_m2_native_read_pilot.py \
  tests/unit/test_m2_alignment_route.py
```

`receipt.json` gives exact SHA-256 values for the frozen source context,
unaltered original FASTQ, original PAF, source and edit receipts, frozen B1 v3
candidate-monomer bundle, frozen M2 prototype, and the two output files. The
fixed route uses its unchanged 4,096-bp input cap, 250,000-alignment budget,
15% monomer edit allowance, 8% length allowance, and two-edit label margin.
The source reads were 14,716, 18,270, and 16,693 bp; only the original
between-flank subsequence enters M2.

The script orients each read by the frozen primary PAF alignment, searches for
two **natural** 1,024-bp assembly flanks at up to 5% edit cost, requires one
best hit per flank, and checks order and <=10-bp agreement with the primary PAF
CIGAR projected array boundaries. The observed trim coordinates agree **exactly**
with PAF for all three reads: `[4769,8324)`, `[3967,7527)`, and
`[3280,6840)` in reference orientation, respectively. Their array spans are
3,555, 3,560, and 3,560 bp. `read_trims.json` records read IDs, source hashes,
strand, edit costs, and both coordinate estimates. The 1,024-bp flank length
was chosen **after inspecting this development material**: shorter 64 to
512-bp natural flanks had multiple optimal hits within the C3 reads. It is a
source-guided technical trim, not independent anchoring or a prospective
parameter validation. Locus discrimination was tested only among three selected
contexts; whole-genome uniqueness and identical donor/haplotype pairing remain
unverified.

M2 receives the **supplied B1 v3 C1/C2/C3 operational 178-bp monomer
templates**, not a de novo native-read catalogue. Their biological monomer
boundaries remain unverified. `per_case.jsonl` contains each assembly and read
decomposition, case-level technical comparison, and the simple read-span
length comparator on the same three eligible reads. The comparator uses the
existing frozen baseline's 5% or 2-bp threshold, whichever is larger.
Applying it to this source-guided trim is a post-hoc development replay. The
runner refuses to overwrite an existing output path. Neither decision is an
assembly-error or native missing-copy call.

## Observed result

| Nine edited C3 cases | Frozen M2 technical path | Simple read-span length |
| --- | ---: | ---: |
| Injected deletions detected | 8/8 | 8/8 |
| Intact `terminal_000` true negative | 1/1 | 1/1 |
| False positives on the intact negative | 0/1 | 0/1 |
| Case-level abstentions | 0/9 | 0/9 |

The extra `C3_source_original` row also compares as supported, but duplicates
the intact assembly sequence and is **excluded** from the nine-case denominator.
All three native intervals decomposed to 20 operational C3 labels. The edited
assemblies decompose into shorter all-C3 paths, so the path method and the
simple span method make the same binary decisions. This pilot shows **no
incremental M2 value over length**. It does not test edit localization or
order-specific discrimination, and eight related edits in one assembly lineage
are not eight independent biological observations.

`biological_audit_state` is `NOT_EVALUATED_PAIRING_UNVERIFIED` for every row.
The strict M2 `audit(..., pairing_status="verified")` route was **never called**.
An exact controlled edit establishes only the injected delta against this
source assembly; it cannot establish that the original read and assembly came
from the same plant or that a native array was collapsed. B2 remains a
previously consumed synthetic hold-out and was not touched or re-scored here.
