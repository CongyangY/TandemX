# YSD56 real-sequence bp edit development case

This is an archived real assembly **sequence** with controlled in silico edits.
It is a development data-integrity test, not same-donor read–assembly validation,
physical array length truth, independent biological replication, or model
performance evidence. No T7 write or new download was made.

## Frozen source and qualification

- Archived source: `paper/evidence/legume_ysd56_candidate_v1/TXF000708_array.fasta`
- Source SHA-256: `0866dbf818cbcd81b7a768c974fb21f221b3de95cfd20e95b9314613b5b04df9`
- Extraction receipt SHA-256: `478233f22a0a3315833abdd01ff6af3df1caeca203db470b7d2838d5503b637c`
- Archive manifest SHA-256: `d0af413da2b5aff24a8d7bdb69391e3b429667c287da26af44cac0fa5fac26c5`
- Confirmation summary SHA-256: `9ad4aea00b73382764a8f559e13ece90153b4b2506e23bb8970e502cdc499ee3`
- Archived coordinate: `CP154579.1:[13966881,14090910)` (zero-based), 124,029 bp.
- TRF and TideHunter each reported full interval coverage and dominant 785-bp
  period in the archived confirmation summary. This supports tandem structure
  in the assembly. It does **not** provide exact per-copy boundaries.
- The fragment length modulo 785 is 784. Complete-unit count and residual-unit
  location therefore remain **unknown**. No 785-bp copy truth is assigned.
- The original full compressed assembly on T7 was not re-read or rehashed in
  this run; the archived extraction and its two local receipts were verified.

## Replay

Run from the repository root:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.controlled_bp_edit \
  --source paper/evidence/legume_ysd56_candidate_v1/TXF000708_array.fasta \
  --archive-receipt paper/evidence/legume_ysd56_candidate_v1/array_extraction_receipt.json \
  --archive-manifest paper/evidence/legume_ysd56_candidate_v1/MANIFEST.tsv \
  --confirmation-summary paper/evidence/legume_ysd56_candidate_v1/array_confirmation_summary.tsv \
  --outdir benchmarks/controlled_collapse/ysd56_bp_development_v1
conda run -n tandemx-dev pytest -q tests/unit/test_controlled_bp_edit.py tests/unit/test_controlled_collapse.py
```

Use a fresh output directory for replay. The editor buffers at most a
1-Mb nominated fragment and deletes a right suffix of `floor(length × p/100)`
bases for p=0,25,50,75,100. Each case has a FASTA, one-row TSV ledger, and
piecewise coordinate chain. `receipt.json` has source, extractor receipt,
manifest, confirmation, generator, and per-output SHA-256 values. Tests
independently reconstruct every output from the archived sequence and recheck
the coordinate chain and hashes.

The TSVs have headers and stable fields. In each ledger, `source_region_*0`
identifies the original genomic half-open interval;
`deleted_fragment_*0` identifies the deleted half-open interval within the
archived fragment; and `deleted_genomic_*0` adds the archive's genomic start.
`injected_deleted_bp` is the exact edit amount, `retained_bp` and
`edited_fragment_*0` describe the resulting fragment, and the sequence hashes
are SHA-256 of uppercase sequence without FASTA formatting. `period_context_bp`
is only a dominant-period observation. `complete_unit_count_truth=unknown`,
`unit_count_status=unknown_no_copy_boundaries`,
`truth_scope=injected_bp_delta_only`, and `read_pairing_status=not_evaluated`
are mandatory interpretation fields. Chain rows use zero-based, half-open
`source_*0` and `edited_*0` coordinates; `delete` maps a source interval to a
zero-width edited junction. All five cases, including 100% deletion, remain in
the ledger and receipt.

| Deleted percent | Injected deletion (bp) | Retained (bp) |
| ---: | ---: | ---: |
| 0 | 0 | 124,029 |
| 25 | 31,007 | 93,022 |
| 50 | 62,014 | 62,015 |
| 75 | 93,021 | 31,008 |
| 100 | 124,029 | 0 |

The 100% case has an empty FASTA sequence record. Some downstream tools may
reject it; that outcome must be recorded as an input/technical failure, not a
zero call. This fragment has no real flanks or paired read evidence in the
editing run. The only truth established is the exact injected bp delta and
the edited sequence/coordinate map.

The focused archive and scorer-boundary tests passed **17/17** on 2026-09-16.
They include independent source-prefix reconstruction for all five edited
FASTAs, coordinate and per-output hash checks, and a test that declared donor
metadata never becomes verified read–assembly evidence.
The complete repository regression then passed **927 tests** in `tandemx-dev`.
