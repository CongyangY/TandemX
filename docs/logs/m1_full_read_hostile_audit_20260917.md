# Independent hostile audit: M1 full-read development comparison

## Verdict and scope

This is a bounded **development feasibility** result. The candidate and the
18-read challenge are toy scale and share one generated 80-bp catalogue. Here
“full-read” means that the whole input record is accepted and every base is
accounted for. Inference still classifies disjoint 40-bp cores using 40-bp
local flanks; it does not jointly optimize a path across the whole read. The
result supports only complete-record coverage with explicit abstentions on
these synthetic reads. It does not establish calibrated genomic copy number,
superiority to production `quantify`, real-read throughput, or assembly-deficit
accuracy. No public backend promotion is justified.

The protocol, candidate, generator, input FASTA, catalogue and truth ledger
were committed together in `9346de2` before the scored run. Results and runner
were committed as `81be93d`; `04f8929` clarified the local-tile wording in
the run log and README. The scoring runner was written after the freeze;
it records its own source hash. This distinction
prevents calling the whole scoring implementation prospectively frozen. I
reviewed `full_read_protocol_v1.json`, `generate_full_read_dev.py`,
`full_read_research.py`, `occupancy_research.py`, `model.py`,
`run_full_read_dev.py`, the worker receipts and the native `copy_number.tsv`.

## Audit findings

| Check | Finding | Consequence |
| --- | --- | --- |
| Input and truth leakage | Candidate inference reads the same FASTA and supplied catalogue as the two comparison workers. Its `classify_read`/`summarize_records` path does not read `truth.json`. Truth is read by the runner after worker execution for scoring. The candidate uses catalogue-only near-family grouping. | No observed direct label leakage. The designed toy and method share 80-bp unit assumptions; this is development, not an independent validation set. |
| Denominator and rejects | Three positive families (`f1`, `f2`, `twin_a`) remain in the unweighted per-family MARE denominator. An unassigned identical twin contributes 100% individual-family error. Across all 11,161 input bp, candidate mass is 7,120 assigned + 1,200 ambiguous + 2,841 unknown; chunked mapper mass is 4,400 + 1,200 + 5,561. | Abstention is not removed from the error denominator. MARE can still mask cancellation of false and missed assignments, so wrong-family and background errors must accompany it. |
| Truth-overlap specificity | Independent overlap of worker traces with the truth ledger gives candidate background unknown 2,755/2,758 bp and positive unknown 86 bp; chunked mapper background unknown 2,758/2,758 bp and positive unknown 2,803 bp. Candidate wrong assignment is 7 bp, including 3 background bp; mapper wrong assignment is 4 bp. | Lower candidate MARE (0.3410 versus 0.5928) is a toy occupancy observation, not unqualified better attribution. Candidate has more erroneous assigned bp. |
| Shared families | `twin_a` and `twin_b` catalogue sequences are exactly equal. All 1,200 truth `twin_a` bp are ambiguous; both individual estimates stay zero rather than being split by truth. | Individual abundance is non-identifiable here. Group evidence cannot be reported as either family's physical copy number. |
| Baseline and output units | The ordinary mapper's existing `align_and_gate` rejects full reads with lengths unlike its 80-bp reference; the runner therefore applies it to nonoverlapping 80-bp chunks, with short tails unknown. Production `quantify` reports genomic-bp estimates with artificial depth 1, whereas both research rows report assigned read bp. | The chunked mapper is a diagnostic adapter. Production is an input-compatible process check, not a common-unit accuracy baseline; its MARE and per-base assignment errors are correctly N/A. |
| Timing and RSS | The candidate timer/RSS stop after a streaming summary pass; its per-base trace is generated afterward. The chunked mapper times its trace, and production times output writing. Catalogue loading and process imports also differ. macOS `ru_maxrss` is whole-process peak while the wall timers cover narrower operations. | Raw worker values are archival only. Equal-input does not mean equal-work resource measurement; speed and peak-memory rankings are N/A. |
| Thresholds and biological scope | Candidate settings (40-bp core/flank, 0.225 edit fraction, 2-edit margin), 80-bp catalogue and scenario seed were in the pre-score commit. No post-score threshold change or reserved final holdout is visible. Three read repetitions per scenario share the same unit templates. | This is one related synthetic lineage with technical error replicates, not 18 independent biological tests. No real assay bias or physical-copy calibration is tested. |

The initial results draft omitted the preregistered `unknown_background_bp`
metric and omitted the imported `occupancy_research.py` source hash. These were
reported to the run owner before finalization. Commit `81be93d` added
`unknown_background_bp`, `unknown_positive_bp`, both missing source hashes,
and `resource_comparison_status=not_rankable_unequal_timing_scopes`. The frozen
protocol still says `equal_input_wall_seconds`; the separate status qualifies
that phrase without rewriting the pre-score record.

I independently recalculated all 18 truth intervals against both per-base
worker traces. Both methods cover every read contiguously and conserve all
11,161 bp. Their assigned-family totals, status mass, and three-family MARE
match the final `results.json`. The source and input SHA-256 values in that
receipt match the committed files, including the protocol, catalogue, reads,
truth, `occupancy_research.py` and the sequence parser. The production MARE is
`null`, and the raw production TSV is retained. Focused `tandemx-dev` tests:
`3 passed` for `tests/unit/test_m1_full_read_research.py`. This verification
does not turn the toy design into independent validation.

## Decision boundary

The method can proceed only as a separate research prototype. A future
production decision needs an independent donor or lineage test, externally
calibrated physical copy-bp units, an eligible full-read baseline with matched
prior information, explicit identical-family refusal, wrong-attribution and
zero-decoy reporting, and equivalent timed work at realistic read/catalogue
scale. Continuous abundance deficit and binary collapse calls require separate
validation. The present development MARE must not reopen the earlier M1 Gate A
NO-GO or be used as a manuscript superiority claim.
