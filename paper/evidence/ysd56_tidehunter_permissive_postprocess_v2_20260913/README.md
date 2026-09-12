# YSD56 recurrent-family TideHunter permissive preflight

Input was the internal recovery copy
`tmp/ysd56_formal_run_recovery_20260913/`, not the mounted T7 volume. Its
`normalization.sqlite`, `execution_summary.tsv`, and `run_manifest.json`
SHA-256 values are recorded in `receipt.json` and match the copied source
receipt.

The frozen population contains 12,070 TandemX families recurring in all three
formal partitions. Of these, 2,099 have exact normalized TideHunter support;
the permissive screen population is the remaining 9,971 recurrent exact-
unmatched families. TideHunter recurrence is not imposed as an extra criterion,
because that was not part of the frozen support screen; matched TideHunter
`chunk_count` would be reported if alignment were feasible.

The exact length-only preflight finds 791,390,883 direct-length pairs and
2,738,327,814 integer-multiple pairs, for 3,529,718,697 unique
length-eligible edlib alignments. This exceeds the predeclared 5,000,000-pair
budget. The 9,971 nonexact recurrent families are therefore technical
unresolved. `permissive_tidehunter_supported_family_count` and unmatched counts
are `NA`; no TandemX-only claim is made.

The q=4 multiset necessary-condition diagnostic retained 27,416 of 35,526
length-eligible sample pairs (77.2%). It was not used for final matching. Even
if adopted after a separate proof-and-regression gate, this sampled retention
does not make the 3.53-billion-pair problem operationally feasible. A future
comparator-only acceleration would need a lossless candidate-index design and
an independently verified no-false-negative test suite; it is not a production
TandemX algorithm change.
