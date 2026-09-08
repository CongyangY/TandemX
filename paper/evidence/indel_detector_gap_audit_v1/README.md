# Indel-rich detector gap audit v1

This audit answers a narrow prioritization question from the already consumed,
frozen seed-2201 guarded-cascade validation: is TandemX still clearly behind
TideHunter in **read-local accuracy** on the three planted indel scenarios?

The answer is no for this frozen synthetic distribution. Both tools had array
recall 1.0 at 0.1%, 1%, and 4% total indels. TandemX array precision was 1.0 in
all three cases; TideHunter precision was 1.0, 1.0, and 0.958904. Boundary MAE
was mixed and remained below 1.7 bp for TandemX. TandemX was, however, 3.34--
3.80 times slower in the three indel scenarios while using 0.32--0.43 times
TideHunter's direct-child peak RSS.

Therefore this evidence does not justify replacing the validated detector with
a new seed-and-chain algorithm as the next scientific priority. The current
priority remains donor-matched collapse validation and broader task-matched
comparators. This is a deferral, not a claim that seed-and-chain cannot help:
real-read interval truth, biological replicates, error profiles beyond the
frozen simulation, and per-condition speed remain unresolved.

`audit.json` records every extracted indel row, the exact source hashes, and the
decision boundary. No seed was rerun and no validation threshold was changed.
