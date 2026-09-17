# Hostile audit of the C3 equal-length structural development challenge

The C3 challenge, generator, six FASTAs, exact-edit receipt and protocol were
committed as `a1d3ad9` before this challenge was scored. The scorer and native
results were committed as `d74bef1`. That freeze followed
the earlier C3 deletion pilot on the **same** assembly lineage, three original
HiFi molecules and supplied operational 178-bp monomer catalogue. C3 is thus
a bounded development extension prompted by a known length-baseline weakness,
not an untouched validation cohort.

I independently checked the frozen protocol/receipt digests, generator and M2
prototype hashes, source config, original native FASTQ, and each of the six
FASTA hashes. All six edited contexts have 13,560 bp, preserve both 5-kb
natural flanks, and retain the 3,560-bp array interval `[5000,8560)`. Direct
reconstruction from the unedited twenty 178-bp tiles reproduced all six
arrays: one intact, a single-tile reverse complement, a four-tile reversed
and reverse-complemented block, two swaps, and a same-length replacement.
The five positive edits change respectively 138, 540, 16, 24 and 8 aligned
array bases while adding/removing zero bases. These are **injected assembly
edits**, not independently observed native-array events or monomer-copy truth.

The scorer uses the same three frozen, source-guided, natural-flank-trimmed C3
read intervals for every case and retains the prior PAF/trim checks. It calls
the unchanged M2 alignment prototype with its default parameters. The span
baseline receives those same read lengths and a fixed 5% or 2-bp threshold
from the earlier C3 pilot. M2 additionally receives supplied operational
monomer templates and the edited assembly sequence, while the baseline uses
lengths alone. This is a same-read **stop-loss baseline** with less prior
information, not an equal-prior peer accuracy ranking. All six edited arrays
and the median read span are 3,560 bp; the span difference is zero and the
threshold is 178 bp in each case.

I checked every hash in the committed score receipt against the current
protocol, challenge receipt, source FASTQ/context/PAF, prior trim files, M2
prototype, scorer and result files. A fresh `tandemx-dev` replay to
`/private/tmp/c3_hostile_replay_20260917` produced byte-identical
`per_case.jsonl`, `summary.json` and `run_log.txt` (matching SHA-256 values
`2dc29de8…`, `b2d16d5e…` and `c23c7e4c…`). Focused generator and scorer
tests passed **5/5**. No frozen input or M2 core was edited in this audit.

| Full frozen denominator | M2 technical label+orientation path | Length baseline |
| --- | ---: | ---: |
| Engineered positives detected | 2/5 | 0/5 |
| Engineered positives called supported | 3/5 | 5/5 |
| Positive abstentions | 0/5 | 0/5 |
| Intact negative supported | 1/1 | 1/1 |

The two M2 detections are **exactly** the single-tile and four-tile
inversions. The intact source and all three native read paths resolve as 20
`C3+` labels. The inversion assemblies resolve with respectively one and four
`C3-` labels. Both swaps and the replacement still resolve as twenty `C3+`
labels, so M2's present label+orientation endpoint cannot distinguish them.
The scorer keeps all five positives in the primary denominator. Its 2/2
“path-identifiable” subset is defined using the same M2 decomposition on the
source and edited assembly, then compared to native read decompositions.
It is an **operational representation-conditioned description**, close to a
tautology when intact reads match the source path. It is not an independent
sensitivity estimate and does not show that the three missed edits are
unobservable from raw sequence or by other methods.

The unchanged molecules come from a pooled sample. Their local C3 context
and flanks support technical comparison, but donor/haplotype pairing to the
assembly, genome-wide unique placement, native tile boundaries, physical
copy truth and biological replication are unverified. The result is a
technical response to engineered equal-length edits only. It cannot be used
as biological collapse accuracy, a calibrated copy-number result, or a new
Methods superiority claim.
