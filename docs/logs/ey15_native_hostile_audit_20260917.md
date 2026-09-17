# Ey15 native interval and span score: independent hostile audit (2026-09-17)

## Scope and reproducibility

I reviewed the source package, the nine exact-bp assembly edits, the frozen
`score_protocol.json` in commit `e4e7760`, the context and whole-assembly PAFs,
the scorer, and the emitted case records. This is a **development** audit of one
selected Ey15-2 (9994) interval. It does not promote the interval to a held-out
benchmark or establish biological copy truth.

The source manifest and original seven-read FASTQ match their archived SHA-256
values. The source package checked the current 1.125-GB `sample_003` gzip
completely; its parent 18.6-GB run was not reread. The exact-edit verifier
reconstructed all nine FASTAs from source context plus recorded deletion
coordinates and checked every resulting hash. A clean scorer replay produced
the same `read_trims.json`, `per_case.jsonl`, and `summary.json` byte for byte.

I independently parsed all seven full-context **primary** PAF CIGARs, required
complete CIGAR syntax and exact query/reference consumption, and projected the
array boundaries. The projections agree exactly with the flank-derived read
boundaries in all seven cases. Three additional local PAF rows are short
repeated-sequence hits and do not span the full context. Each selected read has
one aligner-reported primary hit on the full five-chromosome assembly, at the
claimed Chr1 locus, with MAPQ 60. This is evidence of locus anchoring within
that assembly, subject to aligner reporting and assembly representation.

The seven source read IDs have distinct ZMW suffixes. The 1,024-bp left and
right natural flanks each have one best near-exact read hit at the fixed 5%
edit-distance limit. An additional check of nonoverlapping read segments found
no second hit within that limit for any of the 14 flank searches. The observed
array spans are 3,244, 3,240, 3,246, 3,252, 3,244, 3,241 and 3,248 bp;
their median is 3,244 bp. These are seven molecules from one material, not
seven independent donors or biological replicates.

## Denominator and scoring audit

The denominator is **nine edited assembly cases**: eight injected deletions
and one intact sequence. All cases reuse one 3,244-bp source interval and the
same seven original reads. The read count cannot be multiplied by case count
to create independent observations. The smallest deletion is 811 bp, while the
predeclared intact 5% discordance threshold is 163 bp. This is a large-edit
feasibility check; it does not test behavior near the decision boundary.

The scorer now measures each edited assembly's span directly between its two
natural-flank hits on the edited FASTA, then checks that measurement against
the exact-edit ledger. It does not rely on the ledger's `edited_array_bp` as
its measured predictor. I added a regression test that changes an in-memory
edited FASTA span by 20 bp while leaving its committed receipt unchanged: the
scorer must reject the disagreement. The test passed.

Replayed classifications are eight discordant injected-deletion cases (8/8),
one supported intact case (1/1), and zero abstentions under this selected,
source-guided diagnostic. These labels reflect exact **injected assembly bp
changes**, not naturally missing bp or physical copy count. The diagnostic
uses an already chosen interval and its natural flanks. It is not an
unlocalized, genome-wide read-discovery evaluation or a read-to-each-edited-
assembly mapping benchmark.

## Leakage and scientific limits

The interval was selected from 16 prior TandemX located arrays by maximum
direct period-shift identity before read mapping in this audit. The seven reads
were then selected for source-interval spanning from a seeded 50,267-read
subset. The score's apparent success is consequently conditional on this
high-periodicity, read-supported interval. Although the edited FASTAs and
score protocol were committed before scoring, the case construction is not
independent of the chosen source interval and its known 3,244-bp span.

The assembly and reads are attributed to the same published sample 9994;
exact DNA extraction and haplotype pairing remain unverified. The newer HiFi
assembly shares the technology and sample lineage with the reads, so their
agreement cannot serve as independent physical copy-number truth. Whole-
assembly mapping does not prove the absence of an unrepresented alternative
haplotype. The operational family period is 420 bp, above the frozen M2 route's
300-bp cap; **M2 evaluated 0/9 cases**. Monomer boundaries, physical copy
number, biological collapse sensitivity, and method superiority remain
unevaluated. The score supplies no TandemX runtime or memory benchmark.

## Audit outcome

The exact-edit and source-guided span classifications reproduce within their
stated technical scope. I found no residual error that invalidates the 8/8
and 1/1 large-edit diagnostic after the edited-FASTA measurement correction.
The evidence must remain labeled development-only, one array/one material,
injected-bp truth only, M2-ineligible, and biological accuracy not evaluated.

Focused regression check:

```text
conda run -n tandemx-dev pytest -q tests/unit/test_score_ey15_native_span.py
1 passed
```
