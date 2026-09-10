# Context-adversarial toy fixture contract

**Status:** design specification only. This document defines a small development
negative-control fixture for context-sensitive repeat discovery and quantification.
It contains no measured result, benchmark claim, new holdout, or algorithm change.

## Fixed sequence design

Use a deterministic seed recorded in the fixture manifest and a fixed 120-bp
target monomer (`M120`). The exact sequence should be generated once by the
fixture builder and recorded with its digest. A 24-bp substring (`B24`) is a
fixed contiguous slice of `M120` and is also allowed to occur as its own
24-bp tandem repeat in background. The background sequence must be long enough
to test context, but remain toy-scale and below 1 MB.

Every generated read carries a stable identifier and a source-case label. The
fixture manifest records the seed, sequence digests, read lengths, mutation
parameters, and expected truth intervals.

## Required source cases

The builder must emit these cases independently and in the stated combinations:

1. **Complete array:** one complete six-copy `M120` array (720 bp), with truth
   intervals for all six copies.
2. **Partial high-abundance background:** many copies of only `B24`, producing a
   high-abundance local tandem signal without a complete 120-bp family array.
   This background may itself be a valid 24-bp tandem repeat; it must not be
   labelled as a complete `M120` family by fixture truth.
3. **Coexistence:** the complete six-copy `M120` array and the `B24` background
   occur in the same input collection.
4. **Shared-context families:** two distinct 120-bp monomers, `M120_A` and
   `M120_B`, share a declared internal sequence segment. Their complete arrays
   and shared segment locations are recorded separately.
5. **Non-overlapping read context:** emit reads in which the two near-related
   families do not overlap within the same read. This tests whether a method
   incorrectly infers a relationship from catalogue-level similarity alone.
6. **Reverse complement:** include reverse-complement representations of complete
   arrays and background segments, with strand truth retained.
7. **Substitution stress:** include a declared 1% substitution condition, using
   a deterministic mutation mask and preserving the pre-mutation truth interval.
8. **Indel stress:** include a declared 1% indel condition, with insertion and
   deletion positions/mutation events recorded separately from truth coordinates.
9. **Invalid/empty inputs:** include an `N`-containing read and an empty read as
   explicit input records. Their validation fate must be recorded, not silently
   dropped.

Cases 7 and 8 may be applied to cases 1–5, but the manifest must identify the
parent case and mutation condition for every read. Reverse-complement status is
also an explicit field, not inferred from the algorithm output.

## Ground truth and interpretation boundary

Truth is generated from the fixture construction process: expected repeat
intervals, copy identity, family identity, strand, mutation parent and source
case are written before any TandemX or comparator algorithm runs. Truth must
never be reconstructed from candidate intervals, family assignments, abundance
estimates or other algorithm output.

The 24-bp background is a deliberate partial-context adversary. It can be a
real 24-bp tandem repeat, but it is not a complete 120-bp `M120` family. A
method may report it as a candidate or reject it; either outcome is measured
against the declared truth and does not change the fixture label.

Shared segments between `M120_A` and `M120_B` are context ambiguity, not proof of
homology, higher-order organization, or an evolutionary relationship. Reports
must use candidate/sequence-similarity language and retain an unresolved state
where the evidence cannot distinguish families.

## Fair metrics

Report metrics separately by source case, strand, mutation condition and family:

* interval precision, recall and base-union precision/recall against generated
  repeat intervals;
* complete-array recovery and copy-count error for the six-copy array;
* false-positive calls of `M120` on the `B24`-only background;
* family assignment confusion for `M120_A`, `M120_B`, shared segments and
  unresolved calls;
* reverse-complement invariance (equivalent calls after strand canonicalization);
* sensitivity and false-positive rate under 1% substitutions and 1% indels;
* abundance error only where the fixture supplies a declared abundance truth,
  with zero/absent estimates represented explicitly rather than coerced;
* validation fates for the `N` and empty reads, including clear error categories;
* runtime and peak memory only as descriptive toy execution measurements.

Do not collapse these into one accuracy number. A partial array, a family
confusion, an invalid read, and an unavailable measurement are different fates.
Report denominators for every metric and retain `NA` when a metric is undefined.

## Expected failure modes to preserve

The fixture is intended to expose, not hide, these possible failures:

* high-copy `B24` context may attract partial-array candidates or inflate a
  120-bp abundance estimate;
* shared sequence may produce ambiguous family assignment when reads do not
  contain distinguishing context;
* substitutions or indels may cause a true array to be rejected or fragmented;
* reverse-complement handling may lose strand equivalence;
* an empty or `N`-containing read may trigger validation errors or an explicit
  no-call;
* short or partial evidence may be insufficient for a complete-array claim.

All failures, warnings, rejected reads, no-calls, unresolved relationships and
resource-limited executions remain in the fixture output and receipt. They must
not be removed to make the toy example appear complete. The fixture therefore
serves as a reproducible development negative control and does not establish
performance on real genomes.

## Observed development boundary

In the fixed seeded two-family toy case, separating the families with 20 `G`
bases produced a prototype interval beginning at position 379 rather than the
expected position 380. This is a known one-base background inclusion at the
context boundary. It is retained as a regression observation and must not be
silently removed by changing a threshold or by replacing the case with an `N`
separator. An `N`-separator case remains separately useful for clean mass-
conservation and invalid-context tests. This observation is not a biological
result or a performance claim.
