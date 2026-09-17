# M1 identifiability and final abundance-novelty gate, 2026-09-17

**Decision: NO-GO for an accurate absolute family-abundance estimator as the
current TandemX Methods novelty claim.** This is a decision about the present
evidence and development program, not a proof that no future method can solve
the problem. Do not promote M1 A–F, the occupancy prototype, or a physical
copy-number interpretation into production or the manuscript. Keep the frozen
estimator's output as read-derived family evidence and an assembly-audit score,
with its existing uncertainty warnings. A later quantitative claim would need
a new independent truth design and prespecified calibration gate; no seventh
similar estimator is authorized by this decision.

## Executable evidence and status definitions

`benchmarks/m1_identifiability/diagnostics.py` computes circular 5-mer
catalogue signatures, pairwise weighted shared-word fractions, the complete
signature matrix rank/condition number, each family's fraction of exclusive
word positions, and exact nullspace involvement. Known-source development
reads are independently scored by the same best cyclic-Hamming 18-mismatch
gate used in the frozen tournament, retaining rejection, ties, correct unique
assignment and wrong unique assignment against all source reads. This is a
diagnostic, **not** another abundance estimator.

The per-family state has an explicitly narrow mathematical meaning:

- `NON_IDENTIFIABLE`: the family's coefficient participates in the exact
  nullspace of the supplied error-free catalogue's aggregate 5-mer matrix;
  individual abundance cannot be recovered from those features.
- `PARTIALLY_IDENTIFIABLE`: no exact nullspace involvement, but no catalogue
  5-mer occurs exclusively in that family. Recovery depends on combinations
  of shared features and may be sensitive to model error.
- `IDENTIFIABLE`: full per-family rank and at least one exclusive catalogue
  5-mer. This says **only** that noiseless supplied-catalogue features can
  distinguish it. It is neither calibrated read assignment nor physical
  copy-number identifiability. The matrix condition number and observed read
  assignment remain separate diagnostics; no observed error was used to tune
  a status cutoff.

The reproducible command is `conda run -n tandemx-dev python -m
benchmarks.m1_identifiability.run`. The full per-case output, input identities,
source hashes, all 18 previously formal simulated catalogues and limitations
are in `benchmarks/m1_identifiability/evidence_20260917/results.json` (SHA-256
`27c215d8a537bb4b6cb08f8d387a1973425dc55eded29f8c5fc3d5026e3be26f`
at this checkpoint). Those 18 prior catalogues are **not** a new held-out
validation set. The 12 M1 cases share three toy seeds, not 12 donors.

## What failed, and what remains distinguishable

| Frozen condition, 3 seeds | Catalogue result | Known-source read result | Interpretation |
| --- | --- | --- | --- |
| Related units, 3% substitution error | All three columns full rank; condition number 3.70–3.99; positive families have exclusive-position fractions 0.062–0.188. | Of 1,500 positive reads, 1,496 are uniquely correct, 3 tie, 1 is uniquely wrong; none rejected. | The M1 A–E losses here cannot be attributed to exact catalogue non-identifiability. Their model/feature assumptions and background handling matter. |
| Related units, 12% error | Same catalogue rank/condition as the matching seeds. | 1,470/1,500 positive reads uniquely correct; 3 rejected, 20 tied, 7 uniquely wrong. | Higher errors degrade assignment without changing the catalogue's noiseless rank. |
| `f2` with extra 20% substitutions | Same catalogue rank/condition. | `f2`: 395/540 source reads rejected by the common gate, only 134/540 uniquely correct; `f1`: 952/960 uniquely correct. | Family-specific error and exposure/callability dominate this designed failure. No estimator working only on accepted reads can recover the excluded `f2` observations without an independently validated model for their exclusion. |
| Identical `f1`/`f2` | Rank 2 for 3 families; `f1` and `f2` are `NON_IDENTIFIABLE`, zero exclusive positions; decoy remains individually separable. | All 1,500 positive reads tie between `f1`/`f2`; none has a unique family label. | Individual `f1`/`f2` abundance must be refused or reported only as a combined identifiable group. Assigning a split is unsupported by this input. |

All 18 prior formal simulated catalogues have rank 3/3 and an exclusive 5-mer
for each of their three supplied truth units. Even the three
`shared_fragment` catalogues have full rank (condition 1.12–1.25). This does
**not** contradict their historical shared-read/background errors: simple
catalogue rank is a necessary observability check for this feature matrix,
not a sufficient assay of erroneous reads, genome background, variable
representatives or sampling.

The earlier unified A–F tournament found no point-estimator candidate that
beat ordinary competitive mapping in any of the three distinguishable
development conditions; B's bootstrap coverage was 2/18 and E's 14/18
correlated positive-family cells. The later occupancy prototype did not
consistently reduce error or wrong-family assignment. Those results, kept in
`docs/negative_results/`, are the direct performance basis for the NO-GO.

## Explicit unresolved causes and output boundary

The toy source labels give exact **read-source counts**, not physical array
copy numbers. They contain no genome-wide genomic background, donor-to-donor
coverage variance, library sampling/callability model, or experimentally
measured physical copies. A single monomer representative per family cannot
establish that within-family diversity, nested/HOR units, harmonics and
shared genomic sequence are sufficiently represented. Therefore this analysis
cannot rank those real-world causes, validate a quantitative confidence
interval, or claim a calibrated real-family `IDENTIFIABLE` status.

For a future scientific output gate, `NON_IDENTIFIABLE` objects must have
individual absolute abundance set to NA/unresolved, retaining only group
evidence if group identifiability is demonstrated. `PARTIALLY_IDENTIFIABLE`
objects require explicit ambiguity and no physical copy-number claim.
`IDENTIFIABLE` under this input-only diagnostic still requires independent
read-error, background, representative-diversity and sampling validation
before a calibrated quantitative result can be reported. **This checkpoint
does not change the frozen production TSV or retrofit these statuses to
unvalidated real catalogues.**
