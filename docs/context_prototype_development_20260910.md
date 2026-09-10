# Periodic-context abundance experiment: first development result

This is a mechanism prototype, not a new public command or validated estimator.
It tests whether full-unit periodic context can prevent shared local sequence
from being counted as the abundance of an entire operational repeat family.
Seed matching and chaining are established methods; the proposed research
question concerns their use as an eligibility condition for abundance evidence.
See `context_method_prior_art_20260910.md` for the closest methods.

## Algorithm and current limits

The frozen catalogue supplies each family's period and circular seed positions,
on both strands. Exact read seeds vote for a phase diagonal `(read_start -
template_phase) mod period`. Nearby hits on a diagonal form a chain. The initial
development rule requires two unit lengths, at least 70% of distinct template
start phases, k=11 and a maximum gap of 22 bp between seed starts. These are
explicit prototype parameters, not tuned or calibrated production defaults.
Accepted intervals are unioned within family. A sweep of interval endpoints
separates unique family bases, cross-family ambiguity and unassigned bases.
Different nonoverlapping families in the same read are retained.

The counting unit is a context-qualified read span. It is not a genome-wide
copy count, and the per-family eligible count is not an upper confidence bound.
Unassigned bases can contain missed true repeats. The prototype uses exact
phase diagonals: insertions/deletions can split them and cause undercounting.
Its index stores all seed phases; runtime is O(N*k + H), where H is the number
of seed/index matches, and per-read working memory is O(H). Highly repetitive
templates can make H large. There is no current speed or memory advantage claim.
The production quantifier does not call this code.

## Preserved first experiment

Run from the repository under `tandemx-dev`:

```bash
python -m benchmarks.scripts.evaluate_context_development --output NEW_RESULT.json
```

The script refuses an existing result file. It creates one deterministic 120-bp
unit and explicitly labels every scenario as development. Truth is fixed by
construction; all results, including adverse cases, appear in
`paper/evidence/context_development_v1/result.json`, together with configuration,
unit/input identities and source hashes. No held-out dataset was consumed.

| Development condition | Planted family bp | Raw diagnostic-depth bp | Context-qualified bp |
| --- | ---: | ---: | ---: |
| Six complete units | 720 | 720 | 720 |
| Repeated shared 24-bp fragment and background | 0 | 0 | 0 |
| Repeated shared 80-bp fragment and background | 0 | 7,200 | 0 |
| Complete array plus 80-bp fragment background | 720 | 7,920 | 720 |
| Seven substitutions in 720 bp | 720 | 600 | 720 |
| Seven insertions in 720 bp | 727 | 660 | 419 |
| One complete unit, below the context gate | 120 | 120 | 0 |

The raw-depth comparison is a mechanistic diagnostic: median circular-k-mer
counts times unit length. It is not a complete benchmark against the frozen
multi-k/quality/depth-normalized production estimator. Neither SRF nor ordinary
minimap2 occupancy was run here, so the table does not demonstrate superiority
over those methods. It also cannot attribute TXF000695's empirical instability
to this synthetic mechanism.

The result is informative in both directions: full-unit context rejects the
constructed shared-fragment contamination, but the strict phase model loses
308 of 727 planted bp with insertions and rejects single units. A separate
unit test preserves a one-base flank over-extension. The current prototype
therefore fails the requirement for general abundance replacement.

## Next decision gate

Retain this exact-phase result as the initial ablation. Before widening the
algorithm, compare with ordinary competitive mapping occupancy on the same
fixed catalogue, including background errors, repeat recall and unassigned
fraction. An indel-tolerant bounded phase transition is a candidate extension,
but it must outperform the exact-phase baseline without admitting fragment
decoys; it is not yet implemented. Multi-representative changes, softmax
assignment and parameter searches remain outside this first experiment.
Only after development performance and resource caps are established should
fresh family/genome validation be registered. The initial novelty assessment
remains conditional, not a claim that this prototype is the first of its kind.

## Result fields

Each JSON row records `scenario`, generated `input_sha256`, `read_bp`, planted
`truth_repeat_bp`, `raw_diagnostic_depth_bp`, its absolute error,
`context_qualified_bp`, its absolute error, `ambiguous_bp` and `unassigned_bp`.
All bp are sampled-read sequence quantities. Error is absolute bp difference,
not relative error at zero truth. The single-family experiment has no
ambiguous assignments; multi-family ambiguity is tested separately in unit tests.
