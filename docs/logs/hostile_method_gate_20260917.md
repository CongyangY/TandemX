# Independent hostile review of the M1/M2 method gate, 2026-09-17

Scope: review of committed M1 diagnostics, B1 v2/v3 truth and scoring, two
frozen M2 routes, the paired-flank baseline, and the assembly-only comparator.
I did not change prediction, truth, metric, or method files. The B1 v3 scores
were present before this review. This is a review of evidence, **not** a new
held-out or biological validation.

## Independent score audit

I read each raw `inputs.jsonl`, `truth.jsonl`, and route prediction JSONL,
keyed the rows by `case_id`, applied the frozen 0.5 threshold and non-`ok`
rule directly, then compared the counts with the three recorded
`summary.json` files. The recomputation did not call the B1 scorer.

| B1 v3 route | Input-only prediction SHA-256 | TP / FN / FP / TN / unresolved intact | `ok` / abstain | All-case sensitivity | Intact FPR |
| --- | --- | ---: | ---: | ---: | ---: |
| Monomer-path alignment | `49dcf68261c2c73157183457cdb5ef11e6aca54253f0247fb7cf8830b51c3279` | 24 / 9 / 0 / 6 / 0 | 36 / 3 | 24/33 = 0.7273 | 0/6 |
| Transition graph | `21e490113d9fd8af22238d6e27ef58ca966f3f88a0f0444deda3a722554b8039` | 0 / 33 / 0 / 0 / 6 | 0 / 39 | 0/33 | 0/6, with all negatives unresolved |
| Two-flank length baseline | `2232e2ccab0da61b42b203284314714a96b663d7c8ec291884030ed52b9836cc` | 18 / 15 / 0 / 6 / 0 | 39 / 0 | 18/33 = 0.5455 | 0/6 |

The rows and summaries agree. The graph's nominal FPR of zero must always be
shown with its **6/6 unresolved intact controls** and 100% rejection fraction;
it has no demonstrated specificity. Alignment's eligible-only 24/30 = 0.80
must not replace its intent-to-diagnose 24/33 = 0.7273. By event type,
alignment missed all three rearrangements, all three within-tile compressions,
and three of six duplications. It detected all three inversions, but these are
three edits of the same design on three related windows. The baseline missed
all three inversions and rearrangements, all three within-tile compressions,
all three single-tile deletions, and three of six duplications. Its frozen 5%
strict-inequality threshold puts a 178-bp loss of a 3,560-bp array exactly at
the boundary. This explains the single-tile misses; it does not make them
ignorable.

B1 v2 independently recomputes to 13 cases, 11 injected positives and two
intact controls. Alignment and graph each have 10 TP, one FN from abstention,
zero FP and two TN; the length baseline has nine TP, two FN, zero FP and two
TN. This confirms the earlier decision-count report, while the v2 frozen
scorer treats the route primary metric as null with a non-`ok` case. V2 has
only **10 unique edited assembly sequences** among 13 causal edits. Causal
subtype and breakpoint accuracy cannot be inferred for equivalent sequences.
V3 has 36 unique edited assemblies among 39 rows; the duplicated sequences
are the two intact controls per source window. V3 still has no independent
native monomer/HOR truth, and no route reports a complete breakpoint,
order, period, or calibrated confidence endpoint. Null values for these
metrics are appropriate, rather than evidence of correctness.

## Truth, leakage, and generalization limits

1. V3 uses three 3,560-bp Col-CEN v1.2 assembly intervals, all from **one
   assembly lineage**. The 39 rows are repeated edit conditions on these
   three source sequences. They are not 39 independent biological arrays or
   donors. The intervals were selected for strong 178-bp periodicity and
   exclude known assembly issue regions, so they are a favorable subset of
   native sequence. Sampling uncertainty and cross-genome transfer cannot be
   estimated from these rows.
2. Every case presents three different read IDs but the **same exact synthetic
   source-path sequence** under each ID. There are only three unique source
   read sequence groups across V3. These are neither independent molecules
   nor native HiFi/ONT errors. The native Col-0 read accession mentioned in
   project records has not been verified as the same interval, donor and
   haplotype and is not used here. Therefore the positive findings establish
   detection of engineered assembly-versus-synthetic-path disagreement, not
   original-read support or an assembly error in a biological sample.
3. The public inputs include engineered unique flanks and three supplied
   source-derived consensus candidates. This legitimately narrows a local
   development task, but does not test native anchor discovery, monomer
   discovery, candidate catalogue completeness, partial spanning reads,
   systematic motif errors, homopolymers, mixed haplotypes, or realistic
   coverage. Native biological monomer boundaries and HOR periods are
   unverified. V3 includes no shared-family, nested/compound, confirmed HOR,
   or cross-donor held-out condition. It cannot meet the user's full B1 or
   Gate B requirements.
4. Protocol and metrics were declared frozen before the archived V3 route
   predictions; the route sources were frozen before V3 truth scoring. This
   is useful development discipline. Truth and source arrays coexist in the
   same accessible repository, so it is an **audit split, not secure
   blinding**. Commit time/order alone cannot prove that no author saw the
   causal truth during route or case design. Any later change after looking
   at these failures would consume V3 as development and require an entirely
   new held-out lineage for a final claim.
5. The one-lineage row totals permit a **development route preference**:
   alignment produces usable calls and six more true detections than the
   frozen simple length baseline, whereas graph rejects everything. They do
   not justify the stronger statement that alignment reliably beats the
   baseline on independent arrays, real errors, or native reads. Neither
   precision nor specificity is established by six engineered negatives of
   one kind with no errors; both observed FP counts are zero on this narrow
   control.

## M1 decision and method selection

The M1 `NON_IDENTIFIABLE` state is justified only for the supplied circular
5-mer signature matrix: exact nullspace involvement prevents separate
coefficients from being inferred from those aggregate features. Full rank
plus an exclusive 5-mer does **not** imply representative completeness,
read-callability, calibrated physical copy count, or robust inference from a
real genome. The reported 395/540 rejected `f2` source reads in the
family-specific error condition are direct evidence that the matrix test alone
cannot solve assay bias. The previous A–F tournament and occupancy failures
support the documented **NO-GO for current absolute-abundance novelty**. The
diagnostic's labels should never be presented without their feature/input
scope; the production estimator and TSV have not been retrofitted with a
validated identifiability gate.

The V3 tournament rejects the exact-transition graph route for the present
native-diversity input: 39/39 abstentions, including all six intact controls.
The alignment route is the only current nontrivial M2 development candidate.
It does **not** pass Gate B: no held-out structural-error set exists and the
read evidence is synthetic. Selecting it for the next independent validation
experiment is a development choice, not a scientific winner declaration.

## Competitor and resource audit

The 36,371-bp assembly-only comparator is one engineered 171-bp-monomer,
40-HOR-copy control with four injected variants. TideHunter/TideCluster's two
stages completed and localized one array with base recall 1.0 and precision
0.9999418 against that exact ledger; a reported 855-bp repeat period matches
the canonical five-by-171 length. Those stages do **not** emit a monomer
label path or variant HOR organization. CENdetectHOR failed technically
before final native calls in the no-prior primary run; subsequent 5-kb
window checks also failed and one supplied truth-derived 171-bp prior is
ineligible for a no-prior ranking. TandemX has no corresponding assembly-only,
no-prior run. Thus there is a common FASTA and one valid array-localization
observation, **no completed common HOR accuracy endpoint or three-way
ranking**. Scoring assembly-only tools as negatives for a read-supported
discordance task would change the task and be invalid.

The B1 V3 receipts report alignment 164.211 s and 23.10 MB peak process RSS,
baseline 0.00876 s and 17.50 MB, graph 0.00196 s and 20.50 MB. These are
different local Python route invocations, excluding startup/output; the graph
time measures early abstention. They are not complete workflow or scaling
benchmarks. Comparator profiling instead sampled *process-tree summed RSS*
at 0.2-s intervals and includes distinct external programs/stages; its
TideCluster clustering peak was about 8.2 GiB on this tiny control, with
shared-page double counting possible. These memory values cannot be ranked
directly against Python `ru_maxrss` or extrapolated to 16-Gb genomes.

## Gate disposition at this checkpoint

| Gate | Adversarial finding |
| --- | --- |
| A: stable new quantitative M1 method | **NO-GO** on current evidence; preserve read-derived evidence/audit-score scope. |
| B: held-out M2 structural discordance | **UNTESTED**. V3 is a favorable one-lineage development set with synthetic duplicate read paths. Alignment is a candidate for a future frozen held-out test only. |
| C: biological utility with orthogonal support | **UNTESTED**. No verified same-locus original reads or independent biological prediction is in this benchmark. |
| Assembly-only HOR comparison | **INCOMPLETE**. TideCluster array localization succeeded; CENdetectHOR final result and TandemX same-prior run are absent. |

The next valid Gate B test must freeze alignment and the simple baseline,
enroll a distinct high-quality source lineage with independently verified
original same-donor/haplotype reads, and retain realistic errors, mixed or
partial spanning reads, and intact negatives. It must report all-case
abstention cost, per-array rather than per-edit uncertainty, and the endpoints
that are actually observable. A failure of that test must remain a failure;
it cannot be rescued by retuning on the held-out truth.

## Addendum: B2 synthetic held-out audit after source/scorer freeze

The separately seeded B2 source, input/truth bundle, metrics and scorer were
committed as `6c4df79` before the reported route execution. The B2 route
predictions and score directories were still untracked at the time of this
audit; their files and hashes are explicit below. This is a useful frozen
development-to-synthetic-held-out transition, but repository commit order
alone is not secure proof that no person could inspect the committed truth.
The evaluation's statistical units are **four synthetic founders** at
100/200/300/400-bp monomer lengths, each subjected to the same 13 case
design, not 52 independent genomes. It contains 44 injected positives, eight
intact controls and 48 unique edited assembly sequences. Synthetic read error
realizations are independent in noisy rows; there are zero biological donors,
native raw reads, or independent biological copy truths.

I independently keyed B2 raw inputs, causal truths and four prediction files
by case ID and applied the frozen 0.5/non-`ok` rules without calling the B2
scorer. All counts below match the archived `summary.json` files.

| Frozen route | Prediction SHA-256 used by scorer | TP / FN / FP / TN / unresolved intact | `ok` / abstain | All-case sensitivity | Negative failure rate |
| --- | --- | ---: | ---: | ---: | ---: |
| Monomer-path alignment | `292dce688e4d5c903cf143fb5e0a49c1797d8de35d247e0e733d1dab765ae428` | 6 / 38 / 0 / 1 / 7 | 7 / 45 | 6/44 = 0.1364 | 7/8 = 0.875 |
| Transition graph | `de1f45f9244e47827f9cc6a061b623f4f286b06b105b701b7ec207f1526a47ee` | 0 / 44 / 0 / 0 / 8 | 0 / 52 | 0 | 8/8 = 1.0 |
| Two-flank length baseline | `70fc14fe760970f80240bde8f1a353312df0d892f47aae304dbaba5ac8d19f04` | 20 / 24 / 0 / 4 / 4 | 28 / 24 | 20/44 = 0.4545 | 4/8 = 0.5 |

The baseline's raw prediction file SHA-256
`d46d2f0160a753a75773badf8b9da3353a1ae077b309757736c6e4ff06ecbc27`
has the non-schema status `insufficient_read_support` on exactly 24 rows.
The scored adapter changes **only** those rows' status to `abstain` and
prefixes their reason with `schema_v2_status_normalization:`. Event scores,
decisions, case IDs and other fields are unchanged. Treating the raw status
as non-`ok` gives the same 20/24/0/4/4 confusion counts, so the adapter does
not improve the baseline outcome. It is a schema compatibility repair that
must be disclosed, not hidden as a different method run.

All six alignment true positives and its sole true negative occur in the
100-bp founder. At 200, 300 and 400 bp, each 13-case founder receives no
`ok` call. Of the 45 alignment abstentions, 38 report the frozen 4,096-bp
array ceiling, six report too few resolved anchored reads, and one reports
candidate-monomer length outside the route's 4..300-bp bound. The 24
one-read B2 cases are all non-`ok` for both alignment and baseline, as their
three-read support rules require. The graph route rejects all 52 because its
frozen adapter admits `split=development`, whereas B2 correctly declares
`split=synthetic_held_out`; its zero calls are an **interface eligibility
failure**, not evidence that a graph architecture was evaluated on the B2
sequences. These declared operating limits are legitimate reasons to abstain,
but their cost belongs in the primary denominator.

Alignment's eligible-only 6/6 = 1.0 sensitivity applies to just six positive
cases and 7/52 total coverage. The baseline's eligible-only 20/24 = 0.8333
applies to different eligible cases and 28/52 coverage. Neither conditional
rate overrides the all-case comparison. All three nominal intact FPR values
are 0/8; seven, eight, and four unresolved intact controls respectively make
those values poor specificity evidence. Event subtype, breakpoint, signed-bp,
HOR-period and order summary metrics are null for incomplete route outputs.
The B2 core-time receipts report 12.217 s and 22.74 MB peak process RSS for
alignment versus 0.0123 s and 18.83 MB for the baseline, excluding process
startup and output writing. The graph's 0.000079 s measures only its split
gate and has no RSS receipt. None of these is a completed workflow or
large-genome scaling measurement.

**Stop-loss disposition:** On this frozen synthetic held-out design, the
current alignment route does not beat the simple two-flank length baseline:
6 versus 20 true injected-edit detections at the same zero explicit false
positives, with 45 versus 24 abstentions. The graph route is unusable under
its frozen input contract. The user-defined M2 superiority stop-loss is
therefore triggered for the **present prototypes**. The result is not a
biological Gate B pass or a general theorem that architecture auditing cannot
work. Because B2 truth is now consumed, any redesigned route or expanded
operating range must be treated as new development and evaluated on another
independent, prespecified held-out source; B2 cannot be reused as final
validation. Biological Gate B and Gate C remain untested.
