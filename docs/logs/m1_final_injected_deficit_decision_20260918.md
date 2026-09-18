# Final M1 injected-deficit decision (2026-09-18)

## Scope and frozen sequence

The quantitative missing-bp claim for the current production M1 estimator is
**NO-GO**. Its permitted research position is **Level C scope: binary
discordance/priority evidence only**, with no calibrated bp magnitude and no
validated ordinal severity score. This is a scope decision, **not** an accuracy
claim for a Level C binary detector. That detector still needs independently
frozen structural and biological validation. The public estimator, historical
results, M2 rules and manuscripts were not changed.

Protocol `benchmarks/m1_final_deficit/protocol.json` was committed at `8386133`
before running any panel (SHA-256
`d415717c946a0b5235bc4d7d35c5b4113f5c4ca9f2d5553c3910c06c3586e246`).
The source generator used seeds 73101–73103, three 80-bp families with 25,
40 and 60 planted copies, 0%/12% unit divergence, uniform exact 7-kb reads at
nominal 5/10/20/30×, and fixed reads across five assembly copy fractions
1/.75/.5/.25/0. The source genome is 22,000 bp. Twenty-four independent
**simulation panels** were run; 72 family-panel observations and five related
assembly edits per observation yield 360 method rows. No panel is an
independent biological donor. The edited assemblies were regenerated with
identical catalogue/background and hashed; the exact source and retained family
bp are retained per row. `results.json` records every panel's read hash,
method estimate, copy-edit truth, signed prediction and error.

The first actual run is preserved in
`benchmarks/m1_final_deficit/evidence_v1_20260918/`. Its historical 80-bp
`ordinary_segments` mapping omits reverse-complement search, a concrete
fairness defect discovered after the first score. It must **not** be ranked as
the competitive baseline. Before correction, `protocol_v2.json` was committed
at `7cfb915` (SHA-256
`c0e585559208e5b188fdc4f5e5e13d5927f75bd3b4d9c5dfc79dd9143fc1f5eb`),
with the v1 result hash frozen. The corrected 80-bp periodic-edit baseline
tests both strands and abstains on ties or failed alignment. It regenerated
every read and checked its SHA-256 against v1 before scoring. The corrected
result and combined metrics are in
`benchmarks/m1_final_deficit/evidence_v2_20260918/`. No result was overwritten.

## Endpoint and outcome

For family *f* and assembly edit *e*, injected truth is
`source_bp[f] - edited_assembly_bp[f]`; predicted signed deficit is
`read_estimate_bp[f] - edited_assembly_bp[f]`. Their difference is
`read_estimate_bp[f] - source_bp[f]`, constant across all five edits for a
fixed read panel. Therefore **monotonicity has slope one by algebra**, even for
an inaccurate read estimator. Scoring it as proof of quantitative accuracy or
counting the five edits as independent replicates would be invalid. The primary
calibration unit is the intact-source family-panel error; assembly edits retain
an exact provenance test and demonstrate deficit reporting, not new
independent read evidence.

The prespecified Level A gate required every family at both 20× and 30× in
all three seeds and both divergence settings to have intact error ≤10% of
source bp and edit-panel error ≤20%. No silent NA is permitted. No family-wide
method passed. Level B needs **independent changes in read-side source
abundance**; the fixed-read edit panel cannot establish it. Existing M1
tournament and identifiability evidence also supply no calibrated ordinal
gate. The final position for the current production M1 claim is consequently
Level C scope, and its quantitative missing-bp claim stops now.

| Method | 360 edited-assembly rows with a numeric result | Mean absolute error, bp¹ | Worst error/source bp | Level A gate |
| --- | ---: | ---: | ---: | --- |
| Current production diagnostic-k-mer estimator | 360 | 1747.6 | 97.16% | fail |
| Corrected dual-strand competitive mapping | 360 | 371.2 | 37.50% | fail |
| Frozen full-read research prototype | 360 | 370.7 | 38.97% | fail |
| Exact-flank median read span, one array/family | 310 | 0 conditional on decision | 0 conditional | pass at 20/30× only² |
| Historical forward-only mapping (invalid fair baseline) | 360 | 1878.7 | 92.14% | fail |

¹ The five related edits repeat the same read-estimation error. This MAE is a
descriptive average over 72 family-panels replicated five times, **not** 360
independent accuracy measurements. ² The span arm had decision coverage
10/18, 16/18, 18/18 and 18/18 source-family observations at 5×, 10×, 20× and
30×, respectively (62/72 observations, each repeated over five edits).
Its exact 32-bp flank matching, error-free simulated reads and one planted
array per family make the observed 0-bp conditional error unsurprising. It
does not validate real noisy-read structural accuracy, genome-wide family
abundance, or long unspanned arrays. The three-read minimum is source-local
and is explicitly not a whole-genome coverage floor.

At 20×/30× with 12% unit divergence, the production estimator's mean intact
family error was approximately 3120/3118 bp across the nine family-seed
observations per depth, compared with 264.5/269.4 bp for full-read and
261.6/269.4 bp for corrected competitive mapping. The production method's
diagnostic founder k-mers lose nearly all exact observations under this
within-family divergence; no parameters were adjusted. Even on 0% divergence,
only 12/18 production intact-family observations at 20×/30× meet the 10%
intact-error gate. The full-read prototype does not improve decisively over
corrected competitive mapping on this panel. The prototype and baseline use
the same observed read depth and retain unassigned mass; no ambiguous bp are
redistributed.

On the 72 intact family-panels, mean signed bias was −1543.62 bp for production,
+2.65 bp for corrected competitive mapping and −7.41 bp for full-read; mean
absolute error relative to each source family's bp was 52.76%, 12.16% and
12.19%, respectively. Opposite signed errors can cancel in the latter two
biases. Their 20×/30× absolute errors remained 1696.43, 267.30 and 271.23 bp;
more coverage did not fix production's 12% divergence failure. Across the
three source seeds at 30×, the largest range in an individual family estimate
under otherwise fixed conditions was 486.8 bp (production), 494.8 bp
(corrected mapping) and 528.3 bp (full-read). These seed ranges quantify
synthetic sampling variation, not biological reproducibility. Full family-level
scores and all coverage strata remain inspectable in the archived JSON files.

## Boundaries and disposition

This experiment has **injected assembly deficit truth**, not physical
copy-number truth. Source bp come from the simulator; reads are uniform,
error-free and circular, with a known three-family catalogue and known source
genome length. At a fixed read depth, zeroing the assembly eventually yields a
positive predicted deficit for almost any nonzero read estimate; that is not
binary sensitivity evidence. The native Ey15-2/Macadamia read-span edits are
separate correlated development results with uncertain extraction pairing.
No real-sample Level A, B or C operating characteristics can be inferred.
The future structural-discordance task remains SUPPORTED/DISCORDANT/UNRESOLVED
with explicit reasons and decision coverage; full HOR reconstruction remains
NO-GO. The next scientific work is independent sample-cluster enrollment,
controlled structural edits and prospective orthogonal corroboration. No new
formal runtime claim was made on USB2 storage.

Reproduce v1 with `conda run -n tandemx-dev python -m
benchmarks.m1_final_deficit.run --outdir NEW_EMPTY_DIR`, then v2 with
`conda run -n tandemx-dev python -m benchmarks.m1_final_deficit.run_v2
--v1-results NEW_EMPTY_DIR/results.json --outdir ANOTHER_EMPTY_DIR`.
The v2 script checks the exact frozen v1 result bytes, so a reproduction must
also match the archived v1 SHA-256 before it accepts the input.
The archived result integrity audit confirmed all 24 paired read hashes, all
1800 per-method case rows, exact edit arithmetic and the constant-error
identity. Full `tandemx-dev` pytest passed 1089 tests in 142.24 s. A subsequent
fail-closed duplicate-read-ID guard (which does not affect unique-ID archived
inputs) passed its focused two-test check; the full suite was not rerun after
that guard-only change. No formal runtime or peak-memory comparison was made.
