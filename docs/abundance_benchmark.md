# Conditional copy-number and assembly comparison experiment

This development experiment supplies a **known catalogue** to isolate the
quantification and localization stages. It does not measure de novo catalogue
recovery. The generator/evaluator import no TandemX algorithm code. Source
snapshots, commands, failures, hashes and raw metrics accompany every execution.

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.abundance.run \
  --config benchmarks/configs/abundance_v1.json \
  --split development --outdir /path/to/new/abundance_experiment
```

The three development genomes use seeds 4101–4103, separate from discovery
development. Each contains independent 61/171/421-bp founders with 20/80/200
copies, separated by independently sampled 25-kb flanks (199,100 bp total).
Exact-copy arrays are embedded in a bounded genome retained in memory; the
generator rejects sizes above two million bases. Read output is streamed.
Five assembly versions retain 100%, 75%, 50%, 25%, or 0% of each planted array;
monomers and background remain identical, and truth coordinates are regenerated.

Read starts are uniform over the circular source genome, on either strand.
The circular sampling model avoids chromosome-end coverage depletion. Read
length is 5,000 bp, with target coverages 1/5/20× and independent substitution
rates 0/0.1/1%. Integer read counts make actual depth slightly different from
requested depth. Starts/strands are held constant across error settings and
recorded separately from neutral FASTA IDs. Errors are applied to all bases.
This is not an empirical HiFi error model: no indels, diploidy, biological unit
divergence, coverage bias or unknown background homology is represented yet.

Each genome is one independent simulation realization. Coverage/error/assembly
conditions reuse that genome; its three families and technical runs must not be
treated as independent biological replicates. Held-out seeds 5101–5103 were
reserved until the method/settings were frozen and were used once in the
held-out execution described below. They are now consumed and must not be used
for tuning. The initial design confounds founder length and copy count; a
factorial extension is required for separate effects.

## Endpoints and controls

- Copy number: signed/absolute relative error, interval width and empirical
  truth coverage. Current output intervals are diagnostic-k-mer 10–90% spreads,
  **not sampling confidence intervals**. Their observed coverage is measured
  without assigning them a nominal 80% or 95% coverage claim.
- Sampling oracle: truth repeat bases actually sampled, divided by founder
  length and actual base coverage. Its departure from genomic truth measures
  random sampling fluctuation; estimator-minus-oracle measures additional bias.
  The oracle is evaluation-only and never enters the tool command.
- Localization: per-family interval-union recall/precision, represented bases
  and fragment count. Missing denominators are NA. Empty native BED is valid
  when the stage succeeds; failed commands remain failed.
- Under-representation: planted assembly/genome ratio below 0.6 is positive.
  `possible_collapse` and `reads_only` are counted as detection for this explicit
  endpoint, while native status is retained. Ratios 0.75 and 1 are controls.
  This tests recognition of substantial under-representation, not every loss.
  Report TP/FN/FP/TN, sensitivity, false-positive rate and precision separately.
- Resource records: wall, user/system CPU and direct-child RSS per command.
  Current runs diagnose method behavior; they are not cross-tool rankings.

`genomes/` holds FASTA, known catalogue, full/retained truth and file hashes;
`reads/` holds FASTA, sampling coordinates and observed error counts. `runs/`
retains every CLI receipt/log/output. Root metric TSVs are one row per family
and condition; summaries aggregate these rows explicitly. `validation.json`
reports successful execution separately from scientific acceptance. Raw source,
reads and large outputs stay outside Git; compact tables may accompany the paper.

Tests: `pytest -q tests/unit/test_abundance_benchmark.py` checks independent
base-occupancy truth, reproducibility, reverse strands, circular crossings,
unchanged flanks, score denominators, invalid values, duplicate predictions and
a real five-command CLI workflow including complete assembly absence.

## First executed baseline

`results/abundance_baseline_v1_20260906` at the T7 data root used committed
8218beb source: 177 successful commands, 81 CN observations, 45 assembly-family
and 405 comparison-family observations. Compact evidence is in
`paper/evidence/abundance_baseline`. Only 6/81 k-mer spread intervals contain
genomic truth. At 20×/1% substitution, mean signed error is -18.31%; the sampling
oracle has +1.77% mean error and estimator-minus-oracle is -20.08%. Error-free
estimator-minus-oracle is much smaller (about -0.36% to -0.52%), while 1× random
sampling itself produces large deviations. Both false under-representation calls
and misses occur. Simple exact-copy localization has base recall 1 for every
positive family. These results motivate uncertainty/error-model development;
neither three genomes nor perfect exact-copy localization establish robustness.

The next [streamed factorial experiment](factorial_scale_simulation.md) removes
the fixed length/copy pairing, raises each source genome to 10 Mb, adds variable
GC and biological divergence, a megabase array, empirical read lengths and
indel-aware observed-coordinate truth. This is an additional controlled model;
the original baseline and its uncertainty failures remain evidence.

## First held-out execution

The unchanged configuration was executed once on seeds 5101–5103 using source
`f16596b`; all 177 commands completed. Across 243 positive under-representation
and 162 control family conditions, TP/FN/FP/TN were 208/35/8/154, giving
sensitivity 0.855967, false-positive rate 0.049383 and precision 0.962963.
All 36 positive assembly/family localization rows had base recall 1.0, with
minimum base precision 0.996732. These exact-copy arrays make localization an
easy conditional control rather than representative biological validation.

At 20× with 1% substitutions, sensitivity for 50%-retained arrays fell to 2/9.
At 1×, the fully retained control produced one false call in each error tier.
The aggregate therefore does not justify a robust collapse claim. The compact
archive `paper/evidence/abundance_heldout` contains every metric, verified source
and configuration hashes, a manifest and resource rows reconstructed from all
177 validated receipts. A revised model requires development-only work and a
newly predeclared independent seed set.

## Development-only multi-k decision calibration

The completed development matrix for seeds 4101–4103 was reused without
rerunning its original single-k quantification or localization commands.
`evaluate_multik_collapse.py` recomputed only k=15/21/27/31 read evidence and
compared the original k=21 decision, a multi-k log-linear estimate and a
transparent hybrid. The frozen rule uses the multi-k estimate when available,
falls back to k=21 otherwise, and lowers the ratio decision threshold from 0.6
to 0.5 only when observed haploid depth is below 2.0. All three leave-one-genome-
out development folds selected 0.5 independently.

On the development rows, the original decision had TP/FN/FP/TN=195/48/7/155;
the hybrid had 203/40/3/159. Sensitivity increased from 0.802469 to 0.835391,
false-positive rate decreased from 0.043210 to 0.018519 and precision increased
from 0.965347 to 0.985437. Multi-k alone had 35 unavailable rows, so its
conditional metrics are not a complete-method comparison. This is calibration,
not validation. `benchmarks/configs/abundance_v2.json` froze the rule, the
calibration hashes and untouched held-out seeds 5201–5203 before their one-time
execution from committed source 1231743.

## Frozen-rule validation on fresh seeds

Both hosted CI runs for 1231743 passed before the new matrix was run. The first
and only 5201–5203 baseline completed 177/177 commands. The original k=21 rule
had TP/FN/FP/TN=198/45/14/148: sensitivity 0.814815, false-positive rate
0.086420 and precision 0.933962. Applying the frozen rule without fitting on
these rows changed the counts to 208/35/12/150: sensitivity 0.855967,
false-positive rate 0.074074 and precision 0.945455. Each seed gained true
positives without gaining false positives, although seed 5202 retained 12 false
positives and a 0.222222 false-positive rate.

At 20×/1% substitutions/50% retention, sensitivity improved from 3/9 to 9/9.
At 1× in complete-assembly controls, false calls decreased from five to three
of 27 family conditions. The gain was not uniform: at 1× with 0% or 0.1%
substitutions and 50% retention, sensitivity decreased from 6/9 to 4/9 in each
stratum. The multi-k-only arm had 15 unavailable rows. These
results validate this narrow decision rule on fresh exact-copy simulations;
they do not establish divergent-array, unknown-catalogue, empirical HiFi or
biological under-representation performance. All 5201–5203 seeds are consumed.
Compact paired archives are `paper/evidence/abundance_heldout_v2` and
`paper/evidence/abundance_multik_collapse_heldout`.

## Predeclared divergence and fragmentation extension

`benchmarks/configs/abundance_domain_shift_v1.json` defines a second frozen-rule
validation with unused seeds 5401–5403. It crosses 1%, 3% and 5% independent
founder-to-unit substitutions with one or three same-family array segments.
Three-segment arrays contain 500-bp independently generated interruptions, which
are longer than the localization merger allowance. Each assembly fraction keeps
the same deterministic prefix of unit variants and regenerates exact interval
truth. The original coverage, read-error and five retention levels are retained.

This experiment evaluates domain shift of the already frozen multi-k/depth rule;
it does not recalibrate on the new factors. Substitutions are independent and
length preserving, and fragmentation represents interrupted same-chromosome
arrays rather than contig breaks, graph assemblies or empirical satellite
evolution. The configuration and implementation must be committed and pass CI
before the 5401–5403 baseline is executed. Those seeds are then single use.
