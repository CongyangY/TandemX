# Public quantify depth-calibration ablation

`benchmarks/configs/quantify_calibration_development_v1.yaml` freezes a paired
development experiment on the three existing seed-6301--6303 factorial genomes.
Each of 27 read conditions is passed to the public `tandemx quantify` command
under four methods:

1. total read bases divided by genome size;
2. the same depth with the planted aggregate error rate supplied to the
   independent-error k-mer survival approximation;
3. 1,000 whole-genome-verified, non-array, non-repeat controls; and
4. controls plus the planted aggregate error rate.

The planted error rate is an oracle input used only as an upper-bound ablation.
It is not available in a blind real analysis. FASTQ Phred values provide an
observable but calibration-dependent substitute. Controls are chosen with
simulation truth and verified to occur exactly once after reverse-complement
canonicalization; a real workflow needs independent repeat masking, population
copy-number checks and contamination screening.

The evaluator retains every per-family estimate, truth value, sampling-oracle
difference, diagnostic-spread endpoint, normalization audit field, command,
failure, wall time, CPU time and direct-child peak RSS. Families within a genome
and read conditions generated from that genome are dependent. Resource rows
have one execution per method-condition. The 10th--90th percentile endpoints
remain diagnostic-k-mer spread and are not relabelled as sampling confidence
intervals. The separate joint-read multi-k experiment evaluates approximate
sampling intervals and missingness.

The first complete development run used committed source `c8b3a8b` after both
hosted workflows passed. It completed 108/108 commands, 5,940 family rows and 36
strata. Aggregate mean absolute relative error was 0.401898 for total-bases
normalization, 0.359640 for the planted-error oracle, and 0.365728 for empirical
controls. Controls improved all three independent-genome means but were worse at
nominal 1x. Controls plus oracle generated exactly the same estimates as controls
alone, so the global survival factor is redundant when it corrects both numerator
and denominator. The exported 10th--90th diagnostic spread had very low truth
inclusion and is not a sampling confidence interval.

The original explicit-error implementation also added substantial FASTA scan
time. After the exact-formula all-ACGT shortcut passed hosted CI, the full matrix
was replayed into a new T7 directory. All 5,940 metric rows and 108 copy-number
files were byte-identical; all three control panels and every non-resource summary
field also agreed. Total driver time changed from 610.507 to 349.217 seconds
(-42.80%). Median oracle-error time changed from 4.658 to 1.944 seconds (-58.26%)
and controls-plus-oracle changed from 4.934 to 2.178 seconds (-55.85%). Peak
driver RSS changed from 108.234 to 106.844 MiB. These are one historical run and
one same-machine replay, so they are an engineering regression result rather than
a timing distribution. Compact comparison evidence is in
`paper/evidence/quantify_calibration_fast_fasta_replay_v1`.

The development-selected candidate uses empirical controls only when their mean
depth is at least 2. Its threshold was carried unchanged into the separately
frozen validation below.

## Frozen untouched-genome validation

`factorial_scale_quantify_validation_v1.json` retains the development data-
generating process but replaces its seeds with untouched validation genomes
6401--6403. `quantify_depth_gated_validation_v1.yaml` freezes the candidate and
nine gates before data generation. It also freezes SHA-256 values for the data-
generation JSON, complete Mo17 length histogram and all five generator/helper
source files; every dataset receipt must match them. The evaluator runs two observable public
command modes: total-bases normalization and the same 1,000 simulation-truth-
assisted controls. It selects controls for an entire read condition only when
their mean depth is at least 2; it does not refit the threshold.

Promotion within this conditional known-catalogue IID simulation requires all of
the following: zero failed commands; aggregate MARE reduction >=0.02; candidate
MARE <=0.40; positive MARE reduction for every independent genome; no mean
coverage-stratum regression beyond numerical tolerance; paired improved and
nonworse fractions >=0.40 and >=0.70; and at least one condition routed through
each branch. Individual families may regress and remain in the paired table.
Resource observations are excluded from these scientific gates.

Commit `9eda196` and hosted runs `34088623796`/`34088640260` passed before all
three datasets were generated with `--split validation` as documented in
[factorial_scale_simulation.md](factorial_scale_simulation.md). An independent
audit checked 30 manifests and every one of 93 declared payload files. The
frozen evaluation was then run once with:

```bash
conda run --no-capture-output -n tandemx-dev \
  python -m benchmarks.scripts.evaluate_quantify_depth_gated_validation \
  --config benchmarks/configs/quantify_depth_gated_validation_v1.yaml \
  --datasets \
    /Volumes/T7/Codex/TandemX/data/simulated/factorial_scale_s6401_validation_v1 \
    /Volumes/T7/Codex/TandemX/data/simulated/factorial_scale_s6402_validation_v1 \
    /Volumes/T7/Codex/TandemX/data/simulated/factorial_scale_s6403_validation_v1 \
  --outdir /Volumes/T7/Codex/TandemX/results/quantify_depth_gated_validation_v1_20260907
```

Seeds 6401--6403 are consumed after this once-only evaluation. Seeds 7401--7403
remain reserved and the generator refuses them.

All 54 public commands completed without timeout, producing 2,970 raw and 1,485
candidate family rows. The rule passed all nine predeclared gates. Aggregate MARE
changed from 0.408767 for total-bases normalization to 0.363222 for the candidate,
an absolute reduction of 0.045545. Per-genome reductions were 0.037855, 0.043531
and 0.055248; the 1×, 5× and 20× reductions were 0, 0.050816 and 0.085818.
Among 1,485 paired family conditions, 661 improved, 495 were equal and 329
worsened. Nine conditions used the total-bases branch and 18 used controls.

Ungated controls reached a slightly lower aggregate MARE of 0.359443 in this
validation split. The gate therefore validates robust improvement over the
predeclared total-bases baseline, while remaining conservative relative to the
controls-only result. This observation is retained and was not used to retune the
threshold. Median baseline/controls wall times at 1×, 5× and 20× were
0.500/0.518, 1.616/1.836 and 5.887/6.605 seconds. These are single executions,
not repeated publication timing.

The validated rule is exposed as the opt-in public option
`tandemx quantify --single-copy-kmers controls.tsv --single-copy-min-depth 2`.
Below threshold it records
`total_read_bases_divided_by_genome_size_low_control_depth_fallback`; the default
controls behavior is unchanged pending real-control specificity validation.
The 67-entry compact archive and accepted editable six-panel figure are in
`paper/evidence/quantify_depth_gated_validation_v1/figures_v3`. Figure v1 is
retained as a rejected legend-overlap layout; v2 fixed that layout but lacked
unique runtime/RSS pairing keys in its panel-source table.

The executed command was:

```bash
conda run --no-capture-output -n tandemx-dev \
  python -m benchmarks.scripts.evaluate_quantify_calibration \
  --config benchmarks/configs/quantify_calibration_development_v1.yaml \
  --datasets \
    /path/to/factorial_scale_s6301_v1 \
    /path/to/factorial_scale_s6302_v1 \
    /path/to/factorial_scale_s6303_v1 \
  --outdir /Volumes/T7/Codex/TandemX/results/quantify_calibration_development_v1_20260907
```

Compact tables, per-execution hashes, the recomputed decision and the accepted
six-panel figure are in `paper/evidence/quantify_calibration_development_v1`.
Do not select a default from this development ablation alone. Freeze a separate
family/genome-held-out configuration after reviewing failure modes, then run it
once without changing the rule or silently removing failed conditions.
