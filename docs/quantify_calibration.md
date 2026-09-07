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
depth is at least 2; it requires untouched genome/family validation.

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
