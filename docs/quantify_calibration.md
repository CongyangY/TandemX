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

After committing the source/configuration and passing hosted CI, run:

```bash
conda run --no-capture-output -n tandemx-dev \
  python -m benchmarks.scripts.evaluate_quantify_calibration \
  --config benchmarks/configs/quantify_calibration_development_v1.yaml \
  --datasets \
    /path/to/factorial_scale_s6301_v1 \
    /path/to/factorial_scale_s6302_v1 \
    /path/to/factorial_scale_s6303_v1 \
  --outdir /path/to/new/result
```

Do not select a default from this development ablation alone. Freeze a separate
family/genome-held-out configuration after reviewing failure modes, then run it
once without changing the rule or silently removing failed conditions.
