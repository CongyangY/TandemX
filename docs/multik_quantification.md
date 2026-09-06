# Experimental multi-k attenuation model

`tandemx/quantify/multik.py` is a separately testable prototype. It does not
change the public `quantify` estimator or provide calibrated confidence bounds.
The fixed first experiment uses k = 15, 21, 27, 31, declared before scoring the
three factorial development genomes. Held-out genomes are not eligible for this
replay. No sequencing error rate or planted count is an estimator input.

For each k, retain family-exclusive circular canonical words, filter low
complexity, and correct each word's count by its multiplicity in the founder.
If D_k words remain, G is the supplied genome size and reads have lengths L_i:

```
X_k = sum_i max(0, L_i - k + 1)
theta_k = G / X_k * mean_j(count_j / founder_multiplicity_j)
log(theta_k) = a + b * k
extrapolated_copy_number = exp(a)
effective_word_loss_probability = 1 - exp(b), only when b <= 0
```

Fit ordinary least squares with equal weights to log point estimates. The
intercept extrapolates exponential word survival to k=0; it is not a k=0
sequence measurement. Biological divergence and sequencing error both disrupt
matching and cannot be separated by this fit. Under homogeneous IID substitution
loss, exponential survival is a useful model; indels, heterogeneous units, chance
matches, short arrays and nonuniform selection violate its simple interpretation.
No new mathematical principle or AI novelty is claimed.

Keep maximum absolute log residual and the range of leave-one-k-out fitted log
intercepts. The four k values measure the same reads; regression residuals are
**not** independent sampling errors and do not provide a sampling confidence
interval. Positive slopes are preserved and flagged, never clamped into apparent
improvement. Missing targets, zero counts at any k, missing exposure and numerical
failures yield missing extrapolated estimates; there is no pseudo-count or silent
deletion of unsuccessful families.

Reads are parsed once into batches of at most 512 records or a target 8 Mb
(one unusually long record can exceed the target). Existing exact Rust counters
scan each batch for every k; the Python reference is checked against an independent
naive counter and the native implementation. State scales with catalogue targets,
not total reads. N-containing opportunities remain in exposure and are flagged.
This does not measure background specificity or nuclear depth.

## Paired replay

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.evaluate_multik_factorial \
  --baseline /path/to/completed/factorial_quantification \
  --outdir /path/to/new/multik_replay
pytest -q tests/unit/test_multik.py tests/integration/test_stream_quantify.py
```

Compare the retained native median k=21 estimate, an exposure-corrected mean
k=21, and the extrapolation on identical inputs. This ablation separates the
mean/median and exposure change from fitting attenuation. All per-k targets,
exposures/counts, input hashes and family-condition scores are retained. Missing
estimates stay in availability denominators. Supplied founder catalogues isolate
abundance from discovery. Resource values during other jobs are diagnostic.

Required next tests: composition-matched real backgrounds and absent-family
decoys; related/divergent families, higher-order units and interruptions; independent
error profiles, read selection, genome-size/contamination sensitivity; read-cluster
sampling covariance; held-out material calibration. Reduced IID-development bias
alone is insufficient for public CLI integration.

Prior art already models k-dependent survival/coverage: [KmerStream](https://academic.oup.com/bioinformatics/article/30/24/3541/2422237),
[RESPECT](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009449)
and [Merqury](https://doi.org/10.1186/s13059-020-02134-9).
[Satellite reference models](https://genome.cshlp.org/highwire_display/entity_view/node/1041866/full)
also used several k values, background specificity and GC-matched depth controls.
Any novelty claim must distinguish the validated method from established concepts.
