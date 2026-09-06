# Streamed factorial genome simulation

This independent generator extends the small conditional-abundance experiment.
The `factorial_scale_v1.json` design crosses three founder lengths (61/171/421 bp),
three copy counts (20/80/200), three target GC fractions (0.3/0.5/0.7), and two
within-array substitution probabilities (0/0.02). The 54 independent founders
occupy balanced factorial cells. A separate 171-bp founder has 6,000 copies
(1,026,000 bp) and 1% unit divergence. Each source genome is 10 Mb with IID
45%-GC target background. Actual founder GC and observed biological substitutions
are recorded. These are controlled haploid simulations, not plant genomes.

The generator writes at most a 65,536-bp background chunk or one <=20-kb monomer
at a time, retaining less than one FASTA line between writes. A fixed-width
single-contig index supports bounded random access. It does not load the whole
genome or all reads. Explicit limits: 1 Gb source genome, 200-kb read length,
100x source coverage and the user-supplied worst-case output-base budget.
The budget covers sequence bases, not filesystem overhead/metadata.

Read lengths can be drawn from the complete Mo17 QC histogram with its empirical
frequency distribution. Starts are uniform on a circular genome; strands are
random. Error/coverage settings share a design-seed stream and are paired/nested
technical conditions, not independent genomes. Errors use a separate stream on
both background and repeats. Models: error-free; 0.1% substitution + 0.1%
insertion + 0.1% deletion; and 1% substitution + 0.5% insertion + 0.5% deletion.
Each source base can be deleted or retained, a retained base can be substituted,
and a single insertion can follow either outcome. Realized events are counted.
Only read lengths are empirical: qualities, homopolymer context, initiation/GC
bias, ploidy, TE insertions and biological correlations are not modeled.

Sampling stops after the first whole read reaching requested source coverage.
Actual source and observed-base coverages are separate because indels alter read
length. Truth segments are projected into observed read coordinates and reversed
with the read. Keep partial fragments and `at_least_two_source_units`; they must
not silently enter a detectable-array recall denominator. Later scoring must
declare fragment eligibility and treatment of partially labelled regions. The
generator imports no detector/quantifier code. Founders/truth are never de novo
tool inputs. A generation receipt is not a completed software benchmark.

Run one declared development genome (6301/6302/6303); this controller refuses
reserved seeds 7301/7302/7303. Nine conditions use 1/5/20x source coverage.
One genome and its 55 families are not 55 biological replicates.

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.generate_factorial_scale \
  --config benchmarks/configs/factorial_scale_v1.json \
  --length-histogram /path/to/complete_Mo17_qc/length_histogram.tsv \
  --seed 6301 --max-output-bases 1500000000 --outdir /path/to/new/factorial_s6301
pytest -q tests/unit/test_streamed_factorial_simulation.py
```

Tests create an actual <1-MB example, checking every small random-access
boundary, planted copies/biological substitutions, unchanged background,
independent per-base read-origin masks, strands/circular crossings, event
accounting, reproducibility, malformed histograms and held-out/budget gates.

Outputs: `genome/genome.fa`, `genome_index.json`, `catalogue.fa`,
`truth_copy_number.tsv`, `manifest.json`; nine `reads/condition_*/` directories
with `reads.fa`, `sampling.tsv`, `truth_read_segments.tsv`, `manifest.json`;
root source snapshot, config/histogram copies, `generation.log` and a completion
receipt with per-condition hashes. Later scores must retain condition labels
and these hashes, including when a method fails.

## Conditional copy-number scoring

After generation, run the actual public quantifier with the known founder
catalogue and fixed source genome size. This isolates abundance from discovery:
it is not a de novo family-recall experiment. Only observed reads, catalogue,
genome size and k=21 enter the estimator; planted copy counts, source occupancy
and error probabilities are evaluator-only inputs.

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.abundance.run_stream_quantify \
  --datasets /path/to/factorial_s6301 /path/to/factorial_s6302 /path/to/factorial_s6303 \
  --outdir /path/to/new/factorial_quantification --timeout 900
pytest -q tests/integration/test_stream_quantify.py
```

The controller verifies generation/condition receipts, catalogue/truth and read
hashes; it refuses incomplete, repeated-seed or held-out datasets. It snapshots
software, retains commands, stdout/stderr and direct-child resource measurements,
and requires every execution to finish successfully. Failures stay missing and
are not assigned zero errors. `copy_number_metrics.tsv` retains every family,
seed, requested coverage, error model, biological divergence, period and GC.
`copy_number_summary.tsv` groups by coverage, error, divergence and megabase-array
scope. These are descriptive means across correlated family conditions; they
are not independent biological-replicate confidence intervals.

The evaluator separates signed/absolute error against planted copies from
estimator-minus-source-sampling-oracle error. The oracle divides actually sampled
source repeat bases by founder length and actual source coverage. The estimator
sees observed read lengths; indel-altered exposure is recorded separately.
Native 10th–90th percentile diagnostic-k-mer bounds remain a within-k-mer spread,
not a calibrated sampling confidence interval. Resource measurements made during
data acquisition and other jobs are development diagnostics, not final rankings.
