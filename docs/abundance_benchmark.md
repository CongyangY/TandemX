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
treated as independent biological replicates. Held-out seeds 5101–5103 remain
unused until method/settings are frozen. The initial design confounds founder
length and copy count; a factorial extension is required for separate effects.

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
