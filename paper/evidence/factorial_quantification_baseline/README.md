# Factorial conditional copy-number baseline

Source 6330d5625c102787c7a894822d1b3ba9599cc92d, development genomes
6301/6302/6303, 10 Mb each, 55 known founder families, nine conditions each.
All 27 public `quantify` commands succeeded and all 1,485 family-condition
rows are retained. Input sequence totals sum to 2,340,179,688 observed bases.
Generation receipts and founder/coordinate truth are archived in `../factorial_inputs`.

The native diagnostic-k-mer spread contains planted truth in only 77/1,485
rows. This spread is not a sampling confidence interval. At 20x, each ordinary
factorial divergence/error stratum has 81 family observations across three
genomes. Mean signed errors for unchanged units are +1.246% (error-free),
-4.880% (low errors), and -32.790% (high errors). With 2% biological unit
substitutions they are -34.791%, -38.656%, and -56.586%, respectively.
The corresponding source-oracle differences are in every source row; observed
sampling fluctuation must not be attributed entirely to the estimator.

The separate megabase-array rows stay separate in the summary. No conditions
are omitted to favor the method. All data use IID background and controlled
errors; this is not empirical HiFi chemistry or real-species validation.
Resource receipts are diagnostics during other jobs, not final rankings.
`archive_manifest.json` verifies exact copied input/output receipts and data.
Raw read FASTA and source snapshots remain at the paths in those receipts.
