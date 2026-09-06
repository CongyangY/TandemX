# Fixed multi-k conditional replay

Source 0aad0ab82c50f93f0c80351b3b4b3c1749c1908e. All 27 observed-read inputs
from the three factorial development genomes were scored: 1,485 family conditions,
4,455 paired method rows, fixed k15/21/27/31. The estimator saw only observed
reads, known founder sequences and fixed genome size. No truth/error parameter
was an estimator input. All source/input/baseline hashes and commands are retained.

At 20x and 2% unit divergence with high read errors, the native median k21 mean
signed error is -56.586%, versus +1.590% for the prototype; mean absolute relative
error is 56.586% versus 13.529%. Mean-k21 alone gives -56.552% bias and 56.552%
absolute error. Thus the main change in this controlled stratum comes from the
extrapolation, rather than only changing median to mean or finite read exposure.
Each stratum has 81 families across three independent source genomes, not 81 plants.

Limitations are preserved. At 20x, 384/495 paired absolute errors improve and
111 increase; at 1x, 243 improve, 189 increase and 63 cannot be fitted. At 5x,
324 improve, 168 increase and three cannot be fitted. There are no exact ties
at the declared 1e-12 tolerance. Overall 66 fits have zero support at some k;
two additional positive slopes are retained with a model-violation flag. These
two are included as available point estimates, not declared reliable estimates.
No sampling CI is supplied. Background specificity, real HiFi biases, ploidy,
related families, held-out generalization and calibrated uncertainty remain open.

The six-panel `figures/factorial_quantification.pdf` and editable SVG show both
gains and failures. Source rows, plot script hash and all input/output hashes are
retained; see `figure_legend.md`. Resources are diagnostic during other jobs:
192.175 s for the whole replay, direct-child peak 162.36 MiB. This is not a
matched external-tool timing comparison or a mature software release.

Reproduce the plot:

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.plot_factorial_quantification \
  --replay paper/evidence/factorial_multik_replay \
  --baseline paper/evidence/factorial_quantification_baseline \
  --outdir /path/to/new/figure_directory
```
