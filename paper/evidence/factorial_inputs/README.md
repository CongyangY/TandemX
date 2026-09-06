# Completed streamed factorial inputs

Three independent development source genomes, seeds 6301/6302/6303, each
10,000,000 bp with 55 planted families. Nine paired technical read conditions
per genome cross 1/5/20x source coverage and three independent-error models.
Observed read bases: 780,054,525 / 780,064,624 / 780,060,539. The total is
2,340,179,688 bases, not 27 independent genomes or 165 independent plants.

Each subdirectory retains exact config, empirical length histogram, founder
catalogue, genomic truth/index, generation and read-condition manifests and
hash-checked archive list. Large sequence/origin/segment files remain under
T7 `data/simulated/factorial_scale_sSEED_v1`; regenerate with the documented
streamed controller. Only read lengths are empirical, from the complete Mo17
QC histogram. Held-out seeds 7301/7302/7303 remain unused.

Methods and limitations: `docs/factorial_scale_simulation.md`. A complete
generation receipt alone does not constitute a software validation result.
