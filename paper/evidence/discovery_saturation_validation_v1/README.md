# Discovery-saturation validation v1

This compact archive supports the coverage-saturation paragraph in manuscript
v5 and Supplementary Figure S1/Tables S3-S4. Full generated FASTA files and the
21 native discovery directories remain under
`/Volumes/T7/Codex/TandemX/{data/simulated,results}`.

## Fixed experiment

- validation seeds: 6501, 6502 and 6503;
- nominal depths: 0.5x, 1x, 2x, 5x, 10x, 20x and 30x;
- nested read streams within each seed;
- 55 planted families: 18 low (20-copy), 18 medium (80-copy) and 19 high
  (at least 200-copy);
- one-to-one cyclic global edit matching at identity at least 0.90;
- transition pass: new-family fraction below 0.05 and family Jaccard above
  0.95;
- saturation: two consecutive transitions pass in all three seeds.

The protocol and configuration were committed as `67ca6ba` before data
generation. Commit `fecbb3e` then enforced sequential depths within each seed
and parallel execution across seeds before any TandemX output was examined.
The source snapshot for the 21 runs records Git head
`fecbb3e49c63fd1bf7323c1bf9e82d2e12aea085` and source digest
`9b95495a6d00699892c23d02065430034a7d0930f6e8913323571bf93b124a56`.
The recorded dirty-worktree warning requires use of the source digest rather
than Git head alone; manuscript files outside the run source scope were present.

## Result

All 21 discovery runs completed without a failure or timeout. The pre-specified
rule first declared operational catalogue saturation at 20x: the 5x--10x and
10x--20x transitions passed in all three seeds, and 20x--30x also passed.

All 55 planted families were recovered in every seed at 5x, including the
20-copy tier, and no additional planted family was recovered at higher depth.
Median operational-family count was 34, 46, 53, 55, 56, 57 and 58 across the
seven depths. Median unique-candidate count continued to rise from 36 at 0.5x
to 119 at 30x. Thus 5x truth recovery and 20x operational saturation answer
different questions. A small number of higher-depth operational families were
not additional planted families.

This is controlled simulation evidence from one factorial design and one IID
error model. It does not establish a universal 20x requirement for real plant
libraries or guarantee recovery of rare unmodelled repeat families.

## Archive contents

- `depth_metrics.tsv`: one row per seed and depth;
- `family_recovery.tsv`: every planted family at every seed and depth;
- `transition_metrics.tsv`: per-seed adjacent-depth yield and stability;
- `saturation_decision.json`: primary all-seed decision;
- `run_summary.tsv`, `run_receipt.json`, `evaluation_receipt.json`: completion
  and resource receipts;
- `environment.json`: source, input and execution provenance;
- `figures/`: publication curve, editable SVG, PDF, PNG, legend and panel-level
  source data;
- `config.json`, `protocol.md`: frozen design copies.
