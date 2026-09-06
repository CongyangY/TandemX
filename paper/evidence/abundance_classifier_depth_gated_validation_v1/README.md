# Frozen depth-gated classifier validation

This archive preserves the complete compact evidence chain for the transparent
depth-gated v3 abundance classifier. The rule was developed after the failed
5701-5703 test, frozen in commit `62892a6`, and passed hosted Ubuntu/macOS CI
before seeds 5801-5803 were generated once.

The rule uses the single-k21 estimate and decision threshold 0.6 when the
single-k21 estimated haploid depth is below 2. At estimated depth at least 2, it
uses the alpha-0.5 log-space single/multi-k blend and decision threshold 0.5.
Unavailable or nonpositive multi-k estimates fall back to single k=21. No such
fallback occurred among the 1,620 standard-depth held-out rows.

## Result

| Method | TP | FN | FP | TN | Sensitivity | False-positive rate | Precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Single k=21 | 867 | 591 | 16 | 956 | 0.594650 | 0.016461 | 0.981880 |
| Depth-gated v3 | 953 | 505 | 16 | 956 | 0.653635 | 0.016461 | 0.983488 |

The frozen rule increased sensitivity by 0.058985, left the false-positive rate
unchanged and increased precision by 0.001608. Sensitivity improved in every
held-out seed; the minimum seed-level gain was 0.030864. The maximum seed-level
false-positive-rate delta and minimum precision delta were both zero. All six
predeclared held-out gates passed.

The same baseline run retained full-assembly mean base recall 0.980059,
positive-assembly mean base precision 0.999617 and zero predicted bases in 54
absent-family rows. The weakest displayed stratum was 5% planted monomer
substitution with three segments, with mean base recall 0.951982.

## Evidence map

- `development/`: six-genome post-failure development evidence and frozen rule.
- `baseline/`: one-time 5801-5803 single-k21/localization/comparison outputs,
  all 1,062 command receipts consolidated in `resource_metrics.tsv`, and selected
  source snapshots.
- `multik/`: 4,860 raw single-k21/multi-k rows for 2,430 paired conditions.
- `heldout/`: frozen depth-gated rows, overall and per-seed metrics, acceptance
  receipt and the explicit `no_heldout_fit_or_selection=true` environment field.
- `source_snapshot/`: benchmark scripts required to audit development, multi-k
  and final held-out evaluation.
- `archive_manifest.json`: 72 copied-file records with source path, byte count
  and SHA-256.
- `figures_v1/`: inspected first figure draft, retained for provenance.
- `figures_v2/`: accepted SVG, PDF and PNG, 46-row panel source table and output
  hashes. The SVG contains 126 editable text nodes and zero raster image nodes;
  both the PNG and a 180-dpi render of the PDF were visually inspected.

The 1,062 sequential direct-child receipts sum to 143.257 seconds: 19.487 seconds
for 90 locate commands, 27.697 seconds for 162 quantify commands and 96.073
seconds for 810 compare commands. Maximum recorded child peak RSS was 32.797
MiB. These are small synthetic command-level diagnostics and do not establish
large-genome speed or memory superiority.

This is a known-catalogue simulation with IID read substitutions, three new
simulated genome seeds and no biological collapse truth. It validates the
predeclared rule under this matrix; it is not biological validation, an AI
result or evidence of universal superiority over external tools.
