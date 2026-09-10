# Independent audit of SRF formal unified run (2026-09-10)

Input: `/Volumes/T7/Codex/TandemX/results/srf_formal_unified_v1_20260910`. Read-only audit after `completion.json` reported 72 result cells.

No files on T7 were modified and no native experiment was rerun. Metrics below were recomputed from each cell’s `family_abundance.tsv`, `result.json`, `srf.bed`, `srf.abundance.tsv` and `native_abundance.json`.

## Completion and status

`completion.json`: `{
  "completed_cells": 72,
  "expected_cells": 72,
  "unrun_cells": 0,
  "time_limit_reached": false,
  "failures": 0,
  "elapsed_seconds": 174.33619804214686
}`

| method | cells | status counts | positive-family rows | MARE all rows | MARE successful only | recovered / truth (all cells) | negative predicted / negative total | aggregate negative fraction |
|---|---:|---|---:|---:|---:|---|---:|---:|
| competitive_mapping | 18 | ok=18 | 54 | 0.000304977 | 0.000304977 | N/A; shared catalogue matches 3/3 | 33497 / 3600000 | 0.009304722 |
| srf_k101 | 18 | ok=18 | 54 | 0.068204789 | 0.068204789 | 2/3, 3/3 | 682 / 3600000 | 0.000189444 |
| srf_k151 | 18 | no_catalogue=3, ok=15 | 54 | 0.208782597 | 0.050539116 | 0/3, 3/3 | 0 / 3600000 | 0.000000000 |
| tandemx | 18 | ok=18 | 54 | 0.102612125 | 0.102612125 | 3/3 | 92660 / 3600000 | 0.025738889 |

## SRF bed2abun versus native totals

For each successful SRF cell, I independently summed BED half-open span lengths, summed the second column of `srf.abundance.tsv`, and summed `native_abundance.json`.

| method | successful cells | BED sum vs abundance sum mismatches | BED/native mismatches | maximum absolute difference |
|---|---:|---:|---:|---:|
| srf_k101 | 18 | 0 | 0 | 0 |
| srf_k151 | 15 | 0 | 0 | 0 |

## Interpretation and limits

- The run is complete: 72/72 cells, zero failures, zero unrun cells.
- “MARE all rows” includes positive-truth family rows from cells labelled `no_catalogue`; “successful only” excludes those status cells. This preserves the distinction between a zero estimate and a successful native workflow.
- “Recovered / truth” is copied and cross-checked from per-cell family denominators; it is a family recovery endpoint, not a physical copy-number truth claim.
- Negative fraction is recomputed as pooled `sum(negative_read_predicted_bp) / sum(negative_read_total_bp)`, rather than averaging per-cell fractions.
- BED span totals, `srf.abundance.tsv` totals and `native_abundance.json` totals agree in every successful SRF cell; this checks internal endpoint arithmetic only and does not validate SRF biological accuracy.
- No formal SRF superiority or generalization claim follows from this audit.

## Unified 36-cell accounting for TandemX and competitive mapping

The requested non-SRF controls were also independently checked across all 18
TandemX and 18 competitive-mapping cells (36 cells total). The preferred main
MARE is the pooled mean over every positive-truth family row, including rows
from any cell whose native workflow returned no catalogue; the successful-only
value is retained only as context.

| method | cells | positive-family denominator | pooled MARE (preferred) | matched truth families / truth families | negative predicted / negative total | pooled negative fraction | native catalogue total | unassigned native bp | unassigned abundance bp | abundance excess bp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TandemX | 18 | 54 | 0.102612125 | 54 / 54 | 92,660 / 3,600,000 | 0.025738889 | 259 | 92,407 | 31,584 | 0 |
| Competitive mapping | 18 | 54 | 0.000304977 | 54 / 54 | 33,497 / 3,600,000 | 0.009304722 | 259 | 33,497 | 33,024 | 544 |

The native catalogue totals include 15 three-family cells, two 66-family
shared-fragment cells and one 82-family shared-fragment cell for each method.
The correspondence records were retained separately; unmatched native bases
remain in global accounting and are excluded from mapped family MARE. These
are endpoint bookkeeping results, not a claim that competitive mapping is
biologically superior.

## Unmatched native output retained separately

The unified interval records keep unmatched native output outside the mapped
family MARE numerator while retaining its bases in global prediction accounting.
This matters for the shared-fragment condition: TandemX native catalogue counts
were 66, 82 and 66, with `unassigned_native_bp` 29,847, 35,717 and 26,843 for
seeds 2026091001, 2026091002 and 2026091003; the corresponding negative-read
prediction totals were 29,847, 35,717 and 27,096 bp. Competitive mapping had
the same catalogue counts and `unassigned_native_bp` values 10,966, 19,505 and
3,026 bp. These unmatched masses are not silently discarded, and they are not
included as evidence of family-level MARE or biological error.

Independent de novo recovery for competitive mapping is N/A. Its 54/54 matched
truth count describes the supplied TandemX catalogue and is not a second
independent discovery result.
