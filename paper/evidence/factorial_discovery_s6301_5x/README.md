# Factorial discovery development pilot

Source 57ee82d8c78eb9ced1be135068c8c99567ca723e; source genome6301, 5x clean
and high-error observed-read conditions. Each has3,619 reads, approximately50 Mb
and55 planted founder families. All six actual tool executions and independent
scoring stages completed. This is one simulated source genome with two paired
technical conditions, not six replicates or real-species validation.

| Condition | Tool | Seconds | Peak MiB | Eligible array recall | Eligible array precision | All-input base recall | All-input base precision |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Clean | TandemX | 50.706 | 120.27 | 1 | 1 | .999565 | .999950 |
| Clean | TRF | 46.873 | 183.98 | 1 | 1 | .999631 | .999988 |
| Clean | TideHunter | 14.434 | 345.17 | 1 | .939123 | .999471 | .994767 |
| High error | TandemX | 58.898 | 119.91 | 1 | 1 | .999531 | .999952 |
| High error | TRF | 81.786 | 228.31 | 1 | .698867 | .999443 | .999988 |
| High error | TideHunter | 17.137 | 388.42 | 1 | .920175 | .999403 | .996805 |

Every tool recovers all55 founders under the independent cyclic edit threshold
(both all-planted and observed-eligible denominators). These conditions therefore
do not establish a family-recall advantage. TRF's lower raw array precision in
high-error reads coexists with very high base-union precision; duplicate/overlapping
reports are distinct from biologically wrong sequence. TandemX's base precision
is slightly lower than TRF's and its speed remains lower than TideHunter's.

All timings overlapped other development/acquisition jobs: diagnostic observations,
not final resource rankings. See `inputs.tsv` for complete truth-fragment exclusions,
`summary.tsv` for every metric, per-tool family/array assignments and source receipts.
Array/read intervals are descriptive counts of correlated evidence. Raw output
and source snapshots stay at the T7 paths in `archive_manifest.json`. No universal
superiority, advanced array reconstruction or mature release is claimed.
