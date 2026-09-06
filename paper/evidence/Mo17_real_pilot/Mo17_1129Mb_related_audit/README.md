# Completed 1.129-Gb Mo17 development comparison

Source b53a193. Identical whole-library hash sample: 81,775 reads and
1,128,793,699 bp. All three tool executions and normalization completed.

| Tool | Seconds | Peak MiB | In-scope calls | Positive reads | Union bp |
| --- | ---: | ---: | ---: | ---: | ---: |
| TandemX | 1930.701112 | 710.593750 | 85663 | 51534 | 29573588 |
| TRF | 1019.346266 | 323.515625 | 112081 | 46389 | 26672458 |
| TideHunter | 326.371856 | 390.546875 | 117816 | 57489 | 38562423 |

TandemX retained28,586 operational families; related-mode audit scored259,411
of408,565,405 possible pairs and emitted100,287 related rows. Scanning took
1269.837 s and clustering approximately643.8 s from stage logs. These counts do
not establish family truth, accuracy or missing repeat abundance.

TandemX is slower and uses more measured peak memory in this run. All unfavorable
values are retained. Resource values are one-thread direct-child diagnostics
under concurrent downloads/other jobs, not final publication rankings. This
source predates the exact indexed multiset clustering gate. The newer 1-Gb run
has its own source/output directory and must not overwrite this baseline.
