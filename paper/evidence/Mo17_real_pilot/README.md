# Mo17 random-library development comparison

All inputs derive from the QC-passed complete SRR15447419 batch by seeded
read-ID hash sampling. `baseline` and `native_audit` use the identical 840 reads /
11,680,888 bp; `native_audit_111Mb` uses 8,084 reads / 111,505,681 bp from the
same nested sampling series. This is one technical batch from Mo17, not three
independent materials. Exact input/tool/source hashes are in each environment.

| Input and TandemX audit | TandemX seconds / MiB | TRF seconds / MiB | TideHunter seconds / MiB |
| --- | ---: | ---: | ---: |
| 11.681 Mb, original Python audit | 258.829 / 168.58 | 11.859 / 194.72 | 3.508 / 147.80 |
| Same input, native full audit | 17.348 / 81.70 | 12.019 / 165.17 | 3.644 / 155.03 |
| 111.506 Mb, native full audit | 247.539 / 165.47 | 87.493 / 212.52 | 34.115 / 324.36 |
| Same 111.506 Mb, exact-related audit | 142.919 / 155.45 | 83.620 / 210.84 | 32.653 / 269.34 |

Single one-thread runs overlap acquisition/QC. Memory is direct-child wait4
peak RSS, excluding the controller. They are diagnostic observations, not
isolated repeated resource rankings. Counts, called-base unions and agreement
are not accuracy without an independent denominator. TandemX is slower than
both tools here, despite its substantial audit optimization. All unfavorable
results remain visible. All twelve tool executions and normalizations succeeded.

Six data files are byte-identical before/after the native audit on the 11.681-Mb
input (`output_parity.json`). The 111.506-Mb catalogue contains 4,219 families,
causing 8,897,871 exhaustive audit rows; raw outputs remain on T7. These data
motivated the selectable exact-related policy. Its real rerun preserves all five
primary data files and all 4,595 non-distinct rows byte-for-byte. Only 10,471
pairs require exact alignment; the audit phase takes approximately 0.69 s.
`related_audit_parity.json` and `related_audit_111Mb` retain these checks, all
commands/results and the explicit omitted-pair count. The smaller pair output
is intentional; it is not a claim that the full table remains byte-identical.
