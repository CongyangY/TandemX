# B1 v3: Col-CEN interval development cases

This is a distinct **development** dataset. Three 3,560 bp native assembly
intervals come from the [official Col-CEN v1.2 FASTA](https://github.com/schatzlab/Col-CEN/tree/main/v1.2),
pinned at commit `abb9b614d91c8a0bbd05a199c694fb5eafb5fe30`. The complete
compressed FASTA SHA-256 was checked against the archived source receipt:
`b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8`.
`source_arrays.fa` is a small exact extract; `source_catalogue.json` records
per-interval coordinates, SHA-256 and period checks. The selected intervals
are Chr1 `[15585019,15588579)`, Chr2 `[3965542,3969102)`, and Chr3
`[15060000,15063560)` in 0-based Col-CEN v1.2 coordinates. They avoid the
two known issue intervals listed in `paper/evidence/cohort_screen/Col-CEN_v1.2_issues.bed`.

All three intervals show strong 178 bp periodicity: the minimum identity
between adjacent operational tiles is 0.916, 0.933 and 0.910, respectively.
The 20 tile boundaries per interval were chosen as an operational phase by
maximizing local 178 bp similarity. **Native biological monomer boundaries and
HOR periods have not been independently established.** The selected windows
are deliberately favorable periodic segments and are not a representative
sample of the entire centromere or three independent donors. All three come
from one Col-CEN assembly lineage.

The development generator makes 13 predeclared edits for each source interval,
including two intact controls, terminal and internal loss, complete deletion,
boundary truncation, block/one-tile deletion and duplication, rearrangement,
inversion and per-tile compression. The edited sequence, source and edited
coordinates, source tile mapping and bp delta are exact injected-edit truth.
The full-flank raw-read entries are **three identical exact synthetic copies**
of the unedited source path under separate IDs. Flanks are engineered unique
64 bp sequences. These IDs are not independent molecules or error profiles.
Original Col-0 raw reads are not enrolled here; source accession alone does
not establish extraction, donor, or haplotype pairing. Physical missing-copy
or biological read–assembly accuracy remains blocked.

`protocol.json` and `metrics_frozen.json` were fixed before route predictions.
The public `development_bundle/inputs.jsonl` retains schema version 2 and
includes only assembly sequence, synthetic read sequences, three per-array
178 bp consensus candidates, engineered flanks, source interval metadata and
input hashes. The separate `truth.jsonl` contains causal event labels and
ledger. Both source intervals and results have SHA-256 receipts. This is an
audit split, not secure blinding. Do not tune a route after viewing its truth
and then report its performance as independent validation. There is no held-out
set in v3.

The all-case intent-to-diagnose denominator is 39: 33 injected positives and
six intact controls. An `ok` score ≥0.5 is positive. A positive `abstain`,
`failed`, `not_run`, `unmapped` or missing row counts as FN; an intact non-`ok`
row is `unresolved_negative` and contributes to the negative failure rate.
FPR uses all six intact controls and counts only explicit positive calls as
FP. Coverage and eligible-only metrics are separate. Event subtype and
breakpoint localization are N/A when exact edited sequences have different
causal truths. HOR period error and native repeat-order accuracy are N/A here
because operational 178 bp tiling is not independently verified native unit
truth. Runtime and peak RSS require an external measured route receipt; an
absent receipt yields N/A, not zero. No AUPRC, confidence interval, or
biological collapse accuracy is reported.

Reproduce inside `tandemx-dev`:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.extract_b1_colcen_arrays \
  --source-gz /path/to/verified/Col-CEN_v1.2.fasta.gz \
  --protocol benchmarks/controlled_collapse/v3/protocol.json \
  --out-fasta /tmp/b1v3_source_arrays.fa \
  --out-catalogue /tmp/b1v3_source_catalogue.json
conda run -n tandemx-dev python -m benchmarks.scripts.build_b1_colcen_cases \
  --protocol benchmarks/controlled_collapse/v3/protocol.json \
  --source-fasta benchmarks/controlled_collapse/v3/source_arrays.fa \
  --source-catalogue benchmarks/controlled_collapse/v3/source_catalogue.json \
  --outdir /tmp/b1v3_replay
conda run -n tandemx-dev pytest -q tests/unit/test_b1_colcen_cases.py
```

To score a route after it writes a prediction JSONL, use:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.score_b1_colcen_cases \
  --bundle benchmarks/controlled_collapse/v3/development_bundle \
  --metrics benchmarks/controlled_collapse/v3/metrics_frozen.json \
  --predictions /path/to/predictions.jsonl --outdir /tmp/b1v3_score
```

The 33 TP/6 TN unit-test fixture reads the truth ledger directly to verify
scorer arithmetic. It is a tautological smoke test, not method performance.

## Frozen development run

The v3 source, protocol, metrics, scorer and tests were committed as
`931f0e9` before predictions were scored. The input-only route predictions
were saved under `benchmarks/m2_routes/`; `score_audit.json` fixes their hashes
and the scoring receipt hashes. `route_comparison_per_case.tsv` retains every
case, including refusals. Full repository `pytest -q` passed 1,023 tests
before the scoring pass.

| Route | `ok` / 39 | TP | FN | FP | TN | unresolved intact | Intent sensitivity | Core time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Alignment | 36 | 24 | 9 | 0 | 6 | 0 | 24/33 | 164.211 s |
| Graph | 0 | 0 | 33 | 0 | 0 | 6 | 0/33 | 0.00196 s |
| Simple flank baseline | 39 | 18 | 15 | 0 | 6 | 0 | 18/33 | 0.00876 s |

The graph time measures immediate refusal, not completed inference speed.
These are route core times from different implementations, excluding process
startup and output writing; they are descriptive, not a normalized runtime
comparison. FPR is zero for all three, but graph has a 6/6 negative failure
rate because it refuses all intact controls. Alignment eligible-only
sensitivity is 24/30; its all-case intent sensitivity is 24/33. Alignment
detects inversion and one-tile loss in these cases; the simple flank baseline
misses both classes. Both routes miss all three rearrangements and all three
within-tile compression cases. Alignment explicitly abstains on three of six
duplications. The case families share source intervals and edit rules,
so the 39 rows are not 39 independent biological replicates.

No route supplies complete exact event-type or edited-breakpoint predictions,
so those aggregates remain N/A. The simple flank baseline reports zero signed
bp MAE on exact synthetic full-flank pairs, which is deterministic length
arithmetic here and does not establish raw-read copy-number accuracy.
