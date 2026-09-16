# Competitive occupancy research gate: no public backend

Status: **no-go for public `tandemx quantify` integration**. This remains a
research-only streaming prototype under `benchmarks/m1_shared_signature/`.
The uncommitted CLI, README and public format/algorithm edits made during
exploration were removed. The frozen diagnostic-k-mer default, its output
schema, pipeline and compare consumers were never changed in a commit.

The first attempt used a two-row Python semiglobal edit DP, non-overlapping
80-bp windows, a 12% edit gate and a two-edit family margin. On seed 11 of
the four frozen M1 conditions, positive-family MARE was 0.0314, 0.5250,
0.5194 and 1.0000, versus ordinary mapping 0, 0.0201, 0.3877 and 1.0000.
It took about 7–8 s per 600 80-bp reads. The original four-case receipt,
including its then-current source hashes, is retained as
`evidence_occupancy_research_20260916/initial_python_dp_seed11.json`.
That pre-relocation source snapshot was not committed; this receipt is a
historical development record rather than a claim of byte-identical source
replay. The Python fallback and first-attempt gates remain available in the
current research module for bounded rechecks.

One revision was allowed. It used `edlib` for the same periodic semiglobal edit
distance, aligned the 80-bp gate to the baseline's 18/80 mismatches, used a
unique-best one-edit margin, and made catalogue pairs within one periodic edit
ambiguous by input design. It then replayed all 12 frozen catalogue/read/truth
hashes (three seeds × four conditions) without touching a final holdout.

| Condition | Ordinary mapping mean positive MARE | Revised occupancy mean positive MARE | Revised wrong-family calls, 3 seeds | Mean ambiguous bp / 48,000 | Mean unknown bp / 48,000 | Mean fit s / 600 reads |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Related, 3% substitutions | 0.0034 | 0.0034 | 1 | 80 | 8,000 | 0.0206 |
| Related, 12% substitutions | **0.0144** | 0.0164 | 7 | 613 | 8,080 | 0.0211 |
| `f2` extra 20% substitutions | 0.3796 | 0.3658 | 5 | 560 | 17,920 | 0.0209 |
| Identical `f1`/`f2` | 1.0000 | 1.0000; individual refusal | 0 | 40,000 | 8,000 | 0.0218 |

The revised model made zero calls from the 100 random negative reads per case
and zero calls to the planted-zero decoy. The strong ordinary baseline also
assigned zero to that decoy. At 12% error, occupancy was worse on the primary
error and made seven wrong-family calls across three seeds, matching the
ordinary mapping wrong-family total of seven. Under family-specific error it
improved mean MARE by 0.0138 but made five wrong-family calls versus the
ordinary baseline's three. This mixed development result does not establish
an accuracy or specificity improvement. The revised run's process peak
`ru_maxrss` was 35,700,736 bytes; its ~0.02 s/case speed depends on optional
`edlib` in `tandemx-dev`. The same objective has a much slower pure-Python
fallback. These are toy-process observations, not whole-pipeline or GB-scale
benchmarks. Full per-case provenance, source hashes, timing and confusion are
in `evidence_occupancy_research_20260916/revised_12_case/results.json`.

The decisive scope limit remains the algorithm itself: only eight catalogue
families of at most 1,024 bp are allowed, every 80-bp window is compared to
every family on both strands, and a read is segmented at fixed window
boundaries. Partial arrays, family transitions, uneven read length, repeats
in genomic background and platform-specific indels have no calibrated
exposure/error model. The engine stores one read and two DP rows at a time,
which bounds memory but does not establish production throughput. Its
`assigned_read_bp / haploid_depth` is provisional occupancy, not physical
family copy-number truth or an assembly deficit estimate. Identical and
nearly identical units must remain ambiguous.

The unit test exercises a long single-family read and a synthetic one-base
indel. It does not validate quantitative behavior on long reads containing
mixed repeat families, array boundaries and background, or real HiFi/ONT
error profiles. The small family-error-shift MARE gain is retained in the
table but is not selected as a method, because the other development
condition worsened and specificity did not improve consistently.

Decision: retain the prototype, tests, initial failure and revised replay for
audit; do not expose it through the public CLI, pipeline or `copy_number.tsv`.
There is no further threshold search in this M1 round. Any future public
backend needs a distinct, preregistered full-read/array method, catalogue
scaling, independent material validation and resource gates.
