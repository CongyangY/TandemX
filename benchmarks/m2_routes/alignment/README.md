# M2 route A: local monomer-path dynamic programming

Research prototype, frozen before receiving the B1 v2 development bundle. No
production CLI behavior changes. The input is a pretrimmed array interval or a
B1 v2 case containing exact unique engineered left/right flanks. The route
never establishes biological same-locus anchoring, molecule independence,
haplotype identity, or an assembly error from sequences alone.

The segmenter tiles the full array with supplied candidate monomer templates
and their reverse complements. Each candidate segment has an Edlib global edit
cost. A bounded dynamic program retains up to three distinct label/orientation
paths per endpoint. Local template ambiguity and a competing global path
within two edit costs cause abstention. The read paths are compared to the
assembly path by label-level edit distance. All optimal label alignments define
an interval of possible change locations; repetitive arrays can yield a broad
interval. These are edit costs and binary decisions, **not calibrated
probabilities**.

Frozen development parameters:

| Parameter | Value |
| --- | ---: |
| Minimum informative independently anchored reads | 3 |
| Minimum dominant path fraction | 0.75 |
| Minimum resolvable read fraction | 0.75 |
| Minimum local/global label margin | 2 edit costs |
| Maximum per-monomer edit fraction | 0.15 |
| Maximum monomer length variation | 0.08 |
| Maximum candidate alignment count per sequence | 250,000 |
| Maximum pretrimmed array bases | 4,096 |
| Maximum candidate monomers | 16 |

The adapter treats `event_score=1.0` for `DISCORDANT` and `0.0` for `SUPPORTED`
as an uncalibrated binary decision at the common 0.5 threshold. `AMBIGUOUS` and
`INSUFFICIENT_READ_SUPPORT` become `abstain` with null score. Any other flank
policy, missing or nonunique flanks, invalid input, or unverified pairing
abstains. All case rows are retained. The route does not estimate missing bases
from copy counts and cannot distinguish all biological versus technical causes
of read--assembly disagreement.

```bash
conda run -n tandemx-dev python -m benchmarks.m2_routes.alignment.adapter \
  --inputs /path/to/inputs.jsonl --predictions /path/to/predictions.jsonl
conda run -n tandemx-dev pytest -q tests/unit/test_m2_alignment_route.py
```

Known limitations: candidate catalogue dependence; exact engineered flanks;
no cyclic monomer phase search; no intra-label variant calls; no posterior
calibration; no mixture deconvolution; no whole-array megabase scaling. The
adapter is suitable for B1 synthetic development cases. A real-data result
requires independent native anchor and same-donor/haplotype validation.

## Frozen B1 v2 development input-only run

The method parameters and tests were fixed before reading the B1 v2 input.
Input SHA256:
`7b27e7cdb5b1a0d15a97e6ac08f5947de028c90096d644a8c927a07cc532383a`.
Prototype SHA256:
`043122d1b5987e35dcd957a6471c930f6211fe7a8902055d3a8c37f06d4ed651`.
Adapter SHA256:
`6d601b28746276426b3191768ef2b3150acf5101451f6272421b35715e9a982f`.
Prediction SHA256:
`cac9c50a50d60217f747ae51332f21d42ee5cf3a2d10911797ff75b75380a5d3`.

The 13 opaque development inputs produced 10 `DISCORDANT`, 2 `SUPPORTED`,
and 1 `AMBIGUOUS`/`abstain` outcomes. Summed per-case elapsed time was 0.157 s;
peak process RSS was 17,301,504 bytes on the local macOS host with the
`tandemx-dev` interpreter. This tiny single-lineage run is not a scale
benchmark. The abstaining edited interval contains 352 bases against 24-base
candidate templates, so the full-tiling assumption fails. It was retained and
the algorithm was not adjusted after seeing the input. Truth, sensitivity,
precision and false-positive metrics are assigned by the independent B1
scorer, not inferred from the above call counts.
