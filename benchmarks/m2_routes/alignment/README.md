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
