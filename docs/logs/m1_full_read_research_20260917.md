# M1 full-read competitive assignment: bounded development result

Status: **research-only; no public `quantify` backend promotion**. The
candidate, protocol, generator, 18-read input, truth ledger, and file hashes
were committed as `9346de2` **before** the first scoring run. This is a new
development set; no B2 or previously reserved M1 final holdout was used.
Evidence: `benchmarks/m1_shared_signature/evidence_full_read_dev_20260917/`.

The candidate streams complete reads and classifies disjoint 40-bp core tiles
with 40-bp flanking context against five supplied periodic 80-bp units. It
requires a unique best family in both core and context. Conflicts and the
exactly identical `twin_a`/`twin_b` units abstain; unknown mass is retained.
This removes the unit-length input rejection of the earlier ordinary mapping
prototype, while retaining a fixed-tile boundary limitation. It is still an
occupancy prototype, not a physical-copy model.

Run the frozen comparison with:

```bash
conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.run_full_read_dev --outdir /tmp/m1_full_read_dev_replay
conda run -n tandemx-dev pytest -q tests/unit/test_m1_full_read_research.py
```

The development set contains three technical read repetitions each of pure
`f1`, mixed `f1`/`f2` plus background, reverse-oriented `f2`, an identical
family pair, background-only, and mixed families with substitutions and
indels. It has 18 reads and 11,161 bp: 3,604 `f1`, 3,599 `f2`, 1,200
`twin_a`, and 2,758 background bp. The zero-abundance `decoy_zero` is in the
catalogue. These are not independent donors, real read errors, or a
representative full genome.

| Method | All-positive MARE, refused family = 100% error | Assigned `f1`/`f2`/`twin_a` bp | Wrong assignment bp | Background false bp | Zero-decoy bp |
| --- | ---: | --- | ---: | ---: | ---: |
| Full-read research candidate | 0.3410 | 3560 / 3560 / 0 | 7 | 3 | 0 |
| Ordinary mapping, independent 80-bp chunks | 0.5928 | 2800 / 1600 / 0 | 4 | 0 | 0 |
| Current diagnostic-k-mer quantify, explicit depth 1 | N/A: different endpoint | 1600 / 1680 / 0 estimated genomic bp | N/A | N/A | 0 |

The candidate assigned 7,120 bp, marked 1,200 bp ambiguous, and marked
2,841 bp unknown. Of that unknown mass, 2,755 bp was background and 86 bp
was positive-family sequence. The 80-bp diagnostic adapter marked 2,758
background bp and 2,803 positive-family bp unknown. All 1,200 `twin_a` truth bp were ambiguous; none were
given an individual family estimate. Thus the lower aggregate MARE compared
with chunked mapping does **not** establish better family attribution: the
candidate made more wrong-family assignments (7 versus 4 bp). The ordinary
baseline has larger losses at mixed-read phase boundaries. The production
quantifier has no per-base calls, so its wrong-assignment fields are N/A.

The production row uses `--haploid-depth 1` solely to expose its estimated-bp
field on the same FASTA. That depth is an artificial convention, not measured
single-copy depth or physical calibration. Its MARE is **not computed**:
candidate/baseline family bp are assigned read occupancy, while `quantify`
reports an inferred genomic bp quantity. The production outputs and their
confidence/warning fields remain in the raw worker receipt. Thus the two
occupancy MAREs cannot be used to rank the production estimator.
No abundance-to-assembly deficit was tested.

The three workers ran as separate processes on the same input. Exact raw
wall times, process RSS, and worker receipts are in `results.json`. The candidate times
one streaming inference pass and generates a scoring trace afterward.
Chunked mapping times one pass including its trace, and production quantify
times its whole function including output writing. This toy, unequal work
scope makes the resource **comparison N/A**; the recorded values are
descriptive only. The 80-bp chunked ordinary mapper is a diagnostic adapter,
not a production full-read backend. The
candidate retains one read plus counters in its inference path, bounded by
the explicit 1-Mbp read and eight-family caps; its all-family, both-strand
alignment cost has not been measured on a realistic catalogue.

To promote a public backend, first freeze a complete-read method with
validated boundary behavior and group abstention, then test an independent
donor/lineage set with an externally justified depth calibration. Evaluate
physical copy-bp error and assembly-deficit magnitude separately from binary
collapse calls. Compare all-family errors, wrong attribution, zero-decoy
abundance, uncertainty coverage, throughput and RSS on matched inputs with
no retuning on final holdouts. The present result is a feasibility observation
and an adverse record, not passage of that gate.
