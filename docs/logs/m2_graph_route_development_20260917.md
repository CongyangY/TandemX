# M2 transition-graph route: development-only record

The independent route in `benchmarks/m2_routes/graph/prototype.py` compares
directed, multiplicity-aware monomer transitions between two verified flanks.
It starts with supplied, ordered monomer assignments. It is not the prior
sequence-to-label dynamic program and does not infer assignments or flank
anchors from raw reads. Its status vocabulary is `SUPPORTED`, `DISCORDANT`,
`AMBIGUOUS`, and `INSUFFICIENT_READ_SUPPORT`. `SUPPORTED` means only that the
dominant transition graph matches the assembly graph; it is not proof that the
assembly is correct. A separate adapter emits the common controlled-edit
prediction schema, with uncalibrated binary 0/1 event scores and null bp/path
fields. Abstentions retain null scores.

## Frozen development defaults

- At least three distinct molecule IDs and two distinct named error profiles
  among the dominant graph evidence.
- Both flanks verified and every upstream monomer assignment confidence at
  least 0.9. The route does not calibrate those confidence values.
- A dominant graph must exceed two thirds of qualified reads. At least two
  molecules supporting a minority graph trigger abstention.
- Catalogue sequences are required for a discordance call. If any involved
  monomers are within normalized edit distance 0.15, the label identity is
  ambiguous. These cutoffs were selected for this component trial and are not
  validated on biological reads.
- Multiple supplied haplotype IDs or conflicting full paths with the same
  transition spectrum trigger abstention. The graph does not resolve order
  from identical edge multisets.

## Reproduction and observed result

In `tandemx-dev`:

```bash
python -m benchmarks.m2_routes.graph.run_development \
  --output docs/logs/m2_graph_development_20260917.json
pytest -q tests/unit/test_m2_transition_graph.py
```

Thirteen focused tests passed. The 14 engineered-label cases returned three
`DISCORDANT`, three `SUPPORTED`, four `AMBIGUOUS`, and four
`INSUFFICIENT_READ_SUPPORT`. The three graph-positive examples are copy
multiplicity, transition rewiring, and replacement by a distinct supplied
monomer class. These are direct positive controls, not sensitivity estimates.
An order change with exactly the same transition multiset (`ABACA` versus
`ACABA`) abstained. A one-base monomer subclass difference abstained because
label identity is not separable at the preset cutoff. A same-profile correlated
error scenario abstained. Low support, low-confidence substitution/indel/
homopolymer representations, and a mixed-haplotype scenario did not become
discordance calls. One rare discordant molecule among four matching molecules
was classified graph-supported and retains a nonzero discordant-read count;
this result does not resolve whether it is a rare variant or an error.

The highest single-case elapsed time was under 0.002 s and the largest
`tracemalloc` Python allocation peak was under 100 KiB on these tiny paths.
Exact timings, peaks, all statuses, source hashes and generated input hash are
in the JSON receipt. `tracemalloc` is not process RSS, and 8-bp motifs with at
most six copies cannot predict runtime or memory on centromeric arrays. The
upstream quality/profile metadata are supplied by construction; this trial
does **not** demonstrate robustness to raw sequence errors. A within-label SNP
is invisible to this graph and produces graph support. Systematic errors can
still cross the named profiles; the profile rule reduces one confound but does
not establish sequencing-error independence.

This route has not yet consumed the common B1 sequence-only development inputs,
any realistic noisy reads, a held-out split, or real same-locus molecules. It
must not be called an assembly-error detector from this receipt. All status
and error-type denominators for the planned tournament must be computed from
the frozen common bundle, including abstentions.
