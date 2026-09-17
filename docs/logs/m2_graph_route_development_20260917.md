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

Fifteen focused tests passed after adding the exact synthetic B1 adapter. The
14 engineered-label cases returned three
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

At this initial component check, the route had not consumed the common B1
sequence-only development inputs, any realistic noisy reads, a held-out split,
or real same-locus molecules. It must not be called an assembly-error detector
from this receipt. All status and error-type denominators for the planned
tournament must be computed from the frozen common bundle, including
abstentions.

## Frozen B1 input-only development run

After the B1 input-only bundle and metric contract were frozen, a separate
strict adapter consumed its 13 opaque public rows without reading `truth.jsonl`.
It accepted only the declared exact synthetic error profile, required one
occurrence of each engineered flank, and split the interior into exact 24-bp
catalogue tiles on either orientation. This exact segmentation is possible by
construction; it is not a general noisy-read decomposition method. Because
the declared reads have zero engineered sequencing error, this run allowed one
named `exact_synthetic` profile. The ordinary two-profile requirement remains
in effect for other data. One partial-tile assembly interval abstained.

Command:

```bash
python -m benchmarks.m2_routes.graph.run_common_development \
  --inputs benchmarks/controlled_collapse/v2/development_bundle/inputs.jsonl \
  --predictions benchmarks/m2_routes/graph/common_development_predictions.jsonl \
  --audit docs/logs/m2_graph_common_development_20260917.json
```

The input-only outputs are 10 `DISCORDANT`, two `SUPPORTED`, and one
`COMPONENT_UNRESOLVED`. These counts are not accuracy because the sealed truth
was not inspected in this step. The prediction rows retain every case ID and
the abstention. No exact event type, edited bp count, or full label path is
claimed. The single abstention blocks the frozen protocol's all-case primary
aggregate; a separate selective-coverage descriptive table may be computed by
the scorer but cannot replace the 13-case denominator. The detailed read and
assembly copy counts, graph transitions, timing and Python allocation peaks
are in the audit receipt. The maximum per-case elapsed time was 0.0062 s and
the maximum `tracemalloc` peak was 73,179 bytes on this exact 13-case bundle;
neither is a process-RSS measurement or a scale benchmark.
