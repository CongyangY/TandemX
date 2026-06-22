# Rice Known-repeat Validation Notes

This document describes post hoc interpretation for bounded rice HiFi pilot
runs. It does not define inputs for `tandemx discover`.

## Boundary

The required order is:

```text
rice HiFi reads -> tandemx discover -> monomers.fa -> tandemx annotate-repeats
```

Known rice repeat sequences must not be passed to `tandemx discover`. They are
used only after de novo discovery to interpret candidate families and evaluate
whether a discovered monomer resembles published repeat classes.

## Suggested known-repeat library contents

For rice tests, prepare a FASTA library from documented sources where licensing
and provenance are clear. Useful categories include:

1. CentO / CEN155 centromeric satellite sequences;
2. TrsA / Os48 tandem repeat sequences;
3. X55642.1-like 354/358 bp tandem repeat sequences;
4. 5S rDNA repeat units;
5. 45S rDNA components or repeat units, when appropriate for the analysis.

Do not fabricate these sequences. If `test_data/known_rice_repeats.fa` does not
exist, run only the de novo discovery and keep this as a command template.

## Example command

```bash
tandemx annotate-repeats \
  --catalog results/rice_100k/discover/monomers.fa \
  --known test_data/known_rice_repeats.fa \
  --out results/rice_100k/repeat_annotation.tsv \
  --kmer-size 11
```

The output `repeat_annotation.tsv` reports one best known-repeat match per
discovered family using k-mer Dice, Jaccard, containment and local identity.

## Interpreting the 354/706 bp case

If the 354 bp family matches a published rice tandem repeat such as an
X55642.1-like 354/358 bp repeat, that supports the conclusion that de novo
discovery recovered a real biological signal. It does not mean the known repeat
was used as a template.

If a 706 bp family is related to the 354 bp family, TandemX should report it as
a possible higher-order, dimer-like, partial, or related-family candidate. It
should not be deleted by default. The correct next step is to inspect
`family_similarity.tsv`, `repeat_annotation.tsv`, read support, assembly
localization, and downstream validation evidence.

## Publication-level evidence

Known-repeat annotation is one layer of evidence. For publication-oriented rice
analysis, combine it with:

1. exact discovery command and parameters;
2. read support and copy-number estimates;
3. assembly localization when a compatible assembly is available;
4. comparison against published annotations;
5. FISH probe evidence or experimental cytogenetics when available.

Avoid claiming full array resolution from a 100k read subset. Use language such
as “candidate monomer,” “post hoc known-repeat match,” “possible higher-order
candidate,” and “requires validation.”
