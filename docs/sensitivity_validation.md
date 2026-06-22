# TandemX Sensitivity Validation

## Purpose and boundary

Sensitivity validation asks whether a repeat family that should be detectable is recovered by the de novo workflow. Known repeats are interpretation references, not discovery templates. The required order is:

```text
reads -> tandemx discover -> monomers.fa -> post hoc known-repeat comparison
```

`tandemx discover` must not read the known-repeat FASTA, truth files, or spike-in labels. This separation prevents reference-assisted recovery from being reported as de novo sensitivity.

## Post hoc known-repeat matching

Run:

```bash
tandemx annotate-repeats \
  --catalog results/discover/monomers.fa \
  --known validation/known_repeats.fa \
  --out results/repeat_annotation.tsv \
  --kmer-size 11
```

For each discovered family, the annotation command compares exact k-mer sets against every known repeat in forward and reverse-complement orientations. It reports Dice, Jaccard, containment, best local identity, orientation, and a thresholded post hoc annotation. This method is lightweight and useful for triage, but it is not a replacement for a gapped alignment, repeat-family phylogeny, cytogenetic validation, or manual biological review.

Known-repeat annotation must not be used to tune `tandemx discover` until a separate validation experiment has been completed. For publication-level analysis, report the de novo discovery command, then report annotation as downstream interpretation.

## Synthetic spike-ins in a real-read background

A publication-oriented sensitivity benchmark should preserve a fixed real HiFi background while adding simulated reads containing controlled tandem arrays. The benchmark generator should write spike-in truth separately from analysis inputs and record:

1. background dataset identifier and checksum;
2. spike-in monomer sequence and length;
3. array copies per spike-in read;
4. number and fraction of spike-in reads;
5. simulated error model and random seed;
6. TandemX parameters and software version.

Discovery receives only the combined reads. After completion, truth monomers are compared to `monomers.fa` with the post hoc checker. Background-only negative controls are required to estimate false matches.

## Divergence experiment

Generate monomer variants at several controlled divergence levels, for example 0%, 2%, 5%, 10%, 15%, and 20% substitutions. Use multiple random seeds and retain sequence complexity. Report:

1. recovery at each divergence level;
2. monomer-length error;
3. best post hoc similarity;
4. family splitting or merging;
5. false-positive families in matched background controls.

Insertions, deletions, truncated units, higher-order organization, and mixtures of related families should be separate experiments because a substitution-only curve does not represent those structures.

## Read support and copy-number experiment

Cross spike-in read support with within-read array copy count. A minimal grid should vary:

1. supporting reads below, at, and above `--min-support-reads`;
2. repeat span below, at, and above `--min-repeat-span`;
3. short and long monomers across the configured period range;
4. sequencing depth or sampling fraction;
5. error rate and read length.

Report family recall, candidate-read recall, monomer-length error, similarity score, and runtime for each cell. Threshold-edge cases should be labelled expected non-recovery when they do not meet the configured detection boundary.

## Publication-level validation path

The framework supports a defensible evidence ladder:

1. deterministic toy tests establish implementation correctness;
2. spike-ins quantify sensitivity and threshold behavior in realistic sequence backgrounds;
3. post hoc comparisons test recovery of published repeat families without guiding discovery;
4. matched assemblies assess possible under-representation;
5. independent annotations and experimental FISH results assess biological relevance.

Results must distinguish engineering validation from biological validation. A high k-mer match supports family correspondence but does not prove array location, exact copy number, chromosome specificity, or successful FISH hybridization.

## Optional family collapse during sensitivity analysis

`tandemx discover --collapse-redundant-families` is disabled by default. When
enabled, it collapses only `likely_redundant` pairs from `family_similarity.tsv`
and writes an audit trail in `family_collapse.tsv`. Relationships labelled
`possible_higher_order_or_partial` are retained because they may represent
higher-order organization, dimers, partial overlaps, or genuinely related
families. Publication analyses should report both the original and optional
collapsed catalog decisions when collapse mode is used.
