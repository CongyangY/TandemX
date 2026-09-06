# Indel-aware discovery example

Activate `tandemx-dev` and build the current extension with `pip install -e .`.

```bash
tandemx simulate toy --outdir /tmp/tandemx-elastic-input
tandemx discover --reads /tmp/tandemx-elastic-input/reads.fa \
  --outdir /tmp/tandemx-elastic-output --discovery-method elastic \
  --kmer-backend rust --threads 1
tandemx validate --project /tmp/tandemx-elastic-output
pytest tests/unit/test_elastic_discovery.py tests/integration/test_elastic_cli.py
```

Inputs are FASTA/FASTQ reads (optionally gzip). Outputs are `candidate_reads.tsv`,
`candidate_monomers.fa`, `monomer_membership.tsv`, `monomers.fa`, `families.tsv`, `family_similarity.tsv`, `discovery_summary.json`,
`run.log` and `run_config.yaml`. A valid negative result can contain no candidates.
See [field definitions](../../docs/file_formats.md) and
[algorithm](../../docs/algorithms.md#elastic-read-local-discovery-experimental).

The method allows gaps between adjacent repeat copies and several arrays per
read. It does not reconstruct unobserved sequence, establish higher-order repeat
organization, or calibrate copy number. `--kmer-backend python` runs the reference
alignment implementation. `tandemx run` accepts the same discovery-method option.

The independent challenge can be rerun without comparator dependencies:

```bash
python -m benchmarks.challenge.run --config benchmarks/configs/elastic_development_v1.yaml \
  --split development --outdir /tmp/tandemx-elastic-challenge
```

The output directory must be new/empty. Single-run metrics are measured accuracy;
they do not demonstrate runtime variability or deterministic repetitions. Use
separate validation/held-out splits only as specified in the evaluation plan.

By default elastic mode uses sequence clustering at 95% circular edit similarity.
The membership table preserves below-support and ambiguous assignments. Re-run
into a different directory with `--clustering-method legacy` for the historical
clusterer, or `--cluster-identity 0.98` for a finer operational resolution. Compare
these outputs as a sensitivity analysis, not as evidence of separate ancestry.
`run` accepts both options and includes them in resume fingerprints.
