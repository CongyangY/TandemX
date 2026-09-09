# Orthogonal abundance source audit v1

Frozen source audit for the final TandemX method-science gate. It identifies
complete public sequencing runs capable of testing the read--assembly abundance
deficit without treating the newer assembly as copy-number truth.

Selected complete runs:

- Ey15-2 PCR-free Illumina: `ERR8666067`;
- Macadamia Illumina: `SRR11191912`;
- Macadamia PromethION: `SRR11191910`.

`ena_selected_runs.tsv` records the exact ENA transfer metadata needed by
`benchmarks/scripts/fetch_ena_complete.py`. `source_audit.tsv` records the
biological relationship and evidence limit. `source_queries.json` preserves the
machine query and retrieval date. No FASTQ is versioned in Git.

`ncbi_sra_transfer_objects.tsv` records official NCBI normalized-SRA and SRA
Lite objects for the same runs. Normalized SRA is preferred when an ENA FASTQ
transfer fails integrity checks. SRA Lite omits original quality scores, so it
is permitted only for the frozen sequence-counting and mapping analyses:
Illumina quality correction is disabled, KMC does not consume qualities, and
the ONT check uses alignment identity rather than base-quality filtering. Any
used SRA object must pass its official MD5, `vdb-validate`, and converted
aggregate record/base-count reconciliation before it can enter the validation.
`transfer_failures.tsv` retains failed transfer evidence; failed objects are
never silently admitted.

The Ey15 paper reports a PCR-free extraction from the same ground tissue used
for HMW DNA. The Macadamia paper describes one DNA extraction from leaves of
clonal tree accession 1005 for its sequencing platforms, while the later HiFi
update is only paper-level same-sample because its archive records differ.
Macadamia's later `SRR13480361/2/3` PromethION objects share experiment
`SRX7812077`; they are excluded to avoid double counting possible alternate
partitions of `SRR11191910`.

Sources:

- Rabanal et al., Nucleic Acids Research (2022), DOI
  `10.1093/nar/gkac1115`, BioProject `PRJEB50694`.
- Nock et al., GigaScience (2020), DOI `10.1093/gigascience/giaa146`,
  BioProject `PRJNA609013`.
- Nock et al., GigaByte (2021), DOI `10.46471/gigabyte.24`, BioProject
  `PRJNA694456`.
- ENA Portal API query URLs and selected rows in this directory.
- NCBI SRA Run Selector metadata URLs and selected Lite-object rows in this
  directory.
