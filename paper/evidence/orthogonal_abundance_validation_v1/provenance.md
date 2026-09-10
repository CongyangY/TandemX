# Provenance and execution boundary

The interpretation contract was frozen in
`benchmarks/configs/orthogonal_abundance_validation_v1.json` and
`docs/orthogonal_abundance_validation.md` before the orthogonal family results
were inspected. Finalization used
`benchmarks/scripts/finalize_orthogonal_abundance_validation.py` without
changing families or thresholds.

Key executables were KMC/KMC tools 3.2.4, minimap2 2.31-r1302 and NCBI SRA
Toolkit 3.4.1. KMC sequence counts were run at k=21 and k=31 with minimum count
2 and counter maximum 1,000,000,000. The final target extraction used the read
database as the left operand with `-ocleft`. The Macadamia ONT alignment used
the frozen competitive `map-ont -c --eqx --secondary=yes -N 5` rule; the HiFi
calibration used the corresponding `map-hifi` rule.

Large evidence retained outside Git:

- raw and converted reads: `/Volumes/T7/Codex/TandemX/data/reads`;
- empirical controls, targets and KMC databases:
  `/Volumes/T7/Codex/TandemX/validation/orthogonal_abundance_v1`;
- Macadamia HiFi/ONT PAF files: the `ont` subdirectory of that validation root.

The compact outputs were copied from that validation root after finalization.
The PAF SHA-256 is retained in `ont/ont_occupancy_summary.json`; input-object
MD5/SHA-256 values are retained in `qc/`. Failed inputs and invalid intermediate
counts are listed in `run_fates.tsv` and were not used in final interpretation.

The repository revision immediately preceding this archive was
`d5d5edca370f8fe038ad5a3088f1cb4c5cd9bef2`. The final archive commit is
recorded by Git history rather than recursively embedded into the archive.
