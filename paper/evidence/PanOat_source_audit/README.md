# PanOat primary-source audit

Source: Avni et al., Nature649,131–139 (2026),
https://doi.org/10.1038/s41586-025-09676-7. Original official supplements and
peer-review PDF are retained under the T7 path in retrieval_receipt.json;
SHA-256 and exact URLs are recorded. The original workbook was not edited.
`selected_source_cells.json` preserves cell addresses and native values from
Table S1/S20. The latter worksheet has an old "Table S17" title in cell A1;
we retain its actual worksheet name rather than silently correcting the source.

Victoria is explicitly PacBio HiFi + Hi-C in Table S1!D12, AACCDD in C12,
with10,732,944,064-bp assembly in E12 and10,608,881,087-bp pseudomolecules in H12.
Table S20!F27 matches the archived raw-read sample SAMEA111508775; Q27/R27 link
PRJEB56706 /GCA_947311595. S20!N27:O27 are Iso-Seq read/yield columns, not genomic
HiFi coverage. Do not use those counts as the genomic sequencing denominator.

Public peer-review lessons (PDF pages3–4,23–24), paraphrased: describe sample
sizes, population structure and model diagnostics; connect biological examples
directly to the generated resource; test alternative quantification when related
subgenomes cause assignment errors; supply underlying phenotype/metadata and
explain supplementary columns. An exploratory small-cohort association must
remain exploratory. These are study-specific critiques, not numerical journal
requirements for TandemX. For TandemX they motivate family ambiguity tests,
independent-material confidence intervals, explicit truth/data schemas and a
coherent biological demonstration rather than unrelated panels.
