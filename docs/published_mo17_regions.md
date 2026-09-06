# Published Mo17 region evidence

Source: Chen et al., Nature Genetics 55, 1221–1231 (2023),
[doi:10.1038/s41588-023-01419-6](https://www.nature.com/articles/s41588-023-01419-6).
The CC-BY-4.0 article's original Supplementary Tables workbook is preserved in
`paper/evidence/Mo17_published_regions`, SHA256
`863da22bbc8fa3170035f966f6ced723e418e81541fda96de0e5cc9355569e1d`.
It is read without modifying the workbook. Independent artifact-tool and
standard-library XML readers agreed on all 2,033 cells in the selected sheets.

The extraction contains 64 satellite regions (TR-1 17, knob180 25, CentC 17,
and one each of Cent4/tRNAsat/sat112/sat261/sat268), two rDNA regions, ten
CENH3-defined centromeres and 20 telomeres. These are published annotations and
assay-defined regions, not 96 independently measured repeat families.

Every row has a source sheet, row, coordinate-cell range, component-percentage
cell, source URL and SHA256. Preserve the component percentages: mixed regions
can contain extensive TEs and other repeats. Region occupancy cannot be treated
as all-satellite base truth or as an exact copy-number measurement. CENH3-defined
centromeres are a different endpoint from sequence satellite arrays.

The exact published genome size (2,178,604,320 bp) and all ten long-arm endpoint
coordinates match the NCBI GCA_022117705.1 chromosome report. Short-arm telomeres
start at zero; all selected intervals have `size = end - start`. This supports
an inferred zero-based half-open convention, but the table does not explicitly
label its origin. Raw coordinates remain unchanged with an inference warning.
Table 13 has mixed size arithmetic and is excluded from normalization. Numeric
compatibility does not prove byte identity to another assembly release or exact
donor/extraction identity. A one-base boundary sensitivity is required if these
coordinates later support fine-scale endpoints.

Reproduce without an Excel package:

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.curate_mo17_regions \
  --source-xlsx paper/evidence/Mo17_published_regions/41588_2023_1419_MOESM4_ESM.xlsx \
  --ncbi-report paper/evidence/Mo17_input_qc/reference/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_assembly_report.txt \
  --outdir /path/to/new/regions
pytest -q tests/unit/test_published_region_curation.py
```

`published_regions.tsv` fields: `region_id` (source-derived stable ID),
`region_type`, `published_label`, `chromosome_label`, `reference_contig`,
`reference_length_bp`, `start_raw`, `end_raw`, `published_size_bp`,
`coordinate_note`, `annotated_component`, `component_percent` (0–100),
`source_sheet`, `source_row`, `coordinate_cells`, `component_cell`,
`source_url`, `source_sha256`, and `warning`. The TSV has a header and is not BED.
`curation_receipt.json` records source/report/script/output hashes, type counts,
chromosome compatibility checks, excluded table and inference limits. Changed
source hashes, invalid bounds/sizes, unsupported reference sets or formulas in
selected sheets fail explicitly. Table 14's unrelated expression formulas are
not read by this region extractor.
