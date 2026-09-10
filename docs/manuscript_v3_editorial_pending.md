# Manuscript v3 editorial pending items

This checklist is intentionally outside the manuscript narrative. It records
package work that remains to be verified; it does not claim that the final
figures, release, Bioconda recipe or Zenodo deposition exist.

## Figures and source data

- Build and inspect the editable-vector workflow schematic and its source table.
- Finalize the controlled benchmark composite, old/new assembly comparison,
  orthogonal validation and TXF000695 boundary figure.
- Confirm every panel has a matching legend, source-data table and provenance.
- Keep accepted and rejected layouts archived, but do not describe rejected or
  planned layouts as manuscript results.
- Complete figure numbering and supplementary migration only after the final
  artwork set is frozen.

## Recovery supplement

The full bounded recovery Results and Methods are in
`paper/0910/supplementary_recovery_v3.md`. Preserve the 5,681-bp span, 644 versus
679 repeat-bp comparison, 97/120 shared-span support and the zero validated gain.

## Release and metadata

- Verify reference 12 and dataset-specific metadata in a reference manager.
- Complete author, funding and competing-interest fields from authoritative
  records; do not infer missing values.
- Track software release, Bioconda and Zenodo work here until actual identifiers
  exist; do not put pending status in the manuscript as completed availability.
- Keep commit/CI IDs, hashes, run fates, resource profiles and migration history
  in the technical/evidence package or supplementary materials.

## Threshold sensitivity

The manuscript audit found controlled tables with decision thresholds 0.5 and
0.6, but did not verify a final catalogue-level 0.6 sensitivity table for the
old/new and orthogonal analyses. This remains an evidence-status item and must
not be described as threshold robustness without the corresponding frozen table.

## Archived v2 construction and migration text


### Figure 1. Biological problem and TandemX workflow

**Biological problem and TandemX workflow.** **A,** A long, homogeneous plant
satellite array can be present but under-represented in an otherwise highly
contiguous assembly. **B,** HiFi reads are scanned for tandem arrays and
sequence-supported monomers are organized into an operational family catalogue.
**C,** Diagnostic-k-mer depth estimates read-derived family abundance, whereas
the same representative is localized independently in an assembly. **D,** The
comparison produces represented, possible under-representation and unresolved
states; the continuous read--assembly abundance deficit is not physical missing
bp. **E,** HiFi data generate the TandemX prediction, after which independent
Illumina k-mer abundance and ONT direction can provide orthogonal support. The
latter are validation measurements rather than additional TandemX input modes.
This figure requires a new editable-vector schematic and a compact definitions
source table.

### Figure 2. Integrated analytical validation across controlled benchmarks

**Integrated analytical validation across controlled benchmarks.** **A,** De
novo read-level array and family recovery across repeat periods, divergence and
read-error conditions. **B,** Abundance error across coverage and error strata,
including both improved and worsened family conditions. **C,** Assembly
localization recall and precision, with the divergent and interrupted-array
boundary shown explicitly. **D,** Binary under-representation sensitivity,
false-positive rate and precision, displayed separately from continuous
abundance error. **E,** Task-matched read comparisons with TideHunter and TRF
and assembly comparisons with TRASH/TRASH2 and TideCluster. **F,** Summary of
accuracy, runtime and memory trade-offs. The figure combines only the final
validation and principal adverse boundary; development iterations remain in
Supplementary Figures. Source panels are assembled from the current controlled-
benchmark evidence tables without selecting only favourable strata.

### Figure 3. Repeat-family landscapes across diverse plant long-read datasets

**Repeat-family landscapes across diverse plant long-read datasets.** **A,**
Species, material, genome-size denominator and validated sequence yield for ten
libraries from eight reported plant species. **B,** Number and abundance
distribution of operational repeat families at matched nested sampling scales.
**C,** Monomer-length and sequence-complexity landscape. **D,** Descriptive
fraction of read-derived family abundance represented in available assemblies.
**E,** Representative family profiles from small and large plant genomes.
**F,** Provenance and denominator map identifying pooled material, technical
batches and unresolved donor/reference matches. Called-base fraction is not an
accuracy endpoint, and the libraries are not treated as interchangeable
biological replicates. A consolidated source table must be generated from the
existing cohort-QC and real-diagnostic tables before this figure is rendered.

### Figure 4. Old/new assemblies provide reference-proxy evidence for historical under-representation

**Old/new assemblies provide reference-proxy evidence for historical
under-representation.** **A–B,** Ey15-2 and Macadamia family representation in
historical and newer assemblies, highlighting families independently predicted
from HiFi reads to be depleted in the historical assembly. **C,** HiFi-derived
abundance versus historical-assembly representation. **D,** New-versus-old
assembly gain compared with HiFi-versus-old abundance deficit; binary agreement
and magnitude disagreement are displayed separately. **E,** Sensitivity to the
5-, 15- and 50-kb newer-assembly eligibility denominators. **F,** Assembly-to-
assembly aligned-query coverage as explanatory context. The newer assemblies
are high-quality reference proxies, not absolute or independent copy-number
truth. The accepted Ey15-2 panels and Macadamia family tables will be rebuilt as
one two-species composite rather than retaining the old Ey15-only main figure.

### Figure 5. Orthogonal reads support residual under-representation in three repeat families

**Orthogonal reads support residual under-representation in three repeat
families.** **A,** Frozen HiFi abundance and newer-assembly representation for
six candidates selected before orthogonal inspection. **B,** Independent
Illumina k=21 and k=31 estimates. **C,** Direction-level Macadamia ONT support
across the 616–780-Mb genome-size sensitivity. **D,** Stable support for Ey15-2
TXF000002 and TXF000154. **E,** Cross-platform support for Macadamia TXF000496.
**F,** Final evidence matrix showing three supported and three unresolved
families. The split is not an accuracy rate, and all continuous differences are
read--assembly abundance deficits rather than physical missing-base truth. The
current accepted orthogonal-validation figure supplies the source data but will
be expanded to six focused panels.

### Figure 6. Discordant families and the confidence boundary of read-derived abundance

**Discordant families define the confidence boundary of read-derived
abundance.** **A,** Cross-k and cross-platform outcomes for TXF001517,
TXF000563 and TXF000695. **B,** Complete TXF000695 diagnostic-k-mer depth
distributions at k=21 and k=31, showing the high- and low-depth components and
the shift of the median between them. **C,** Diagnostic-set size and observed
fraction: both k values contain 467 diagnostic words, of which 381 and 299 were
observed at least twice. **D,** Median depth differs by approximately 763-fold
and estimated abundance by 670-fold. **E,** ONT and HiFi mapping occupancy are
closer to the low-abundance direction, while the mechanism remains unresolved.
**F,** Reporting-level QC summary showing cross-k fold, observed fraction,
MAD/median and the heuristic labels `unstable_across_k`,
`depth_median_boundary_sensitive` and `reliability=low`. These warnings are not
a validated classifier and do not alter the estimator or family interpretation.
The derived QC values and their raw-input mapping are retained in
`paper/0910/source_data/figure6_txf000695_qc.tsv` and its accompanying notes.

## Supplementary migration map

The technical master, all source data and every adverse or incomplete result are
retained. Commit identifiers, CI run identifiers, exact seeds, hashes,
byte-identity checks, failed transfers, operator-stopped runs, build histories,
development versions and evidence E1–E42 are moved out of the narrative Results
and remain available through `paper/manuscript.md`, Supplementary Methods,
Supplementary Results, Source Data and `paper/evidence/`. No technical record is
deleted or relabelled as successful evidence.

| Technical-master material | Destination in the new article package |
| --- | --- |
| Original Figures 1–2 (multi-k and joint-read uncertainty) | Condensed in new Figure 2; full panels and all unavailable intervals retained as Supplementary Figures S1–S2 |
| Original Figures 3–7 (conditional comparison, domain shift, localizer and classifier sequence) | Final validation and principal boundary condensed in new Figure 2; complete development, failed held-out and final depth-gated results retained as Supplementary Figures S3–S7 |
| Original Figures 8–11 (cascade, calibration and performance gates) | Essential task-matched trade-offs condensed in new Figure 2; complete gates, runtime distributions and adverse cells retained as Supplementary Figures S8–S11 |
| Original Figure 12 (Ey15-2 proxy analysis) and Macadamia proxy evidence | Rebuilt together as new Figure 4; full species-specific panels retained in Supplementary Figures |
| Original Figure 13 (orthogonal validation) | Source for new Figure 5; accepted and rejected layouts remain archived |
| New TXF000695 depth-distribution analysis | New Figure 6; exact diagnostic-word table and summary become Source Data |
| Current Supplementary Figures S1–S6 | Renumbered after the above transfer without dropping input-QC, real-diagnostic, TideCluster or failed-run panels |
| Original Supplementary Tables S1–S39 and Evidence E1–E42 | Retained in full and renumbered only after main-text figures and tables are finalized |

The main text will cite the biological or methodological conclusion; the
supplement will carry the exact frozen configuration, complete unfavourable
rows, run fate, reproduction record and resource context. This separation
changes article hierarchy, not evidence content.
