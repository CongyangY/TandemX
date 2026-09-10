# Supplemental soybean and wild-soybean repeat-source audit

**Date:** 2026-09-10
**Scope:** public, sequence-resolved tandem/satellite/centromeric-repeat records relevant to *Glycine max* and *Glycine soja*. This is a bounded supplement to [the existing 13-clone audit](soybean_known_repeat_sequence_sources_20260910.md). It does not query, upload, or compare TXF candidates or any other unpublished TandemX sequence against external services.

## Result

The existing local 13-clone set is incomplete at the **accession/source** level. This supplement preserves seven additional directly obtainable published sequences:

- four historical *G. max* GenBank records: `Z26334.1`, `AF297983.1`, `AF297984.1`, and `AF297985.1`;
- three exact monomer strings published in Supplementary Table S2 for SBRS1–3. The study reports signal evidence in both cultivated soybean and wild soybean (*G. soja*).

They are stored in [`benchmarks/inputs/soybean_known_repeats_supplement_v1_20260910/`](../benchmarks/inputs/soybean_known_repeats_supplement_v1_20260910/). They are **post hoc known-repeat exclusion references**, not TandemX discoveries, abundance estimates, or a new family catalogue.

## Evidence and retrieval

| Records | Evidence for being in scope | Exact source | Sequence handling |
|---|---|---|---|
| `Z26334.1` | GenBank annotation: *G. max* satellite DNA | [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/Z26334.1) | Retrieved directly from NCBI E-utilities; original bytes retained in `genbank_historical_tandem_records.fasta`. |
| `AF297983.1`–`AF297985.1` | GenBank annotations: *G. max* TRS1–3 tandem repetitive-repeat regions | [AF297983.1](https://www.ncbi.nlm.nih.gov/nuccore/AF297983.1), [AF297984.1](https://www.ncbi.nlm.nih.gov/nuccore/AF297984.1), [AF297985.1](https://www.ncbi.nlm.nih.gov/nuccore/AF297985.1) | Retrieved directly from NCBI E-utilities; original bytes retained. |
| `SBRS1`–`SBRS3` | Chen et al. report tandem repeats with lengths 48, 124, and 201 bp; FISH/Southern evidence covers *G. max* and *G. soja*. | [Open-access article and Supplementary Table S2](https://pmc.ncbi.nlm.nih.gov/articles/PMC5430363/) | Exact strings transcribed only from the authors’ Table S2 and labelled `direct_published_supplement`. The table strings have lengths 92, 124, and 196 bp respectively; the discrepancy for SBRS1/SBRS3 is retained rather than edited. |

The historical-record search is additionally consistent with the primary report that characterized a soybean satellite family ([PubMed 8900841](https://pubmed.ncbi.nlm.nih.gov/8900841/)). This audit uses direct sequence records rather than relying on the patent's family labels.

## What this does and does not establish

The current local 13-clone set does not contain any of the seven accessions/identifiers above. An exact forward/reverse-complement full-string containment check against that set found no full query record contained in an existing clone. This verifies an accession-level omission; it does **not** establish seven independent repeat families, homology, centromeric localization, or copy-number relevance. Tandem units can differ by rotation, truncation, divergence, or clone context, so any family merge needs an explicitly recorded sequence analysis under the relevant study protocol.

The previous set already contains `U11026.1` (SB92) plus `CentGm-1`/`CentGm-4` records. This supplement does not relabel any historical TRS or SBRS record as SB92/CentGm based on a name, presumed provenance, or a partial pattern.

## Unresolved boundaries

- `CentGm-2`, `CentGm273`, and `CentGm444` remain literature names without an independently located, direct public sequence record in this audit. They are not represented by proxy sequences.
- The YSD56 assembly-derived `trf*` annotations remain outside this known-repeat source set; no assembly sequence was extracted here.
- `AB536713.1` remains excluded as a mixed retrotransposon/CentGm-4 fragment under the prior audit.
- No assertion about a candidate's novelty, homology, or exclusion outcome follows from this inventory alone. Such calls require a locally recorded sequence comparison and the bounded known-library exclusion protocol.

## Files and integrity

| File | Purpose | SHA-256 |
|---|---|---|
| `genbank_historical_tandem_records.fasta` | Original four-record NCBI retrieval | `963c4c93443e0ab798349602acc7a74ca9e060f28f55536785fad10c7b81e9ce` |
| `soybean_wild_soybean_published_repeats.fasta` | Seven normalized headers with original sequence strings | `9c47db01bede8adb8043e549b5d901d66bf62bb32237c1c6a4f282fe0051b784` |
| `source_inventory.tsv` | Field-level provenance, availability, and unresolved states | generated locally from the sources above |

No raw reads, assemblies, or unpublished candidate sequences were downloaded, uploaded, or sent to a remote service.
