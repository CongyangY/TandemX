# Peanut known tandem/satellite repeat exclusion library — 2026-09-10

## Purpose and boundary

This compact library is a pre-discovery exclusion aid for planned V14167,
K30076 and tetraploid-peanut analyses. It is not a comprehensive repeat
catalogue, an annotation database, or a source of new sequence inference. A
match may identify a known or related repeat for follow-up; a non-match cannot
support novelty.

The sequence-bearing inputs are in
`benchmarks/inputs/peanut_known_repeats_20260910/`. They contain five directly
retrieved GenBank records, with the raw NCBI response, structured FASTA,
unit manifest and SHA-256 receipt. No raw reads or assemblies were downloaded,
and TandemX was not run.

## Resolved direct sequences

| Repeat | Accession | Species | Length | Sequence status | Evidence and intended interpretation |
| --- | --- | --- | ---: | --- | --- |
| H-b-Ah | [KF957858.1](https://www.ncbi.nlm.nih.gov/nuccore/KF957858.1) | *A. hypogaea* | 317 bp | `resolved_direct` | The 2016 FISH study calls this a heterochromatin-band repeat, with mainly A-genome pericentromeric signals; it is not by itself CENH3 evidence. |
| H-b-Ad | [KF957859.1](https://www.ncbi.nlm.nih.gov/nuccore/KF957859.1) | *A. duranensis* | 318 bp | `resolved_direct` | PCR-cloned H-b homolog. |
| H-b-Ai | [KF957860.1](https://www.ncbi.nlm.nih.gov/nuccore/KF957860.1) | *A. ipaensis* | 317 bp | `resolved_direct` | PCR-cloned H-b homolog. |
| B-c-Ad | [KF957861.1](https://www.ncbi.nlm.nih.gov/nuccore/KF957861.1) | *A. duranensis* | 164 bp | `resolved_direct` | B-c homolog; the paper reports a shared 95-bp conserved region among B-c sequences. |
| B-c-Ai | [KF957862.1](https://www.ncbi.nlm.nih.gov/nuccore/KF957862.1) | *A. ipaensis* | 105 bp | `resolved_direct` | B-c homolog; the retrieved GenBank record is 105 bp, whereas the paper describes a 103-bp PCR unit. |

The five records were retrieved directly from NCBI E-utilities as one small
FASTA response. The receipt records its SHA-256 and the derived FASTA SHA-256.
No record was trimmed, rotated, concatenated or converted into a putative
monomer. Consequently the B-c-Ai file length remains the retrieved 105 bp,
rather than being adjusted to the paper's reported 103-bp unit.

## Historical and 2026-paper items not represented by sequence

| Item | Traceable evidence | Sequence status | Handling |
| --- | --- | --- | --- |
| B-c-Ah, historical cultivated-peanut B-genome repeat | Zhang et al. 2012 reports a 115-bp tandem repeat mainly localized to B-chromosome centromeres ([PMID 22797674](https://pubmed.ncbi.nlm.nih.gov/22797674/)); the 2016 study identifies it as B-c-Ah and reports B-c-Ad/Ai accessions. | `accession_unresolved` | No B-c-Ah sequence was guessed from the 95-bp shared block, figure alignment, or B-c-Ad/Ai records. It is listed unresolved in `source_units.tsv`. |
| CentO, 2026 T2T peanut study | Bian et al. report phylogenetic analysis of the CentO repeat unit across V14167, K30076 and A(t)/B(t) subgenomes ([DOI](https://doi.org/10.1038/s41588-026-02577-z); [PMC13175896](https://pmc.ncbi.nlm.nih.gov/articles/PMC13175896/)). | `sequence_not_retrieved` | The accessible article establishes the label and centromere context, but this audit did not locate a direct CentO sequence accession or a supplementary sequence file. No 116-bp sequence is inferred or equated with CentO. |

The 2026 paper contributes the six target materials and centromere context,
but not a sequence entry to this compact library. It states that the two
diploid progenitors and four tetraploid accessions were assembled T2T and uses
CentO repeat-unit monomers in its centromere analysis. That does not establish
the sequence, functional CENH3 binding, or a universal monomer for a future
TandemX catalogue.

## Use in a future discovery run

Use this library only after de novo discovery to label exact or high-confidence
related sequence candidates. Preserve the full discovered catalogue and report
which accession, orientation, alignment scope and threshold generated any
known-repeat label. A library non-match remains `not_compared_comprehensively`
or `sequence_not_retrieved` as appropriate, never `novel`.

The unresolved items should be revisited only when an author-supplied repeat
file, a stable primary accession, or an unambiguous supplementary sequence
resource becomes available. The library currently contains **5 resolved
sequences** and **2 documented unresolved items**.
