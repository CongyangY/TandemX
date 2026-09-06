# Observed reference concordance and compartment context

Whole-library seed6101 subsets of ERR6210723 were aligned to the independently
checked Col-CEN v1.2 reference, including ChrM and ChrC. Native minimap2
2.31-r1302 used fixed map-hifi/base-alignment settings, one thread,50-Mb query
batches and one3-Gb index batch. Native hashes, source snapshots and full-query
denominators are retained. Source was an explicitly hashed development snapshot
on top of d2aa4d6, not claimed to be that commit's exact implementation.

The11.766-Mb control passed:751/752 reads have primary alignments,182 reads have
at least one organellar primary alignment. The118.497-Mb follow-up passed with
7,557 reads,9,712 native alignment rows (7,800 primary and1,912 secondary).
7,554/7,557 reads map (99.9603%);7,089 have primary MAPQ20–254 alignments.
Primary query-span union covers99.9396% of input bases. Three unmapped reads
remain in the denominator; they are not classified as contaminants.

Primary query spans on ChrM/ChrC total26,290,920 bp, **22.1871% of all input
bases**. Primary spans on the other five contigs total92,139,003 bp;4,936 bp
overlap these compartment-specific spans. Thus these categories are not blindly
summed as disjoint sequence origins.1,693 reads have organellar primary
alignments. Most organellar aligned target bases map to ChrC, whose summed
primary alignment-base depth is166.672×, versus0.671–0.766× on the five
non-organelle contigs. ChrM depth is1.423×. Native secondary output is capped;
MAPQ or primary status is not proof of unique genomic origin.

This observation makes total-library normalization a concrete concern: dividing
all118,496,530 bases by the five-chromosome131,559,676-bp size gives0.9007×,
although a substantial fraction aligns to organellar reference sequence.
That is not yet a measured22.2% copy-number bias or a validated correction.
Reference divergence/collapse, integrated organellar sequences, multimapping,
read lengths and empirical nuclear single-copy depth must be considered before
changing the denominator of a repeat-copy estimator. The current research
quantifier has not automatically removed organellar reads.

The118-Mb native child took48.840 s and1,291.563 MiB peak RSS; concurrent
diagnostic timings are not comparator rankings. SQLite-backed QC validates
query/target lengths, all PAF coordinates/types and CIGAR consumption, excludes
deletions from aligned reference blocks, unions query spans, and keeps missing
MAPQ255 distinct from high confidence. Source/helper snapshots for both runs
are archived because the compartment-span fields were added between controls.
PAF and SQLite stay at the exact T7 paths in the receipts; output hashes enable
verification. File/mapping QC does not establish family recall, FISH success,
assembly completeness or contamination status.

The Nipponbare 11.418-Mb control used frozen438d078 and the exact checked
GCA_034140825.1 reference. All626 reads had primary alignments;611 had primary
MAPQ20–254 alignments. Primary span union was11,408,400 of11,418,016 bp
(99.9158%). The native output contained639 primary and66 secondary rows.
No organellar contigs are present in this GCA assembly, so its zero organellar
counts cannot establish absence of organellar reads. The matching cultivar/
project does not resolve different raw/reference BioSample dates and ages.
This control supplies reference concordance only, not validated nuclear depth.
