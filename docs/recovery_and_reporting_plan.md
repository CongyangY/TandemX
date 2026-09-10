# Targeted recovery proof of concept and user workflow

Authorized 2026-09-10 by the user's attached development request and explicit
follow-up. Current narrative: `paper/0910/manuscript_v2.md`. This bounded mandate
supersedes the earlier feature freeze only for recovery feasibility and the
requested user interface/reporting work. Detector, abundance, localization and
comparison thresholds remain frozen; unitFinder remains permanently stopped.

## Decisions and division of work

1. Astra owns scientific design, recovery, evidence interpretation, integration,
   review and manuscript consistency. Terra implements the simplified run
   interface and independently the family data/HTML report. Luna implements
   static SVG/PDF/PNG exports. Subtasks cannot change scientific thresholds.
2. Extend the existing pipeline, sequence readers and family audit. Preserve
   original stage products. Add a summary-facing `families/`, `figures/`,
   `summary.tsv`, `summary.json` and `report.html` rather than replacing evidence.
3. `run --reads ... [--assembly ...] -o ...` uses existing numerical defaults,
   strict optional YAML configuration and fingerprint-validated automatic resume.
   Assembly length may supply a clearly labelled provisional normalization
   denominator; it is not an independently estimated haploid genome size.
   Without assembly, genome size or depth, discovery remains available and
   absolute abundance stays unavailable, with an actionable explanation.
4. Architecture views preserve actual sequence-supported pair evidence. Period
   multiples do not establish unit order or a biological HOR. Missing values,
   failed stages and absent analysis remain distinguishable from zero results.

## Recovery enrollment fixed before new locus/read inspection

The selected families are Ey15-2 TXF000002 and TXF000154, the two families with
stable orthogonal support in the completed independent-abundance archive. No
third family is promoted from unresolved. The previously analysed family-level
newer-assembly quantities are already known: this is recovery-input separation,
not a claim of a wholly unobserved validation set.

Inputs are the enrolled historical CLR-Canu assembly, frozen discovery catalogue,
old-assembly localization, original full HiFi ERR8666125 and its existing QC.
No newer-assembly sequence, newer locus coordinates, alignments to the newer
assembly or newer-based single-copy controls may guide recruitment or recovery.
Record input hashes, exact commands, software hashes and all failure outcomes.

Use every existing historical localized interval for the two selected families.
A family with no localized historical interval retains
`unresolved_no_assembly_locus`; do not manufacture a locus from newer coordinates.
For each interval examine 2-kb flanks at offsets 0, 2, 5 and 10 kb on both sides.
Incomplete, ambiguous or out-of-contig flanks are retained as rejected anchors.
Map flank sequences against the historical assembly with pinned minimap2,
retaining secondary hits. An eligible anchor must have >=90% query coverage,
>=98% identity, its expected historical hit and no competing >=90% coverage /
>=95% identity hit. Uniqueness is relative to this assembly, which may itself
omit copies, and must be checked again for ambiguous read placements.

Recruit reads competitively to the selected tandemized representatives and
candidate flank sequences with minimap2. Preserve repeat-only, single-flank,
dual-flank and ambiguous assignments. No family-wide abundance deficit is
allocated to a particular locus. Require matching strand, order and distinct
read IDs for dual-flank support. A direct spanning-read extraction is the first
lightweight candidate reconstruction; it does not require writing an assembler.
At least three independent dual-anchor reads and concordant span lengths
(max[100 bp, 1% median length]) are required before selecting the median-span
observed sequence as a supported candidate. Retain conflicting paths rather
than averaging repeat copy counts. An unpolished read extraction cannot be
labelled a fully resolved assembly; use `partially_resolved` and explain the
sequence-accuracy/haplotype limits. `resolved` remains reserved for independently
supported complete reconstruction and validation, not endpoint availability.

If no pair of eligible anchors exists, stop that locus as
`unresolved_no_unique_anchor`. If reads anchor only one side and do not span,
report the observed read-span bound and
`unresolved_array_exceeds_read_information` only when supporting evidence
justifies that bound; otherwise use `insufficient_read_support`. Never infer
an array's physical length solely from family-level k-mer abundance.

Lock candidate FASTA/table hashes and a completion receipt before allowing the
newer assembly into a separate post-recovery evaluation. Compare matching
flank-supported intervals, sequence representation and length, describing the
newer assembly as a reference proxy. If all candidates lack reliable anchors
or sufficient read information, retain a negative PoC and stop algorithm
expansion. Formal module promotion requires at least one clear sequence gain
with independently supported placement and improved post-hoc proxy agreement.

## Required outputs and verification gates

Recovery writes candidate and validation TSVs, 0-based BED, recruited-read
evidence, candidate-only FASTA and an offline report. Unresolved runs have
header-only or empty sequence outputs with explicit statuses, never synthetic
recovered sequence. Original assemblies are never patched automatically.

Each stage follows implementation -> focused tests -> toy/real execution ->
visual/source review -> full regression -> reviewed commit. Large data, mapping
outputs and retained failed attempts stay under `/Volumes/T7/Codex/TandemX`.
Only selected compact evidence is versioned. Figure QA checks visible layout,
source values, editable SVG text, PDF rendering and offline HTML behavior.
Report software delivery and manuscript readiness separately.

## Implementation audit and bounded continuation

The primary anchor pair is the nearest eligible flank on each side, selected
using the historical assembly alone. More distant offsets contribute to the
uniqueness audit but are not selected according to read outcomes. This is a
bounded primary-pair experiment, not an exhaustive search over all possible
anchors. `maximum_read_span_bp` is the maximum observed dual-anchor interval;
the distinct `maximum_recruited_read_length_bp` is the longest full recruited
read. Candidate span includes intervening non-repeat sequence when offsets are
nonzero, so `candidate_span_bp` is separate from `recovered_bp`, which remains
unavailable until repeat-specific validation.

The first preparation-only v1 is retained because its implementation did not
yet seal the prepared state. v2 added complete preparation locks without changing
flank or read thresholds. The v2 long-template mapping was operator-stopped for
the bounded resource budget after its first 2,173-read batch required 119.606 s;
this is not an OOM, natural tool failure, formal timeout or biological negative.
Its partial PAF and command receipt remain on T7.

The separately frozen v3 changes only the repeat recruitment template minimum
from 30,000 to 1,000 bp. These short tandemized sequences test family presence,
not array occupancy or length. All identity, aligned-span, anchor and candidate
support rules remain unchanged. It imports the locked v2 flank audit without
rerunning that mapping, verifies the full source-read hash again and first runs
a non-inferential 1,000-read cost pilot. No recovery candidate or newer sequence
was inspected to choose the shorter template. The result must retain both the
failed execution path and the limited short-template recruitment semantics.


The v3 cost pilot completed in 1.146 s for the first 1,000 reads and retained
5,984 qualifying alignment rows from 190 distinct reads. This pilot is only a
resource check. Before full-read candidate outcomes were available, v4 changed
the collector to stream evidence rows and raised only the explicit resource
caps to 1,000,000 distinct reads and 10,000,000 rows. It imports the identical
completed v3 mapping after verifying commands, target/input hashes and mapper
identity; no biological threshold or anchor selection rule changes. If the
original v3 collector hits its smaller cap, retain that technical failure
separately from the v4 biological outcome.
