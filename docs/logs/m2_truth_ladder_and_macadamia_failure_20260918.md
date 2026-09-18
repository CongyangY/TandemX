# M2 truth ladder and Macadamia abstention diagnosis, 2026-09-18

## Frozen route and evidence separation

The unchanged M2 alignment prototype is
`benchmarks/m2_routes/alignment/prototype.py`, SHA-256
`043122d1b5987e35dcd957a6471c930f6211fe7a8902055d3a8c37f06d4ed651`.
Its supplied-template global tiling allows per-copy edit distance <=15%,
length variation <=8%, and at most 4,096 array bases. Scores are edit costs,
not posterior probabilities. The five-level ladder protocol was committed at
`4329c5e` **before** Level 1–3 results were generated. Levels 4 and 5 use
earlier frozen real-context runs and were rerun to new temporary directories;
every archived per-case, summary and receipt/log byte matched by SHA-256.
The replay receipt is
`benchmarks/m2_routes/truth_ladder_v1_20260918/results/upper_level_replay.json`.
All five levels remain development diagnostics; they do not constitute
independent biological validation.

| Level | Frozen material and truth | Result | Interpretation |
| --- | --- | --- | --- |
| 1 | Seeded unrelated 96-bp A/B/C/D monomers; exact `ABCABC` and `ABCDABCD` | 2/2 exact paths, resolved | Basic segmentation and label/orientation implementation works on easy rule truth. |
| 2 | Same motifs with two substitutions, one insertion, one deletion per copy; one reversed B | 2/2 exact paths, resolved | Declared small noise and a reversal are handled. |
| 3 | Three archived real Col-CEN operational 178-bp consensus units arranged as synthetic `C1-C2-C3` twice | 1/1 exact path, resolved | A supplied plant-derived catalogue works when an artificial HOR has exact copy boundaries. |
| 4 | Real Col-CEN C3 assembly context; five exact equal-length structural edits, one intact | 2/5 edits detected, 3/5 supported incorrectly; intact supported | Current label/orientation representation misses three known edits; this is one development lineage. |
| 5 | Seven original Macadamia SRA records and nine exact assembly-bp edits | 0/9 technical decisions, 9/9 abstain | No biological accuracy can be computed; all seven read paths lack a full tiling with the frozen template. |

Level 1–3 input sequences, truth labels, recovered labels and hashes are in
`benchmarks/m2_routes/truth_ladder_v1_20260918/results/`. The Level 4/5
original and replay SHA-256 matches are archived with them. This is not an
independent cross-sample performance curve: Levels 1–3 are supplied-template
synthetics; Level 4 is prior Col-CEN development; Level 5 is prior Macadamia
development.

## Exact Macadamia 9/9 abstention path

All nine cases retain the same seven original SRA records, selected before
assembly editing. The full context PAF had 2,550 reported primary alignments:
231 span the entire nominated array, 1,451 overlap it only partially, and
868 do not overlap it. Of the geometric spanners, **seven** distinct records
passed the frozen >=99% whole-alignment identity plus >=1-kb natural-flank
requirements; their seven read arrays are 3,145–3,157 bp. The 1,451 partial
overlaps are not eligible same-locus structural support, and SRA record IDs
do not establish distinct PacBio ZMWs.

Every one of the seven frozen read decompositions was `AMBIGUOUS` with
`no_full_monomer_tiling`. Each nonempty edited assembly case had the same
reason; the complete deletion produced the formally resolved **empty** path,
but it still abstained because all seven read paths were unresolved. No
posterior/confidence probability is defined by this prototype. Per-case and
per-read machine-readable details, including nominal copy count, orientation
availability, HOR candidate scope, exact rule, and partial-spanner counts,
are in `benchmarks/controlled_collapse/macadamia_bp_provisional_v1/m2_diagnosis/`.

The first 144-bp segment of all seven read arrays and eight nonempty edited
arrays has best unrotated-template edit distance **27**, above the frozen
per-copy limit **21**. The left-boundary edit has distance **64**. With a
diagnostic cyclic phase rotation, the first-segment distance falls to 10–12
for the seven original records and to 10 or 11 for the nonempty assembly
cases. Local alignment finds a 3–4-edit template match within the original
arrays. The immediate failure is therefore a **template/array-boundary phase
mismatch** under the current no-rotation global tiler, not absence of any
monomer-like sequence. All unrotated longest tilable prefixes are zero bp.

The extra post-result probe chose one 130-bp rotation from the intact assembly
start and applied it without lowering any edit threshold. It made 6/7 read
arrays and 4/9 edited assembly arrays individually resolvable; one read and
five edited arrays remained unresolved. Because this phase was selected after
seeing the Macadamia source, these are **exploratory diagnostics only** and
must not replace the frozen 0/9 M2 score or become held-out accuracy. The
probe archive is
`benchmarks/controlled_collapse/macadamia_bp_provisional_v1/m2_phase_probe/`.

## Cause classification and scientific decision

| Candidate cause | Evidence-bounded decision |
| --- | --- |
| A. Too few reads | Not the immediate cause: seven qualified records exceed the frozen minimum of three. Molecule independence remains unverified. |
| B. Reads too short | Not the immediate cause: all seven span the selected array and natural flanks. |
| C. Monomer decomposition failure | **Observed primary failure:** no full path; first unrotated tile fails the fixed emission limit. |
| D. Family non-identifiability | The supplied Macadamia catalogue has only one 144-bp family template, so multi-label HOR order is not inferable from this input. No claim of biological absence of HOR follows. |
| E. Path ambiguity | No competing complete path was reached; the failure occurs before a resolved-path comparison. |
| F. Confidence threshold too strict | The frozen 15% emission limit excludes the unrotated start. Phase rotation lowers the distance without relaxing it, so merely reducing a confidence threshold is unsupported. |
| G. Intrinsically absent assembly/read signal | Not established. Local template matches and phase-probe recovery show some sequence information, but full structure remains unresolved. |
| H. Implementation defect | No simple implementation failure appears on Levels 1–3. Missing phase invariance at arbitrary real-array boundaries is a concrete **method representation defect**. |

The present full-HOR reconstruction route is **NO-GO as a TandemX core
contribution**. The frozen structural Gate B remains NO-GO, and no post-result
parameter change has been used to reverse it. A narrower research target is
read-supported structural **discordance detection** with explicit abstention,
not full HOR reconstruction or calibrated copy-count reconstruction. Even that
narrower target has not passed independent real-material validation or
competitor-matched evaluation and is not promoted to production. Any
phase-aware algorithm would be a new version requiring a new held-out source
cluster; Ey15-2, Col-CEN and Macadamia remain development materials.
