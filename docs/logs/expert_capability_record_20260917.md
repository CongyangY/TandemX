# TandemX expert capability record, 2026-09-17

This record describes the checked repository and archived evidence as of the
current development branch. It distinguishes implemented software, completed
technical experiments, independent validation, and open claims. `README.md`
describes usage; `docs/current_status.md` contains the long chronological
handoff. Historical measurements are indexed in `performance_scope_20260917.md`.
This document is an evidence map, not a publication-ready Methods claim.

## Implemented software and operational maturity

The repository is Python package version **0.1.0** with optional Rust native
components and a dedicated `tandemx-dev` conda environment. Public CLI
commands include `run`, `discover`, `quantify`, `locate`, `compare`, `probe`,
`visualize`, `cohort`, `simulate`, `import`, `validate`, and post hoc
`annotate-repeats`. The core workflow starts discovery directly from reads;
the main catalogue is not seeded with a truth monomer library. With an
assembly, subsequent stages estimate read-derived family evidence, localize
array candidates, compare read and assembly abundance, rank probe candidates,
and produce an offline HTML report and static SVG/PDF/PNG figures. Each stage's
primary tables, run configuration, and logs remain inspectable. Output schemas
and warnings are documented in `docs/file_formats.md`.

`tandemx run` reuses only validated completed stages whose input and effective
command fingerprints match. Direct `discover` and `quantify` have optional
in-stage checkpoints with changed-input and corrupt-checkpoint refusal.
Discover checkpointing currently requires one input file, full sampling, and
no automatic discovery budget; it snapshots all accumulated candidates and
rehashes/replays on resume, which can add substantial I/O. The orchestrated
`run` command does not expose those in-stage options. Interrupted partial
outputs are not eligible for validated reuse. These controls reduce silent
continuation risk but do not prove fault tolerance on every multi-hour,
multi-file plant-genome run.

The code has been exercised on toy data, synthetic planted-truth challenges,
and real-read diagnostics above 1 Gb. The exact-output native alignment
optimization reduced both time and RSS in two paired Mo17 replays; all core
products matched byte for byte. A later workspace-reuse idea had a time/RSS
tradeoff and was rejected. Full details and the unfavorable real-tool
comparisons are in `performance_scope_20260917.md`. The latest full local suite
after the rice v2 audit and five-source eligibility reconciliation passed
1,075 Python tests; hosted Ubuntu/macOS Python, Rust, and wheel validation at
commits `c3e8624` and `16538ca` passed.
Tests and CI establish source behavior under covered inputs, not biological
accuracy or 7–20-Gb production readiness.

## Scientific functions and claim limits

| Function | Evidence for present capability | Boundaries and unresolved claims |
| --- | --- | --- |
| Read-first repeat discovery | Direct read input; deterministic monomer/family catalogue; real-read runs, planted-truth synthetic recovery, exact-output engineering replays. | Real families, boundaries, low coverage, nested/shared arrays, and large-genome completeness are not established by the toy/synthetic truth alone. TideHunter is faster in several same-input diagnostics. |
| Read-derived family abundance | Diagnostic k-mer and read-evidence output; uncertainty and warnings; synthetic depth experiments archived. | The frozen estimator is not calibrated absolute physical copy number. Shared/identical signatures can be non-identifiable; genome-size normalization may use a provisional assembly-length proxy. M1 absolute abundance Gate A is NO-GO. |
| Assembly localization and comparison | Array candidate tracks, read-versus-assembly summaries, controlled injected edit frameworks. | Continuous discordance is an estimated read–assembly abundance deficit, not measured missing physical bases. A newer assembly is not absolute copy truth. Binary collapse classification and missing-base magnitude require separate evaluation. |
| Architecture/HOR | Candidate family hierarchy plus separately implemented research monomer-label path and transition graph; controlled structural errors. | Production hierarchy edges are heuristic, not validated HOR calls. Frozen M2 alignment route detected 6/44 B2 synthetic injected events and abstained on 45/52 cases; Gate B is NO-GO. Research M2 was not promoted. |
| FISH probe ranking | CLI ranks candidates and writes in-silico specificity outputs on toy and development inputs. | No prospectively validated FISH success rate or universal off-target guarantee is established. |
| Visualization/reporting | Offline HTML and static/editable vector products with underlying source tables. | Figure availability does not validate the biological interpretations rendered in them. |

The M1 tournament evaluated six research directions against strong baselines
under development conditions. None established a stable independent gain in
abundance calibration with safe shared-family assignment. A later complete-read
local-tile prototype lowered assigned-read-bp MARE from 0.5928 to 0.3410 on
one small synthetic development input versus an 80-bp chunked mapping adapter,
while wrong-family assignment increased from 4 to 7 bp. It is a different
endpoint from production genomic-bp output, and resource timers cover
different work. It cannot justify backend replacement or a speed/memory claim.
See `m1_full_read_research_20260917.md` and its hostile audit.

The M2 C3 original-read development test used supplied operational 178-bp
monomers and source-guided flank trimming. Its label-path comparison detected
8/8 engineered large deletions and supported one intact control; a simple
read-span rule made the same decisions. In a separate equal-length six-case
challenge, M2 detected two inversions but missed two swaps and one tile
replacement (2/5), while the length rule detected 0/5. The structural result
is a diagnostic of what this representation can and cannot encode, not an
equal-prior or biological superiority claim. The B2 synthetic gate remains
NO-GO, as documented in `method_gate_decision_20260917.md` and
`c3_equal_length_hostile_audit_20260917.md`.

## Native-source and controlled-truth status

Col-CEN v1.2 has three 3,560-bp nominated arrays in 13,560-bp contexts and
six unedited original HiFi molecules: C1 has two local spanners, C2 one, C3
three. All six reported one primary alignment over the expected locus and
natural flanks on the complete 132.081-Mb reference under the bounded mapping
protocol. Twenty-seven exact terminal/internal/boundary assembly edits were
independently reconstructed; the reads were unchanged. The 178-bp unit is
operational. Exact individual/extraction/haplotype pairing and physical copy
truth remain unverified; all are **development** inputs. Source and edits are
under `benchmarks/controlled_collapse/native_pair_audit_20260917/` and
`native_edit_development_v1/`.

Ey15-2 (9994) contributes a second source lineage: one selected 3,244-bp
array and seven original HiFi molecules whose full-assembly reported primary
alignments cover the locus. Nine exact edited assemblies have unchanged read
hashes and independent sequence reconstruction. The preregistered,
source-guided 5% read-span rule classified all eight large injected deletions
as discordant and the intact control as supported. Its smallest deletion was
811 bp versus a 163-bp threshold. These nine outcomes reuse one chosen array
and the same seven reads; they do not estimate independent-donor sensitivity
or performance near the threshold. The 420-bp operational period exceeds the
frozen M2 route's 300-bp cap, so M2 evaluated **0/9** Ey15 cases. Exact DNA
extraction pairing, physical copy count, and natural collapse truth remain
unverified. See `ey15_native_span_score_20260917.md` and independent
`ey15_native_hostile_audit_20260917.md`.

The complete controlled-collapse design covers precise injected bp edits,
but the two qualified current lineages are both Arabidopsis development
materials. They do not establish transfer to wheat, barley, rye, oat, maize,
or macadamia. No native final biological held-out cohort has been enrolled;
the existing B2 holdout is synthetic and was consumed once. Same-study,
same-accession, or same-published-sample evidence must not be counted as
independent biological replicates or absolute molecular truth.

## Competitor and resource interpretation

TRF and TideHunter have actually run in several same-input real-read discovery
diagnostics. TandemX was faster than TRF and slower than TideHunter on the
listed Chinese Spring, Victoria, and Morex examples; memory ordering changes
by dataset and competitor. The 1.170-Gb Morex comparison was concurrent,
not a final isolated throughput rank. `performance_scope_20260917.md` gives
all numeric values, units, and conditions.

TideCluster/TideHunter and CENdetectHOR have both completed one equal-FASTA
5,700-bp all-array interface control and returned an all-array window with a
60-bp period. Neither produced the designed separate 30-bp A/B order in
that control. On the more realistic 36,371-bp flanked 171-bp assembly input,
TideCluster has prior final localization evidence. CENdetectHOR has now been
rerun under its actual default and a fixed documented diagnostic setting;
both stopped before final decomposition/HOR output, so its accuracy endpoint
is **technical N/A**, not zero. The failures and failed-stage RSS are retained
in `cendetecthor_171_precondition_20260917.md`. A missing endpoint is not a
false negative, and small control resource figures cannot be scaled to full
chromosomes.

## Work still open or deliberately stopped

The principal missing evidence is orthogonal physical abundance and
independent-donor structural truth. Neither current source lineage provides
it. Final held-out real-array validation, large-genome end-to-end isolated
performance, cross-hardware reproducibility, experimental FISH validation,
and a task-equivalent successful HOR comparison remain uncompleted. A
non-identifiable family must be reported ambiguous rather than assigned a
precise copy count. `unitFinder` is stopped. Research M1/M2 prototypes remain
outside the production CLI. Formal tagged release, Bioconda, Zenodo, and
manuscript promotion are paused under the current project instructions.

The completed bounded Mo17 source screen selected one 3-kb CentC interval
from six assembly-only candidates before read mapping. Its current 81,775-read
subset passed complete hash, gzip and FASTQ checks. Six distinct original
reads geometrically spanned the nominated array and 1-kb flanks, but none
reached the frozen 95% whole-alignment identity cutoff (best 90.06%).
Flank identities were also low. The interval is therefore ineligible for the
three-read controlled-collapse development case, and no alternate locus was
selected after seeing the reads. Assembly and reads have different archive
BioSamples, and full-reference mapping/edit generation were not run. See
`benchmarks/controlled_collapse/mo17_native_pair_audit_20260917/README.md`.

The completed Nipponbare/AGIS1.0 bounded source screen verified the current
385.711-Mb assembly and 62,345-read original HiFi subset by complete hashes
and structure. Its precommitted 155-bp/4-kb interval selection examined
771,185 assembly windows and the top 1,000 ranked candidates; none met the
frozen 0.8 within-context flank 31-mer uniqueness threshold (maximum 0.638).
No array was selected and reads were not mapped in v1. A separately frozen
**development v2** used assembly-only information from the v1 top 1,000 to
evaluate 235,500 length/period combinations. It committed the selected
CP132242.1 5-kb context before read mapping. Of 62,345 original reads, only
one met the fixed >=99% identity and >=1-kb flank criteria, below the required
three distinct spanners. It therefore did not run full-reference anchoring or
make edits, and it did not reselect an interval. The distinct assembly/read
BioSamples and dates also leave exact donor pairing unverified. These outcomes
apply to their declared screens, not all rice arrays. See
`benchmarks/controlled_collapse/rice_native_pair_audit_20260917/` and
`benchmarks/controlled_collapse/rice_native_pair_audit_v2_20260917/`.

Macadamia jansenii has a published assembly/read relationship but no existing
bounded original-HiFi subset was found. The two full compressed read streams
total 18,150,429,537 bytes; only historical full-file hashes are available
and exact donor/extraction pairing is unresolved. Under the current USB2/EIO
condition, no full reread, scan, native interval, or edit was attempted. See
`benchmarks/controlled_collapse/macadamia_native_pair_audit_20260917/`.

The unified, machine-readable eligibility snapshot is
`benchmarks/inputs/native_source_registry_20260917.json`; it identifies two
technical development lineages and three negative/not-enrolled source screens.
It is a post-experiment reconciliation, not benchmark preregistration. Current
T7 access is through a 480-Mb/s USB2 path with prior I/O errors. Large
unverified source files and large T7 writes are excluded from current
benchmark input until integrity and storage reliability are established.
The previously identified K30076 block-2 attempt-001 gzip/chunk is invalid:
same byte count but changed SHA-256 and failed gzip, followed by T7 I/O
errors. It is excluded from every current benchmark. Reacquiring and
revalidating that input requires a reliable storage path; its failed bytes
must not be reclassified as biological zero coverage.
