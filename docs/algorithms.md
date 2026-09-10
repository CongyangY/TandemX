# TandemX Algorithm Design

For `backend=rust`, the representative candidate gate now uses an injective
ACGTN word encoding and native exact multiplicity accumulation. Python remains
the reference gate. Complete monomer sequences now cross the native boundary,
so circular word counting no longer creates a Python string `Counter` for every
candidate. Identical rounding and fixed representatives preserve the existing
alignment and membership decisions; no biological threshold changed.
See the native-index ablation in `clustering_resource_replay.md`. Measured
performance remains separate from algorithmic equivalence.

Source-backed repeat queries are curated by accession/version, complete-record
checksums, species/material and explicit monomer coordinates. Multiunit clones
are excluded from monomer scoring. These post hoc queries are not a complete
truth catalogue; see `known_monomer_sources.md`.

Representative postings now use contiguous 64-bit items with exact 32-bit ID and
multiplicity fields instead of one Python dictionary per word. Their numerical
meaning and all candidate selection/alignment rules are unchanged; oversized
fields fail explicitly. The index is released before constructing output records.
See [isolated resource replay](clustering_resource_replay.md) for source-matched,
fresh-process parity and memory checks. No whole-pipeline speed/memory claim
follows from compact storage alone.

Real-input benchmark normalization now defaults to a disk-backed read-length and
interval index (`benchmarks/scripts/real_disk.py`). Native adapters yield rows;
coordinate validation and interval unions preserve the previous endpoint, call
order and duplicate policy. Python state is bounded by one record/interval;
SQLite uses a 16 MiB page-cache target and on-disk indexes. This is evaluator
scaling, not a change to discovery or biological interpretation. See
[real comparator replay](real_comparator_pilot.md).

Sequence clustering now accumulates canonical circular 9-mer **multiplicities**
through its inverted representative index. For identity >.9 and query length>=20,
it rejects length-incompatible pairs and canonical multiset overlap below the
existing oriented q-gram bound. Collapsing a word with its reverse complement can
only increase overlap, so a pair rejected here would fail both unchanged oriented
tests. Integer edit-budget rounding, candidate order, witnessed alignments,
representatives, ties and ambiguity records are preserved. Lower identity/short
queries keep exhaustive selection. The index stores multiplicity per representative;
its memory cost and dense-graph scaling still need real measurement.

`benchmarks.scripts.replay_clustering --previous-run /path/to/completed/real_run
--outdir /path/to/new/replay` compares original snapshot clustering with the new
index on identical exported candidates. Exported scores are rounded, so this
tests two algorithms on serialized evidence, not complete live-pipeline parity.
Exact family/membership JSON, source/input hashes and diagnostic stage times are
retained. No per-stage memory ranking is inferred from this single-process replay.

The experimental [read-cluster ratio model](read_cluster_quantification.md)
handles within-read k-mer dependence and finite read-end opportunities using
streaming family moments. It is a Python research reference, separate from the
default quantify estimator, with explicit missing intervals for sparse support.
The separate [multi-k attenuation prototype](multik_quantification.md) fits
word-survival trends across fixed k values, with mean-k21 and native-median
ablations. It does not infer sampling intervals from correlated k values.
The [joint-read extension](joint_multik_uncertainty.md) propagates read-level
cross-k moments through the fitted intercept; its approximate interval requires
independent reads and does not include extrapolation or catalogue bias.
Full-file [QC and nested sampling](complete_data_qc.md) support unbiased selection
within included libraries and record the limits of library/material sampling.
The [reference-mapping QC](reference_mapping_qc.md) adds a disk-backed PAF audit
and CIGAR-based coverage summaries; mapping concordance is not accuracy truth.

The [real-input comparator pilot](real_comparator_pilot.md) verifies selected
FASTQ hashes, converts identical inputs for all tools and reports descriptive
interval unions. It does not use cross-tool agreement as biological truth.

This document describes current MVP algorithms and planned future algorithms. The repository currently implements the toy simulator, toy-scale `discover`, `quantify`, `locate`, `probe`, and `visualize` MVPs. The default workflow is de novo: reads are passed to `tandemx discover`, and its output catalog is reused by downstream commands.

## Candidate Periodic k-mer Discovery

Valid nonempty input with no candidates or no families above the support
threshold now produces a completed zero-result run. A `discovery_summary.json`
receipt records processed input counts, zero output counts and output SHA-256
hashes. Validation allows empty discovery tables/catalogs only against an intact
receipt. Empty/malformed reads and a selection that processes no reads remain
errors. Pipeline dependencies are skipped with an explicit no-family reason.

The challenge evaluation is independent of these analysis algorithms. Array
true positives require the same read, interval IoU >=0.5 and period error <=
max(2 bp, rounded 2% of the truth period), using maximum-cardinality one-to-one
matching. Read detection ignores exact period/boundaries and is a separate
endpoint. Family recovery uses one-to-one matching between distinct, strand/
rotation-canonical predicted consensuses and truth monomers, with equal length
and exhaustive circular ungapped identity >=0.90. Unequal-length consensuses
fail this deliberately strict endpoint; it is not an indel-aware homology test.
TandemX is scored from its final catalog, and per-array finders from their
consensuses. No family precision is inferred from catalog size or length alone.

MVP goal: identify simple candidate tandem repeat monomers de novo from toy HiFi-like sequence reads.

Default `--discovery-method legacy` implementation (retained for ablation):

1. stream FASTA or FASTQ reads, including gzip-compressed inputs;
2. apply `--max-reads`, `--max-read-bases`, reproducible `--sample-rate`, minimum length filters and the automatic large-input discovery budget when enabled;
3. directly scan bounded short periods so STR-like 2-19 bp repeats can be retained;
4. extract canonical, non-low-complexity k-mers for one read at a time with a rolling 2-bit encoder in the selected Python or Rust backend for longer periods;
5. retain only repeated within-read k-mer positions with strict position/pair caps;
6. build a spacing histogram in the configured period range;
7. retain only `--top-periods` supported spacing peaks;
8. add plausible fundamental divisors of spacing peaks and refine periods against the strongest local shifted-identity interval in linear time;
9. require evidence spanning at least two units and report inferred local array boundaries instead of assigning the whole read;
10. build a strand- and phase-canonical monomer consensus from complete units in that interval;
11. append accepted candidates immediately to `candidate_reads.tsv`;
12. cluster candidates only when both period and circular sequence sketches are compatible, preventing unrelated same-length repeats from being merged;
13. build a deterministic cyclic consensus for each supported cluster and write `monomers.fa` and `families.tsv`.

`monomers.fa` is an output of de novo discovery. It is not a required input to `tandemx discover`.

The previous full period scan has been removed from the execution path. `--kmer-backend auto` is the default and uses the Rust backend when the compiled extension and requested k-mer size are supported, otherwise it falls back to Python. The Rust backend implements the same single-read seed/spacing/refinement boundary and returns only a compact result to Python; parsing, clustering and output remain in Python. With Rust, read-local scanning can run across multiple threads while preserving deterministic output order. Neither backend makes TandemX a full large-plant-genome production workflow. See `docs/performance.md`.

For large real-read inputs, `discover` can stop early by design, but only when the user opts in. If `--enable-auto-discovery-budget` is set and the user has not set explicit `--max-reads` or `--max-read-bases` limits, TandemX enables an automatic discovery budget in either of two cases:

1. `--genome-size` is provided, in which case discover targets approximately `--target-discovery-coverage` genome equivalents, bounded by `--auto-discovery-max-bases`;
2. the counted input exceeds `--auto-discovery-trigger-bases`, in which case discover uses a conservative absolute cap of `--auto-discovery-max-bases`.

This bounded mode is intended to keep monomer discovery in the subset regime rather than exhaustively scanning every read in 30-50X HiFi data. When multiple files are supplied and the automatic budget is active, reads are pulled in round-robin order across files so early stopping is less biased toward the first file.

MVP constraints:

1. toy data only;
2. simple tandem arrays only;
3. FASTA/FASTQ/gzip parsing, multiple read input files and streaming progress are supported;
4. one best candidate period per read;
5. short-period and low-complexity candidates are retained with warnings, not interpreted as high-confidence long satellite arrays by default;
6. local boundaries are an identity-based estimate, not a complete array reconstruction;
7. consensus is a cyclic column majority, not a full indel-aware multiple alignment;
8. no higher-order repeat inference;
9. no multiprocessing, intra-step checkpointing or production-scale full-workflow backend yet.

Future work:

1. robust clustering for related satellite families;
2. strand-aware consensus refinement;
3. multiprocessing or distributed chunk processing for non-Rust backends;
4. uncertainty modeling for ambiguous monomer periods.

## Elastic read-local discovery (experimental)

`discover` and `run` accept `--discovery-method elastic`. This opt-in method
uses the same input-only discovery contract and downstream schemas. The legacy
mode remains the default pending independent evaluation and resource profiling.

1. Select distinct supported spacing bands from bounded per-read seeds. Rust mode
   now reuses the native bounded seed extractor and spacing histogram already
   used by the legacy kernel; Python remains an independent reference. Returned
   histogram counts and overflow counts are tested for exact equality, including
   pair caps, Ns, low complexity, reverse complements and short-period-only scans. Nearby
   peaks covered by an existing band do not consume another `--top-periods` slot.
   Short periods 2–19 are scanned directly when requested.
2. Align a read to itself around each positive offset. Band half-width is
   `min(period-1, max(3, ceil(0.08*period)))`. A diagonal pair scores +2 for a
   match, -3 for a mismatch; each gap base costs 4. Ambiguous bases break paths.
   A path resets at score <=0 or a drawdown >40 from its own peak. This is a
   banded, drawdown-pruned local alignment heuristic, not an exact unbounded
   Smith-Waterman optimum.
   The native kernel updates score and path-peak rows in place, left to right:
   the current band and its right neighbour still hold the previous row, while
   the left neighbour already holds the current row. Rejected cells are reset;
   the inactive right tail is cleared after the last upward dependency is used.
   This removes per-row allocations without changing scores, tie order, traceback
   or acceptance rules. Differential tests check complete paths against Python.
   Native traceback directions are packed into two bits per cell because the
   recurrence has only stop, diagonal, left and up states. Candidate periods for
   one read are submitted together, allowing the native code to uppercase and
   copy the read once before evaluating each period in the original order. These
   are storage and call-boundary changes: scoring, bands, tie order, traceback
   and accepted paths remain unchanged. Complete discovery replays require
   byte-identical core outputs before a resource difference is reported.
3. Trace candidate local paths in descending score order. Require alignment
   identity >=0.75, the configured span, and a comparison span >=0.8 times the
   median aligned offset. The minimum alignment score is
   `max(min(40,min_span), ceil(0.7*max(period,min_span-period)))`.
4. Filter composition-driven matches using `(identity-q)/(1-q) >=0.7`, where
   `q=sum(f_b^2)` from A/C/G/T frequencies in the local interval. All-one-base
   intervals fail this filter. This score is not a calibrated probability and
   does not account for dinucleotide dependence or local-hit selection.
5. Keep multiple spatially distinct arrays. If two hits overlap by at least
   half of the shorter interval, keep the higher total alignment score, then
   longer span, then smaller period. This suppresses redundant harmonic calls;
   it does not establish that an array has no higher-order organization.
6. Follow aligned adjacent-copy coordinates to cut approximate unit boundaries.
   Omit short terminal fragments; a terminal segment within 10% of the median
   offset can be included. Use at most 32 evenly spaced observed units. Select
   a median-length template with high cross-unit 7-mer support, then perform up
   to two deterministic rounds of banded global alignment and majority voting,
   including insertions and deletions. This is an observed-unit consensus, not
   ancestral sequence inference or full partial-order assembly.
7. Canonicalize strand/rotation and, by default in elastic mode, apply the
   operational sequence clustering described below. `period_bp` is consensus length; a differing median aligned offset
   is retained in `warning`. `score` is matches divided by all aligned columns,
   including gaps. Confidence labels remain uncalibrated.

The independent Python reference and Rust kernels are tested for identical
alignment paths and consensus. Rust releases the GIL during alignment; seed
histogram construction, local/global alignment and bounded edit comparison
release the GIL in native Rust. Period selection, consensus voting and cluster
orchestration remain in Python. Each trace is
limited to 32 million cells and fails explicitly if exceeded. Peak trace space
is O(read_length × band_width), plus read-local alignment pairs; no whole read
collection or genome is loaded for the alignment. Candidate/family state still
grows across reads and is a separate production-scale limitation.

Development findings are archived separately from the baseline; both rejected
first attempts and corrected results remain available. Perfect recovery on
development seed 1101 is not a held-out accuracy or speed claim. The expanded
comparator/metric contract is in [comparator_matrix.md](comparator_matrix.md).

`--discovery-method cascade` first calls the same native seed/period screen. A
gap-free proposal can bypass elastic traceback only when its interval spans at
least 30% of the read, shifted identity is at least 0.95, at least 95% of shifted
columns contain A/C/G/T, its span is within 2% of an integer number of proposed
units, and composition-adjusted identity is at least 0.7. The unit-span and
valid-column guards were added after an intermediate development rule worsened
0.1% indel boundary MAE from 1.364 to 14.221 bp. Reads outside all five gates use
the unchanged elastic path, preserving short-array, phase-drift and multiple-array
behavior. These constants were selected on seed 1201. After commit `054b935`
and both hosted workflows passed, frozen seed 2201 was run once and passed all
14 predeclared validation gates. The result supports the guarded split on this
synthetic distribution; it does not establish real-data or universal superiority.

## Operational monomer clustering (experimental)

`--clustering-method auto` selects `sequence` for elastic discovery and the
historical `legacy` clusterer for legacy discovery. Either can be requested
explicitly for ablation. `--cluster-identity 0.95` is an operational sequence
resolution, not a universal definition of a biological satellite family.

Exact canonical candidate sequences are grouped and ordered by distinct-read
support, total supporting span, mean read-local alignment identity, then sequence.
Each group is compared to **fixed observed representatives** in abundance order.
Circular global unit-cost edit similarity must be at least the configured value:
`1 - edit_distance / max(lengths)`. Both strands and every rotation of the query
are eligible. Matches require known A/C/G/T bases; N is never matching evidence.
Length difference and shared circular q-gram counts provide necessary-condition
rejections. Seed-supported rotations are tried first; a thresholded, banded
Python/Rust global edit distance verifies a merge. Every remaining rotation is
tried before rejecting pairs that pass the q-gram bound. Acceptance can stop at
its first valid witness: reported distance is an **upper bound** and similarity
a **lower bound**, not a claim to the optimal alignment.

Representatives do not drift and clustering is not transitive. When multiple
existing clusters pass, choose the highest witnessed similarity, then earlier
representative rank; preserve the alternatives. Thus this is not an all-pairs
nearest-neighbour assignment. Keep clusters meeting distinct-read support;
several arrays from one read count as one supporting read. The observed consensus
representative is used as the output monomer, avoiding cross-read consensus
chimeras. Broader families can contain several operational monomer clusters.

`candidate_monomers.fa` preserves every read-local candidate sequence even when
its cluster fails support filtering. `monomer_membership.tsv` records assignment,
representative hash, similarity bound, alternatives and support-filter fate.
Candidates with an N fraction above the allowed edit fraction remain unassigned
with `unresolved_sequence`; they are retained in the candidate FASTA. Lower N
fractions remain penalized, including comparison with a representative itself.
Confidence labels are heuristic and uncalibrated. Optional historical redundancy
collapse is a separate output and must not replace this evidence silently.

Native comparisons release the GIL and use band-width memory; comparisons above
32 million potential cells fail explicitly. Candidate and cluster state still
grows with the number of candidates. This does not establish production-scale
memory bounds for whole libraries.

## Representative-pair redundancy audit

After clustering, `family_similarity.tsv` examines every distinct pair in stable
catalogue order. Cache each representative's canonical k-mer set once; the
selected Python/Rust backend then computes the same all-offset ungapped local
identity on both strands, with minimum overlap `min(50, len(a), len(b))`. Ties
prefer greater overlap and then the first orientation/offset. This historical
heuristic counts equal ambiguous characters as matches; it is not the gapped
circular clustering identity (which penalizes N), nor a homology truth label.

Rows stream to a `.partial` file that is renamed only after a successful audit.
Retain related-pair warnings, and only the `likely_redundant` edges when optional
collapse is requested. This avoids retaining every pair and constructing a full
output string in memory. It does not eliminate quadratic time/output size or the
potentially growing warning/related-edge state. The native kernel releases the
GIL. Randomized/edge Python parity checks cover scores, overlaps, orientation,
table bytes, warning order and optional collapse. Rebuild the extension before
running this version; a stale extension fails explicitly.

The first complete-library Mo17 random pilot exposed this bottleneck: 840 reads
scanned in 14.514 s, followed by 179,101 pair comparisons for 599 families; total
TandemX wall time was 258.829 s. Native/streamed audit rerun took 17.348 s and
81.70 MiB RSS versus 168.58 MiB initially; all six outputs were byte-identical.
These single runs overlapped acquisition/QC and are diagnostic, not final ranking.

`--family-audit related` builds an inverted index of exact canonical k-mer sets.
For one representative at a time, count intersections with later representatives
in stable order. All existing non-distinct rules require a nonzero intersection.
Evaluate the rules with identity/overlap fractions set to one: if even this upper
bound is distinct, ungapped alignment cannot produce a related pair. Otherwise
compute the ordinary exact alignment, and emit only actual non-distinct pairs.
There is no approximate sketch or sequence-identity cutoff added to this filter.
Changing relationship rules requires rechecking this monotonicity condition.

The primary catalogue, its warnings and optional collapse remain identical to
full mode. The pair table intentionally omits distinct rows and has a companion
`family_audit_summary.json` with possible/scored/pruned/emitted counts. The full
mode remains default. Index and candidate state grow with the catalogue; highly
similar catalogues can still require quadratic work/output. Empty catalogues
write a header-only table and a zero-count receipt. This mode does not change
read-local detection or gapped sequence clustering.

The audit simultaneously streams `family_hierarchy.tsv`. Each
`possible_higher_order_or_partial` relationship becomes an edge directed from
the shorter to the longer representative. A nearest integer length multiple of
at least two, with absolute ratio error no greater than 0.05, is labelled a
`putative_period_multiple`; other related edges remain unresolved. This preserves
all supported 171/342/684-style alternatives as a graph instead of forcing one
parent or collapsing period multiples. The graph is a catalogue QC result. It
does not reconstruct HOR order, prove direct ancestry, or distinguish a true HOR
from harmonic detection and partial representatives without additional context.

## Diagnostic k-mer Copy-number Calibration

MVP goal: estimate repeat family copy number from diagnostic k-mer depth on toy reads using a repeat catalog discovered upstream.

Current MVP implementation:

1. read the monomer catalog generated by `tandemx discover` or a supported
   provenance-labelled import from `--catalog/--catalogue` or `--monomers`;
2. enumerate canonical k-mers in circular tandem context from every monomer phase, including when `k` is longer than the monomer;
3. count k-mer multiplicity within each monomer;
4. remove low-complexity k-mers and k-mers shared by multiple families;
5. optionally read independently selected background controls from
   `--single-copy-kmers`; each row declares a canonicalizable k-mer and expected
   haploid copy number, and overlap with a repeat diagnostic k-mer is rejected;
6. stream canonical target k-mers from FASTA or FASTQ reads, including
   gzip-compressed inputs;
7. estimate the probability that a complete k-mer survives base-calling error.
   FASTQ uses the mean product of Phred-derived per-base correctness across all
   valid windows; FASTA can use `(1 - --read-error-rate)^k`. This independent-
   error approximation is explicit and can be disabled;
8. divide observed target counts by this survival probability, then correct each
   diagnostic depth by its multiplicity within the circular monomer;
9. summarize corrected diagnostic k-mer depth with the median;
10. use `--haploid-depth` when supplied; otherwise use the arithmetic mean of
    all corrected single-copy control depths, including zero observations, then
    fall back to total read bases divided by `--genome-size`. The mean avoids the
    zero-median discontinuity near 1x coverage; median, MAD and zero fraction are
    retained as control-QC fields. Optional `--single-copy-min-depth` applies a
    condition-level observed-control-depth gate: below the threshold it uses the
    total-bases fallback while retaining control statistics and an explicit
    fallback method/warning;
11. estimate copy number as `median_kmer_depth / haploid_depth` and repeat bp as
    `estimated_copy_number * monomer_length`;
12. report raw/corrected depth, normalization method, survival probability,
    control count/mean/median/MAD/zero fraction, diagnostic median absolute
    deviation and empirical 10th-90th
    percentile copy-number endpoints;
13. batch Rust target counting while releasing the Python GIL; both backends
    retain only requested diagnostic/control target counts, not all read k-mers.

MVP constraints:

1. toy reads only;
2. no complex ploidy model;
3. TandemX does not infer that user-supplied controls are truly single-copy;
   their selection and assembly/population provenance must be recorded;
4. no external k-mer counter;
5. without explicit depth or controls, normalization uses a rough total-bases/
   genome-size estimate and is labelled with a warning;
6. confidence labels are based on diagnostic k-mer count, depth dispersion and whether haploid depth was provided;
7. repeat diagnostic uniqueness is checked within the supplied catalog, not
   against an independent genome background, so the output emits
   `genome_background_uniqueness_not_verified` and does not assign high confidence;
8. Phred values are treated as calibrated independent error probabilities and
   the 10th-90th percentile endpoints are diagnostic-k-mer dispersion, not a
   calibrated biological confidence interval.

The downstream `--catalog` input reuses discovery results. It does not mean TandemX requires repeat sequences before de novo discovery.

The paired development design for the optional normalization paths is documented
in [quantify_calibration.md](quantify_calibration.md). It scores public-command
outputs and retains the oracle-input boundary rather than treating planted error
rates as observable data. The value 2.0 was frozen after development and passed
nine predeclared gates on three additional simulated genomes; this supports the
explicit opt-in mode but does not validate real control specificity or make it a
default.

Future work:

1. contamination-aware k-mer filtering;
2. automatic, population-aware single-copy control selection and QC;
3. depth modeling across multiple samples;
4. independently calibrated bootstrap, negative-binomial or Bayesian intervals;
5. scalable k-mer counting backends.

## Assembly Density Localization

MVP goal: locate simple repeat evidence on a toy assembly and summarize density in windows.

Current MVP implementation:

1. stream assembly FASTA in bounded sequence chunks and retain `k-1` bases across chunk boundaries so cross-chunk k-mers and coordinates remain exact;
2. enumerate circular, non-low-complexity canonical k-mers and build one code-to-family index;
3. exclude k-mers shared by multiple catalog families, then scan each assembly contig once with a rolling 2-bit encoder;
4. merge matching k-mers online into intervals instead of retaining every hit; the default permits a `2k` unhit gap, while `iid_base` permits at most one monomer length to bridge mutation-broken anchors;
5. finalize nearby intervals for each family as scanning advances;
6. filter very short intervals and enforce `--min-identity` using either the default exact diagnostic-k-mer fraction or the optional assumption-limited `iid_base` proxy;
7. write candidate arrays as 0-based half-open `arrays.bed`;
8. compute the union of family intervals, then use prefix coverage and binary search for sliding windows so overlaps are not double counted and density remains in `[0, 1]`;
9. write `repeat_density.bedgraph`;
10. if `copy_number.tsv` is provided, call `tandemx.compare.mvp` to write a backward-compatible `assembly_vs_read_cn.tsv`.

MVP constraints:

1. toy assemblies only;
2. k-mer evidence only, no read mapping;
3. `--min-identity` is a k-mer-derived proxy, not alignment identity or an exact per-copy placement claim; `iid_base` converts an exact-k-mer fraction `f` to `f^(1/k)` under the independent-substitution relation `P(exact k-mer) = P(base match)^k`, which overlapping windows, shared-k-mer removal, indels and structured variation can violate; its one-monomer gap bridge can join nearby homologous tracts and is reported in each array warning;
4. simple threshold classification for `possible_collapse` and `possible_overexpansion`;
5. bigWig output is future work.

Future work:

1. alignment-backed localization;
2. chromosome-scale streaming;
3. family-specific tracks;
4. confidence labels for ambiguous or multi-family hits.

## Output Validation

MVP goal: check that TandemX core outputs match documented schemas before downstream use.

Current MVP implementation:

1. scan a project directory with `tandemx validate --project`;
2. validate required TSV columns for candidate, family, copy-number, comparison, probe, and FISH tables;
3. validate numeric fields as integers or floats while streaming files line by line;
4. require documented confidence, status, and warning fields where defined;
5. validate `arrays.bed` and `repeat_density.bedgraph` as 0-based half-open intervals;
6. validate TandemX FASTA headers for `monomers.fa` and `probes.fa`;
7. fail clearly on empty recognized output files.

MVP constraints:

1. schema validation does not prove biological correctness;
2. validators check documented field structure, not full cross-file consistency;
3. recognized filenames are fixed for the current MVP.

## Assembly-vs-read Comparison

MVP goal: compare toy read-based and assembly-based repeat abundance.

Current MVP implementation:

1. read `copy_number.tsv` and use `estimated_bp` as the read-based family abundance;
2. read `arrays.bed`, merge overlapping intervals per family and contig, and sum interval-union lengths without double counting;
3. calculate `assembly_read_ratio` as assembly-estimated bp divided by read-estimated bp;
4. classify status as `consistent`, `possible_collapse`, `possible_overexpansion`, `assembly_only`, `reads_only` or `low_confidence`;
5. write `assembly_vs_read_cn.tsv` from `tandemx compare`.

`repeat_density.bedgraph` is not the main input because it lacks `family_id` and cannot support family-level comparison.

MVP constraints:

1. simple threshold-based classification;
2. no claim of proof of collapse;
3. no complex assembly quality model;
4. no multi-sample population comparison.

Current cohort integration:

1. `tandemx cohort` reads each sample's monomer catalogue, copy-number table,
   and optional assembly/read comparison;
2. it clusters circular monomers against fixed, abundance-ordered observed
   representatives with witnessed edit-similarity bounds and no transitive
   membership propagation;
3. it aggregates copy-number interval endpoints in bp and assembly/read ratios
   only after local families are mapped to a pan family;
4. it emits explicit `NA`, `not_observed`, and `not_evaluated` states rather
   than treating missing evidence as zero abundance;
5. it rejects duplicate, unknown or incompletely quantified catalogue families,
   hashes every supplied input and writes a per-sample QC audit;
6. it writes point, lower-endpoint and upper-endpoint abundance matrices and
   propagates the least confident local evidence into each aggregate;
7. it ranks pan families by maximum sample read-estimated bp with a stable ID
   tie-break, then renders the requested top families as a four-panel editable
   SVG/PDF overview with a machine-readable plot-source table and hash receipt.

The interval panel displays the sum of local marginal endpoint widths relative
to the point estimate. It is a QC diagnostic, not a calibrated joint interval
across families or samples. Display clipping of assembly/read ratios does not
change values in `cohort_plot_source.tsv` or the long-form source table.

Remaining work includes chromosome/subgenome summaries and calibrated joint
uncertainty across biological samples.

## TideHunter catalogue import

`tandemx import tidehunter` separates detector choice from TandemX's catalogue
and downstream analysis. It performs the following steps:

1. stream all original read files into a temporary SQLite index containing only
   the normalized ID, native header ID and sequence length while calculating an
   ordered ID/sequence digest;
2. strictly parse every 11-field TideHunter `-f 2` row, including late rows,
   and require finite copy/identity values, valid consensus DNA, matching
   consensus length, increasing subunit starts and valid 1-based coordinates;
3. match the full TideHunter header ID to the source index, verify native read
   length, and convert the interval to 0-based half-open coordinates;
4. retain native fields in `tidehunter_import.tsv` and create standard
   TandemX candidates whose score is the native average-match percentage and
   whose unit count is the native copy-number field;
5. cluster native consensuses with the existing fixed-representative circular
   edit model, then run the normal family relationship audit and write standard
   catalogue outputs and hash receipts.

The temporary source index bounds ID/length memory, but imported candidate
consensuses and clustering state still grow with the number of TideHunter calls.
Confidence labels remain operational and uncalibrated. The importer enables
TandemX quantification, localization, assembly comparison, probe prioritization
and cohort analysis after TideHunter detection; it does not convert external
calls into TandemX native detector evidence.

## FISH Probe Scoring

MVP goal: rank simple candidate probes for toy repeat families.

Current MVP implementation:

1. read the de novo monomer catalog, assembly FASTA, `copy_number.tsv`, and `arrays.bed`;
2. generate deterministic windows using `--min-len` and `--max-len`; when a monomer is shorter than the requested probe, generate the probe from tandem circular context;
3. exclude probes with high single-base low-complexity ratio;
4. compute probe length, GC content and a salt/formamide-adjusted long-oligo Tm approximation;
5. index all probe k-mers together and stream each assembly contig once to find predicted target and off-target regions;
6. treat hits overlapping same-family `arrays.bed` intervals as target-array hits;
7. estimate `arrayiness_score` as target hits divided by all predicted hits;
8. estimate `specificity_score` as `1 / (1 + off_target_hits)`;
9. combine normalized copy number, specificity, arrayiness and GC balance into `probe_score`;
10. write ranked probes and toy in silico FISH predicted signal regions.

MVP constraints:

1. score is a prioritization heuristic;
2. no guarantee of experimental FISH success;
3. the Tm value is an approximation, not a thermodynamic hybridization model;
4. no full off-target alignment;
5. no probe tiling optimization;
6. no experimental calibration.

Future work:

1. off-target search against full assemblies;
2. probe tiling and multiplex design;
3. empirical calibration against published and new FISH data.

## Static Visualization

MVP goal: generate basic publication-oriented static summaries from existing TandemX toy outputs.

Current MVP implementation:

1. use matplotlib with a non-interactive backend;
2. generate `catalogue_summary.svg/pdf` from copy-number, assembly comparison and probe-ranking tables;
3. generate `assembly_vs_read.svg/pdf` from `assembly_vs_read_cn.tsv`;
4. generate `in_silico_fish.svg/pdf` from `in_silico_fish.tsv`;
5. keep plots simple and deterministic.

MVP constraints:

1. no web dashboard;
2. no interactive plots;
3. no pixel-level tests;
4. no seaborn dependency;
5. figures summarize toy outputs only.

### Authorized reporting and recovery extension (2026-09-10)

The current reporting extension reads existing family/abundance/localization
evidence and exports an offline HTML view and static SVG/PDF/PNG. It introduces
no detector, classifier or biological HOR rule. Summary deficits are zero-bounded
differences of existing estimates; composite confidence retains the least
available stage confidence. Candidate graph edges preserve sequence evidence
and ambiguity. Figure selection affects display only, with full tables retained.

The simplified run interface uses the existing fixed diagnostic k=21 and host
thread cap. Advanced YAML is strictly validated, with explicit CLI precedence.
Automatic resume requires validated stage outputs and matching input/command
fingerprints. An assembly total may provide a visibly provisional normalization
denominator when genome size is omitted; this is not a genome-size estimator.
Reads-only discovery remains available without a denominator, while absolute
abundance remains unavailable. Existing explicitly normalized analyses retain
their numerical algorithms and thresholds.

The recovery PoC uses frozen historical loci and old-assembly flanks. Flanks are
tested for assembly-relative uniqueness, then full original HiFi reads are
mapped to those anchors and tandemized family targets. The closest eligible
left/right pair is selected from old evidence before read outcomes. At least
three distinct, correctly ordered, same-strand dual-anchor reads with concordant
span lengths can nominate an observed median-length fragment. This is an
unpolished candidate, not a consensus assembly, locus-specific copy truth or
automatic patch. Candidate span and repeat-specific recovered bases are separate
fields. No unique anchors, missing old loci, sparse support and conflicting paths
remain explicit unresolved states. The newer assembly is prohibited during
candidate generation and reserved for a separate post-lock reference-proxy
comparison. See `docs/recovery_and_reporting_plan.md` for the fixed rules,
retained resource failures and the stop criterion.


### Exact Python selected-k-mer counting (2026-09-10)

For k <= 31, Python selected-target counting uses rolling canonical 2-bit codes
and converts only enrolled target keys back to their existing string identities.
The >31 and non-ASCII sequence cases retain the string iterator. The count,
quality correction, normalization and output fields are unchanged. This avoids
constructing a substring and reverse complement for every non-target window;
it is an exact implementation optimization, not a new abundance estimator.
See `performance_opportunity_audit_20260910.md` for the small-input parity and
resource experiment. Rust behavior is unchanged. The experimental periodic-
context gate under `benchmarks/challenge/` is not called by public commands.
