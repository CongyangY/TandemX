# TandemX: read-level evidence and sampling uncertainty for plant satellite repeat analysis

**Evidence-backed development manuscript, 6 September 2026.** Author information
is not assigned. This draft contains completed results, with unresolved release
and biological-validation requirements listed in `submission_readiness.md`.
It is not a submission-ready manuscript or a claim of universal superiority.

## Abstract

Satellite repeat analysis requires distinguishing repeat detection, family
identity, copy-number estimation and representation in an assembly. TandemX
combines streaming read analysis with sequence-supported family catalogues,
diagnostic k-mer quantification and explicit evidence records. We evaluated
development versions using independent simulated genomes and public plant HiFi
libraries. Three 10-Mb simulated genomes generated 27 coverage/error conditions
and 1,485 family conditions. In the 20× stratum with 2% unit divergence and high
read errors, fixed multi-k extrapolation reduced mean absolute relative
copy-number error from 56.59% for the median-k21 baseline to 13.53%. A joint-read
sampling calculation retained dependence across k values: at 20×, 424 of 495
family conditions supported intervals, of which 401 contained truth (94.58%).
Low-coverage missingness remained substantial. File-level QC covered 262.731 Gb
in ten libraries across eight reported plant species. Reference concordance in a 118.497-Mb
Arabidopsis subset identified 22.19% of input bases with organellar primary
alignment spans, highlighting a potential total-library normalization bias.
On three fresh predeclared conditional genomes, a frozen multi-k/depth rule
improved under-representation sensitivity from 81.48% to 85.60%, false-positive
rate from 8.64% to 7.41% and precision from 93.40% to 94.55%. The 20×/1%-error/
50%-retention stratum improved from 3/9 to 9/9, but one genome retained a 22.22%
false-positive rate. On a second predeclared divergent-array test, bounded anchor
bridging achieved 98.10% full-assembly recall and 99.96% positive-assembly
precision. The frozen multi-k rule raised classification sensitivity from 60.01%
to 91.84%, while false-positive rate rose from 1.23% to 8.54% and precision fell
from 98.65% to 94.16%.
Exact native optimizations shortened a 1.129-Gb maize replay by 23.61% relative
to its preceding indexed version while preserving seven output files byte for
byte; peak memory increased 0.81%. External comparisons showed task-dependent
trade-offs and did not establish universal speed, memory or accuracy superiority.
These results support further evaluation of read-based repeat evidence and
conditional uncertainty, while independent biological truth, held-out
generalization and production-scale validation remain necessary.

## Background

Tandem arrays pose several distinct analytical questions. A method can detect
repetitive sequence without identifying its smallest repeat unit; recover a
family without correctly delimiting every array; or localize an array whose
assembly copy number differs from the underlying sample. Treating these
outcomes as a single accuracy measure obscures both useful performance and
important failure modes. TRF established de novo detection of approximate
tandem repeats without supplying a known pattern [1]. Subsequent methods
address long-read repeat detection, assembly annotation and satellite-family
discovery through different representations and input assumptions [2–4].

High-contiguity plant references provide new opportunities to compare raw-read
and assembly evidence. They do not remove the need to inspect difficult-region
provenance. For example, the complete Nipponbare study describes a modelled
terminal 45S rDNA sequence; that region cannot automatically be treated as exact
copy-number truth [5]. A benchmark must therefore distinguish assembly-derived
coordinates, estimated copy counts, independently supported repeat families
and orthogonal biological observations. Technical library batches and nested
read samples also do not constitute independent plants.

TandemX is being developed to connect these evidence layers. Its intended scope
is candidate monomer discovery, repeat-copy estimation, assembly localization,
possible under-representation assessment and probe prioritization. It does not
assume that HiFi reads reconstruct all megabase satellite arrays. The current
study evaluates identifiable components of that programme and preserves
unfavourable results alongside improvements. No learned model is included in
the reported estimator, and no AI-based novelty claim is made.

## Results

### An auditable read-first implementation

TandemX exposes discovery, quantification, localization, probe, comparison and
visualization commands, together with validation and example workflows. The
development discovery path uses local periodic evidence and elastic alignment,
followed by sequence-supported family clustering. It preserves candidate
sequences and membership instead of equating matching period lengths with
homology. Operational clustering thresholds define analytical families; they
do not establish evolutionary ancestry. Run configuration, logs, stable tables
and source/input hashes connect results to their implementation.

The native clustering gate rejects candidate comparisons that cannot satisfy
the existing length and canonical-word multiplicity bounds. A subsequent
alignment workspace change retains only the required score/peak rows. These
changes preserve alignment and assignment rules. On 81,775 maize reads totalling
1,128,793,699 bp, the updated implementation returned 85,663 candidates and
28,586 operational families. All seven principal products were byte-identical
to the preceding indexed run. Elapsed time decreased from 1,625.652 to 1,241.876 s
(23.61%); observed peak RSS increased from 691.500 to 697.109 MiB (0.81%). Both
changes are reported because this is a speed improvement, not an all-resource
improvement. Concurrent workload and combined code changes limit causal
attribution of the timing difference (Table 3; Evidence E1).

In a narrower exact-output interface ablation on 7,094 fixed Morex candidates,
passing complete sequences into the native index rather than Python word lists
produced the same 4,380-family payload. Clustering time decreased from 13.611 to
12.198 s and peak child RSS from 61.203 to 58.047 MiB (Evidence E15). This was
one fixed-order run with local native build artifacts, so the magnitude is an
engineering diagnostic rather than a final performance estimate.

### Comparator conclusions depend on the measured endpoint

In a 5× read simulation from one 10-Mb development genome, TandemX, TRF and
TideHunter each recovered all 55 planted founders in both clean and high-error
reads. These conditions therefore do not establish superior family recall.
For high-error reads, eligible-array precision was 1.000 for TandemX, 0.699 for
TRF and 0.920 for TideHunter, whereas base-union precision was 0.999952,
0.999988 and 0.996805, respectively. TRF's overlapping or duplicate array
reports reduced its array precision while its base precision remained very
high. TandemX ran in 58.898 s, compared with 81.786 s for TRF and 17.137 s for
TideHunter. A different ranking for another metric is not a contradiction
(Table 2; Evidence E2).

Assembly-oriented comparisons used the independently generated 10-Mb genome
itself. TRASH2 recovered 55/55 families. Its approximate array summary matched
49/55 truth intervals under the joint interval/period criterion, but bounds
derived from its explicit native unit-to-array memberships matched 55/55, with
3.155-bp conditional boundary error. All 11,436 native unit sequences matched
the reference under standard one-based interpretation. Unit-base precision
and recall were 0.999899 and 0.999944. Thus approximate summary boundaries
cannot fairly be presented as missed repeat families (Evidence E3).

TRASH1 produced 49 regions. Primary consensus sequences recovered 43/55 families;
including its native secondary consensus sequences increased recovery to 48/55.
All 11,257 unit sequences matched the reference with an independently audited
−1-bp coordinate adjustment. The adjusted unit-base recall was 0.992706 and
precision 0.999371, despite only 26/55 matches under the combined period and
region criterion. The native fractional periodicity 171.5 was preserved rather
than rounded or omitted. These evaluations retain alternative native endpoints,
parser failures and coordinate interpretations (Evidence E4). Assembly and
HiFi-read run times are not compared across different inputs and platforms.

### Multi-k extrapolation reduces a controlled quantification bias

The quantification experiment supplied the true founder catalogue to isolate
copy estimation from discovery errors. Three independent 10-Mb genomes each
contained 54 factorial arrays and one 1.026-Mb array. Nine read conditions per
genome crossed nominal 1×, 5× and 20× sampling with three error tiers. The
estimator received observed reads, the supplied catalogue and fixed genome
size; it did not receive planted copy counts, error rates or source coordinates.

At each k, the prototype averages multiplicity-corrected counts of diagnostic
words and accounts for finite read-end opportunities. It then fits a fixed
log-linear trend across k = 15, 21, 27 and 31. At 20× with 2% unit divergence and high
read errors, the median-k21 baseline had −56.586% mean signed relative error
and 56.586% mean absolute relative error. Multi-k extrapolation gave +1.590%
signed error and 13.529% absolute error. Replacing the median alone by the
mean-k21 exposure estimate did not explain this improvement: its absolute
error remained 56.552% (Figure 1; Evidence E5).

The improvement was not universal. Across all 495 20× family conditions,
absolute error improved in 384 and increased in 111. Across the full experiment,
66 fits lacked support at one or more k values, and two additional positive
slopes were retained with model-inconsistency flags. An extrapolated word-loss
trend cannot distinguish biological divergence from sequencing error. Related
families, background word sharing and discovered-catalogue error remain
outside this conditional result.

### Joint-read uncertainty retains correlated k values and sparse-support failures

Measurements at different k values use the same reads. We therefore retained
their joint read-level contributions rather than treating the four k values
as independent replicates. Streaming first and cross moments propagated a
read-cluster sampling variance through the fitted log-copy intercept. A fixed
minimum of 20 effective reads on every k axis guarded against sparse support.
The native and Python collectors agreed with the previous point estimator,
and explicit per-read calculations independently verified the variance formula.

At 20×, 424/495 conditions yielded intervals and 401/424 contained truth (94.58%).
The fraction of all conditions with an available truth-covering interval was
81.01%. At 5×, 66 intervals were available and 60 contained truth; at 1×, only 9
were available, all covering truth. The latter conditional fraction does not
demonstrate useful low-coverage calibration. Overall 499 intervals were
available, 920 conditions had insufficient effective support and 66 lacked a
finite point fit. Dependence among families and shared read conditions means
these fractions are descriptive; genome-specific summaries are retained
(Figure 2; Evidence E6). The interval excludes extrapolation bias, genome-size
uncertainty, catalogue error and biological pooling.

### Held-out conditional tests expose failure modes and validate a frozen correction

The predeclared seeds 5101–5103 were used once after the original conditional
comparison settings were frozen. Each 199.1-kb genome contained three exact-copy
arrays, and the supplied catalogue isolated read quantification, localization
and the fixed 0.6 assembly/read-ratio rule from de novo discovery. All 177
commands completed. Across 243 positive retained-copy conditions and 162
controls, TP/FN/FP/TN were 208/35/8/154, corresponding to sensitivity 85.60%,
false-positive rate 4.94% and precision 96.30% (Evidence E11).

All 36 nonzero assembly/family localization conditions had base recall 1.0 and
minimum base precision 0.996732. Comparison errors therefore arose after this
simple exact-copy localization control. At 20× with 1% substitutions, the
50%-retention condition recovered only 2/9 positives; at 1×, each error tier
produced one false call for a fully retained assembly. These observations are
consistent with the measured error-dependent copy-number underestimation and
finite-read sampling, but do not uniquely assign cause. The aggregate metrics
do not support a robust general collapse claim, and the consumed seeds cannot
be used to tune a revised model.

We therefore calibrated a transparent rule only on earlier development seeds
4101–4103. It used k=15/21/27/31 extrapolation when available, fell back to k=21
otherwise, and changed the decision threshold from 0.6 to 0.5 only below
observed haploid depth 2. All three leave-one-genome-out development folds
selected 0.5. The rule, calibration hashes and fresh seeds 5201–5203 were
committed in 1231743; both hosted CI runs passed before one-time execution.

On the fresh 405-condition matrix, the original rule yielded TP/FN/FP/TN=
198/45/14/148 and the frozen rule yielded 208/35/12/150 (Evidence E12).
Sensitivity improved from 81.48% to 85.60%, false-positive rate from 8.64% to
7.41% and precision from 93.40% to 94.55%. Every seed gained true positives
without gaining false positives, but seed 5202 retained 12 false positives.
The 20×/1%-error/50%-retention stratum improved from 3/9 to 9/9; 1× complete-
assembly false calls decreased from five to three of 27. However, 1×/0 or 0.1%-
error/50%-retention sensitivity decreased from 6/9 to 4/9 in each stratum.
Multi-k alone had 15
unavailable rows. This supports the frozen rule under this narrow exact-copy
model, not a general or biological collapse claim (Figure 3).

We next applied the unchanged rule to a predeclared domain-shift matrix with
fresh seeds 5401–5403. The matrix crossed 1%, 3% and 5% independent
founder-to-unit substitutions with one or three same-family array segments,
while retaining the original coverage, read-error and assembly-retention tiers.
All 1,062 baseline commands completed before multi-k replay. Across 2,430 paired
family conditions, baseline TP/FN/FP/TN were 1230/228/625/347 and frozen-rule
counts were 1263/195/663/309 (Evidence E13). Sensitivity increased from 84.36%
to 86.63%, but false-positive rate increased from 64.30% to 68.21% and precision
decreased from 66.31% to 65.58%. Thus the exact-copy gain did not transfer.

Full-assembly localization mean base recall decreased from 0.754905/0.692814
for one/three segments at 1% unit divergence to 0.000793 for both structures at
3% and zero at 5%. Under 3–5% divergence, scenario-level false-positive rates
were 85.80–88.89% for the baseline and 88.89% for the frozen rule. These adverse
results identify the exact diagnostic-k-mer support filter as a primary current
failure boundary; the consumed domain-shift seeds cannot be used to tune its
replacement (Figure 4).

We therefore separated localizer development from a second fresh validation.
The first IID-proxy development run on seeds 5301–5303 retained high precision
but failed its 0.95 full-assembly recall gate (0.903608), with as many as 119
predicted fragments. A second version reused only those development seeds and
bridged exact anchors across at most one monomer length. It passed all three
development gates: full-assembly mean recall 0.975728, positive-assembly mean
precision 0.999441 and zero predicted bases in 54 absent-family rows. The first
failure was retained and the selected result hashes were frozen before any new
seed was evaluated (Evidence E16).

Commit 4e662db and its hosted Ubuntu/macOS checks preceded the one-time run of
seeds 5501–5503. All 1,062 commands completed. Full-assembly mean base recall was
0.981033, positive-assembly mean precision was 0.999565 and the absent-family
false-positive rate remained 0/54. The weakest full-assembly stratum, 5%
divergence with three segments, had mean recall 0.949274. Mean absolute relative
localized repeat-bp error across positive assemblies was 0.024861, although one
short-array condition reached 0.596721 (Figure 5; Supplementary Table S13).

The downstream assembly/read endpoint remained a trade-off. The single-k21
baseline produced TP/FN/FP/TN=875/583/12/960, with sensitivity 0.600137,
false-positive rate 0.012346 and precision 0.986471. Replaying the unchanged
multi-k/depth rule produced 1339/119/83/889: sensitivity 0.918381,
false-positive rate 0.085391 and precision 0.941632. Thus successful localization
did not make the older classifier uniformly better; its sensitivity gain added
71 false positives and lowered precision.

### Eight-species file QC and reference concordance expose normalization concerns

Complete archived FASTQ files from maize Mo17, Arabidopsis Col-0N/Col-0R/Ey15-2R, rice
Nipponbare, barley Morex, rye Lo7, wheat Chinese Spring, oat Victoria and wild
soybean YSD56 passed source checksum and full-file validation, totalling
262,731,255,175 bp in 14,937,608 reads (Table 1).
Checks covered record structure, gzip integrity, exact duplicate archive IDs,
length distributions, GC/N content and reported base-quality distributions.
These are ten included libraries from eight reported species, not eight completed
biological accuracy validations. Several accessions represent one technical batch of a
larger study. Col-0N and Ey15-2R derive from pooled plants, whereas Col-0R is
reported as a single plant. These units are not interchangeable replicates. Nested samples preserve
whole-library selection and exact read IDs. All ten sampling ladders are complete.
The YSD56 ladder contains 11.670-Mb, 110.438-Mb, 1.099-Gb and 11.050-Gb
samples, corresponding to nominal total-base coverages of 0.0116×, 0.1095×,
1.0897× and 10.9570× against the 1,008,523,555-bp assembly denominator
(Figures S1–S3; Evidence E7). These ratios are sampling denominators, not
measured nuclear depth.

Whole-library random-sample comparator diagnostics on Morex, Nipponbare,
Victoria, Chinese Spring and Lo7 successfully normalized all three read-tool
outputs. TandemX was faster than TRF and slower than TideHunter in these
observations. Its measured RSS was often lower at smaller inputs but this did
not persist across larger samples. More calls or more covered bases do not
establish greater accuracy without independent truth. The 1.129-Gb maize
experiment also showed that smaller-input rankings cannot be extrapolated
to larger inputs (Evidence E8).

The subsequently completed 110.203-Mb Victoria run provides a further resource
counterexample: TandemX used 250.453 MiB versus 176.719 MiB for TRF and 621.453 MiB
for TideHunter. Its elapsed time was 164.738 s versus 362.689 and 81.982 s,
respectively. Thus its memory advantage in the smaller oat input did not hold
against TRF in this larger nested sample.

The 126.731-Mb Chinese Spring comparison completed in 139.529 s/223.828 MiB
for TandemX, 332.287 s/311.938 MiB for TRF and 48.230 s/413.500 MiB for
TideHunter. All three methods also completed the 11.640-Mb Lo7 control.
These are concurrent development diagnostics with one repetition; their resource
observations do not replace isolated scaling experiments or biological accuracy
(Figure S4).

The 1.170-Gb Morex diagnostic completed in 1909.525 s/913.891 MiB for TandemX,
2409.145 s/286.172 MiB for TRF and 566.196 s/580.375 MiB for TideHunter.
TandemX was 20.74% faster than TRF, but used 3.19-fold its peak RSS and required
3.37-fold the TideHunter wall time; it also used 1.57-fold the TideHunter peak
RSS. This is a concrete current deficit rather than evidence of multi-metric
superiority. The scan and clustering stages remain optimization targets.

Post hoc sequence checks used actual consensus outputs and one historical query
per material. At 90% cyclic global edit similarity, all three tools recovered
CentC from the 111-Mb Mo17 input and Rice358 from the 114-Mb Nipponbare input.
TandemX's best family representatives had similarities 0.9872 and 0.9358,
respectively. In Morex, TRF and TideHunter recovered an HvT01-like 118-bp
consensus at 0.9068, as did TandemX's candidate stage, whereas the corresponding
TandemX family representative was 0.8983 and fell below the threshold. At 0.89
all methods/stages recovered it and at 0.95 none did (Evidence E10).

The three historical queries derive from materials other than the test donors.
These results therefore measure selected source-query recovery and expose a
family-representation boundary; they do not provide complete real-family recall,
false-negative rates or donor-specific truth.

The original IPK MorexV3 pseudomolecule FASTA was also acquired with its
published SHA-256 fixed in advance. The 4,296,032,540-byte file matched that
checksum; streaming QC found eight unique records and 4,225,605,719 sequence
bases, including 1,353,994 `N` bases and no other ambiguity codes (Evidence E14).
This verifies file identity and parseable reference content. Same-cultivar study
context does not establish an identical raw-read donor or exact satellite copy
truth.

Reference QC aligned the Col-0N 11.766-Mb and 118.497-Mb samples to the checked
Col-CEN v1.2 reference with its declared mitochondrial and chloroplast contigs.
In the larger sample, 7,554/7,557 reads mapped and primary query spans covered
99.9396% of input bases. Organellar primary query spans totalled 26,290,920 bp,
or 22.1871% of input, with 4,936 bp overlapping other-reference primary spans.
The chloroplast contig accumulated 166.672× summed primary aligned-base
depth, compared with 0.671–0.766× on the five other contigs. All-read nominal
depth based on the five-chromosome size was 0.9007× (Evidence E9).

This provides a concrete reason to test nuclear-depth normalization separately
from total-library exposure. It does not establish an exact 22.2% copy-number
bias or justify automatically discarding all organelle-aligned reads. Integrated
organellar sequences, uncertain placements, genotype/reference differences
and sequence-dependent sampling require additional evaluation.

In a separate 11.418-Mb Nipponbare control, 626/626 reads had primary mappings,
611 had primary MAPQ 20–254 alignments and primary query spans covered 99.9158%
of input bases. This exact GCA reference contains no organellar contigs, so zero
reported organelle mappings cannot establish absence of organellar reads. Raw
and reference BioSamples also differ; exact donor identity remains unresolved.

## Discussion

The completed experiments identify three separable improvements: exact engineering
changes reduced a measured discovery run time, a conditional multi-k model
reduced a specific copy-estimation bias, and a frozen transparent rule improved
all three assembly-comparison confusion metrics on fresh exact-copy simulations.
The independent domain-shift test then showed that this gain failed when
founder-to-unit divergence reached 3–5%, chiefly because assembly localization
lost nearly all true repeat bases and produced many false collapse calls. Bounded
anchor bridging subsequently restored high localization recall on a second
predeclared simulation, including fresh seeds, but the unchanged multi-k
classifier exchanged substantially higher sensitivity for more false positives.
None establishes that TandemX dominates existing tools. In the current read simulation, all three
methods recovered every founder. TRF retained slightly higher base precision,
and TideHunter remained faster. TRASH2 was a strong assembly baseline after
its primary unit output was interpreted correctly. Maintaining these results
is necessary for a defensible comparison and identifies where algorithmic
work is still required.

The uncertainty analysis also distinguishes usable inference from a plausible
looking numerical interval. Correlated k values cannot supply independent
replication, and a high conditional coverage rate can conceal substantial
missingness. A practical estimator must report interval availability, model
diagnostics and the sampling unit together. The present guard is conservative
at low depth, while real-genome specificity and model bias are not yet calibrated.
The observed organellar contribution adds a denominator issue that idealized
haploid simulations do not reproduce. Consequently, the experimental estimator
has not replaced the public quantification default.

The current biological evidence is incomplete. Additional species and materials
must be evaluated with documented technical/biological replication, true
genomic depth rather than file size alone, and relevant hard negatives.
Matched-donor evidence, known repeat-family recovery, collapse validation beyond
exact simulated arrays and probe/FISH concordance are necessary before biological
under-representation or experimental-success claims. Multi-gigabase files have
been acquired and sampled, but complete production-scale discovery across the
cohort and isolated cross-platform resource comparisons remain outstanding.
AI would be justified only by reproducible gains over transparent baselines on
held-out families and species; it is not a substitute for these evidence gaps.

## Methods

### Software and reproducibility

Development used the dedicated `tandemx-dev` Python 3.11 environment and a
PyO3/Rust extension. The current source passed 451 local Python tests. Its Rust
source is unchanged from commit 1231743, whose hosted Linux/macOS workflows
passed Python tests, 15 Rust tests, formatting, clippy with warnings denied and
distributable-wheel builds.
Per-run manifests and compact evidence archives preserve
the exact source, input and output hashes used for each result and take
precedence over a manuscript-level version label.
Public inputs and reuse paths are described in the repository documentation.

### Independent simulation and split discipline

Development seeds 6301–6303 generated separate 10-Mb genomes. Factorial arrays
crossed periods 61, 171 and 421 bp; copies 20, 80 and 200; target GC 0.3, 0.5 and 0.7;
and unit substitution rates 0 or 0.02. Each genome also contained a 171-bp,
6,000-copy array with 0.01 unit divergence. Background sequence was IID with
GC 0.45. Read lengths followed the observed Mo17 distribution, with uniform
circular starts and random strand. Error tiers were error-free, 0.001 each for
substitution/insertion/deletion, or 0.01/0.005/0.005, respectively. Coverage/error
conditions shared their source genome and read-start design. They are technical
conditions, not independent species. Held-out seeds 7301–7303 were not used.

The separate conditional assembly experiment used three previously reserved
seeds 5101–5103, periods 61/171/421 bp, copy counts 20/80/200 and five assemblies
retaining 100/75/50/25/0% of each exact array. Uniform circular 5-kb reads crossed
1/5/20× coverage and 0/0.1/1% substitutions. Under-representation truth was an
actual integer-copy assembly/genome ratio below 0.6; native `possible_collapse`
or `reads_only` status was positive. The unchanged matrix was run once, after
which these seeds were marked consumed. It supplies a known catalogue and omits
indels, unit divergence, ploidy and empirical sequencing bias.

Development seeds 4101–4103 then calibrated a fixed multi-k/depth decision rule.
The predeclared model and calibration artifact hashes were committed before
fresh seeds 5201–5203 were executed once. Held-out evaluation accepted only the
frozen k values, depth cutoff, thresholds, fallback rule and calibration hashes;
it did not search thresholds. Every method was scored on identical assembly,
coverage, error and family keys. Multi-k unavailable values remained explicit.

The domain-shift extension used fresh seeds 5401–5403 and the same frozen model.
It crossed unit substitution rates 0.01/0.03/0.05 with one or three array
segments separated by 500-bp independent sequence, plus all original coverage,
read-error and retention tiers. Copy variants were deterministic prefixes across
assembly fractions. Multiple non-overlapping truth intervals were summed by
family for copy number and assembly ratio and union-scored for localization.
Substitutions were independent and length preserving; the design does not model
indels, empirical satellite evolution or contig breaks.

Localizer development used seeds 5301–5303 only. The IID proxy transformed an
observed exact k-mer fraction `f` to `f^(1/k)` under an explicit independent-
substitution approximation and applied a minimum proxy identity of 0.90. After
the first version failed the predeclared full-assembly mean-recall gate, the
second version allowed exact-anchor gaps up to the larger of `2k` and one monomer
length. Its selection gates were full-assembly mean recall at least 0.95,
positive-assembly mean precision at least 0.95 and absent-family false-positive
rate zero. Development validation, metrics, environment and configuration hashes
were embedded in `abundance_localizer_heldout_v1.json`; the runner re-hashed them
and recomputed the gates before creating held-out output. The configuration and
guard were committed and passed hosted CI before seeds 5501–5503 were executed
once. Held-out reporting retained the same six divergence/fragmentation
scenarios and all original coverage, read-error and assembly-retention tiers.
Exact-k-mer IID identity is not alignment identity; overlapping k-mers,
insertions/deletions, shared-word filtering and structured satellite variation
violate the approximation.

### Comparator execution and scoring

Read comparators used identical FASTA input, one thread, period scope 30–1000 bp
and minimum repeat span 100 bp. Truth fragments, array matches, family recovery
and interval unions were evaluated separately. Array matching was one-to-one,
requiring interval IoU ≥ 0.5 and period error ≤ max(2 bp, rounded 2% truth period).
Cyclic sequence recovery used the documented exact edit-identity threshold 0.9,
with reverse-complement treatment and native duplicate handling retained.
Failure, timeout or malformed output was unavailable evidence, not zero recall.

TRASH1 and TRASH2 ran de novo without supplied templates in an immutable
ARM64 Linux container with R 4.4.3, one CPU and 8-GiB memory allocation. Native
commits were f7d53a0 and f290a5e. Monomer sequences were checked by strand-aware
reference extraction under explicit coordinate alternatives. TRASH2 unit-array
membership was taken from its own identifiers; truth overlap did not define
the predicted boundaries. Native primary/secondary consensus alternatives and
fractional periodicity were preserved. Container time/RSS, cgroup memory and
available cgroup CPU counters were distinguished from host-client measurements.

### Conditional multi-k estimator and read-level variance

For family f and k, the target bank contains catalogue-exclusive canonical
circular words after the documented low-complexity filter. A word's contribution
is corrected for its multiplicity in the founder. For read i, X_ik=max(0,L_i−k+1)
and Y_ifk is the mean corrected diagnostic-word count. The point estimate is
theta_fk=G Σ_iY_ifk/Σ_iX_ik. OLS fits log(theta_fk)=a_f+b_f * k, with estimated
copies exp(a_f). There are no truth-derived pseudocounts or slope clamps.

The intercept weights are w_k=1/K−mean(k)(k−mean(k))/Σ(k−mean(k))².
Read influences are g_if=Σ_k w_k(Y_ifk/Σ_iY_ifk−X_ik/Σ_iX_ik).
Variance is N/(N−1)Σ_i g_if², and the approximate sampling interval is
exp(a_f±1.959964√variance). Streaming cross moments retain all k correlations.
Minimum effective support is min_k[(Σ_iY_ifk)²/Σ_iY_ifk²]. Values below 20,
missing diagnostic support and numerical-resolution failures yield unavailable
intervals. The genome size and catalogue are fixed conditioning inputs.

### Raw-data and reference QC

Complete ENA files were checked against official source sizes/MD5 values and
independently hashed. A streaming FASTQ parser checked all records through EOF;
an on-disk exact ID index detected duplicate archive identifiers. Phred+33
summaries describe reported quality, not empirically measured error. Fixed-seed
hash selection scanned complete libraries and retained nested subset IDs,
achieved sizes and distribution histograms. No prefix-download sample was
substituted for this whole-library design.

Reference concordance used minimap2 2.31-r1302 with `map-hifi -c --cs=short`,
secondary output retained with cap 5, one thread, 50-Mb query batches and a
single 3-Gb index batch [6]. Reference hashes and all query/target lengths were
validated. SQLite sorting supported query-span unions and CIGAR target-block
coverage without loading the complete read collection. `P/I` tags were treated
as primary and `S/i` as secondary; MAPQ 255 was missing. Target M/= /X blocks
contributed to aligned-base coverage, excluding deletions/skips. Summed
alignment-base depth may include overlapping split alignments from one read
and is not unique-copy depth. Compartment spans were separately unioned with
their overlap retained. Unmapped reads were not classified as contamination.

## Data and code availability

Code and committed evidence are at [CongyangY/TandemX](https://github.com/CongyangY/TandemX).
Accession-level source receipts, checksums, environments, command lines and
compact native outputs are supplied in the evidence directories below. Large
FASTQ/FASTA/PAF files and scratch databases remain under the configured data
root and are regenerated from recorded public accessions. A clean versioned
release, complete archive deposition and all manuscript source-data packaging
remain release requirements; a current development checkout is not a final
publication release.

## Figure legends

**Figure 1. Conditional multi-k copy estimation across simulated divergence and
read errors.** Six-panel development figure and full legend:
`evidence/factorial_multik_replay/figures/factorial_quantification.pdf` and
`evidence/factorial_multik_replay/figure_legend.md`. All paired source rows are
retained. Gains and worse-performing conditions are shown together.

**Figure 2. Joint-read sampling intervals and explicit unavailable estimates.**
Six panels show availability, seed-specific conditional coverage, divergence/
error strata, 20× point estimates, interval width against effective support and
missingness by abundance. Dashed lines in coverage panels denote nominal 95%.
The exact panel legend is in `evidence/factorial_joint_multik/README.md`.
The inspected version 2 layout moves the panel-A legend clear of the data;
its source rows are identical to version 1.

**Figure 3. Predeclared conditional assembly-comparison validation.** Six panels
show the fresh three-genome design, overall sensitivity/false-positive rate/
precision, confusion counts, per-seed sensitivity and false-positive rates, and
baseline-to-frozen-rule sensitivity across nine 50%-retention coverage/error
strata. The heatmap includes the two adverse 1× cells. Inputs, panel values and
hashes are in `evidence/abundance_multik_collapse_heldout/figures_v2`. These are
exact-copy known-catalogue simulations, not biological collapse truth.

**Figure 4. Frozen-rule failure under unit divergence and interrupted arrays.**
Six panels show the predeclared domain-shift design, aggregate classification
metrics, confusion counts, sensitivity and false-positive rate by unit
divergence and structure, and full-assembly localization recall. The frozen
rule gains sensitivity while increasing false positives and lowering precision;
at 3–5% divergence, localization recall approaches or reaches zero. Inputs,
panel values and hashes are in
`evidence/abundance_domain_shift_multik_heldout/figures_v3`. This adverse
known-catalogue simulation is not biological collapse truth.

**Figure 5. Divergence-aware anchor bridging and fresh held-out validation.**
Six panels show the paired failed/selected development runs, full-assembly recall,
fragment-count error, held-out recall and repeat-bp error across assembly
retention, and downstream classification metrics. The localizer passes all three
predeclared aggregate gates on seeds 5501–5503. The unchanged frozen multi-k
classifier gains sensitivity while increasing false-positive rate and lowering
precision. Inputs, panel values, hashes and the complete legend are in
`evidence/abundance_localizer_multik_heldout/figures_v3` and
`evidence/abundance_localizer_multik_heldout/figure_legend.md`. This known-
catalogue substitution simulation is not biological collapse truth.

**Figure S1. Complete Mo17 input QC.** Four-panel source-backed distributions,
with input and plotting receipts, in `evidence/Mo17_input_qc/figures_checked`.

**Figure S2. Complete Col-0N input QC.** Four-panel source-backed distributions,
with pooled-material and reported-quality limits, in
`evidence/Col0N_input_qc/figures`. These plots do not establish biological truth.

**Figure S3. Cross-cohort input QC.** Six panels show validated sequence volume,
median and N50 read length, whole-file GC fraction, and row-normalized read
length, per-read GC and reported mean-quality distributions for ten complete
libraries from eight reported plant species. The final bins include values at or
beyond the labelled bound. Reported quality is not empirical accuracy, and the
libraries are not interchangeable biological replicates. Source rows and hashes
are in `evidence/multispecies_input_qc/figures_v2`.

**Figure S4. One-thread real-read comparator diagnostics.** Six panels show wall
time, throughput, peak RSS, paired wall-time and RSS ratios to TRF, and called-
base fraction across ten nested inputs from five materials. Lines connect nested
sizes only within a material/tool. These single executions overlapped acquisition
work and are not final isolated rankings. Called-base fraction is output extent,
not accuracy. The inspected version 2 and source rows are in
`evidence/multispecies_real_diagnostics/figures_v2`.

Additional method, biological validation and resource-scaling multi-panel
figures remain required; their absence is tracked in `submission_readiness.md`.

## Tables and evidence index

- Table 1: ten-library, eight-species input/QC table in `tables/input_cohort.tsv`.
- Table 2: executed six-run read comparison in
  `evidence/factorial_discovery_s6301_5x/summary.tsv`; scope/denominators in its README.
- Table 3: exact 1.129-Gb replay in `evidence/Mo17_alignment_workspace/full_1129Mb/validation.json`
  and the preceding run archived under `evidence/Mo17_clustering_index`.
- Supplementary Table S1: all joint-interval family conditions in
  `evidence/factorial_joint_multik/metrics.tsv`; S2: per-k inputs in its `per_k.tsv`.
- Supplementary Table S3: all conditional point-estimator comparisons in
  `evidence/factorial_multik_replay/metrics.tsv`.
- Supplementary Table S4: native TRASH/TRASH2 region, catalogue, unit and
  coordinate-audit tables in the corresponding evidence directories.
- Supplementary Table S5: per-read, composition and reference-mapping summaries
  in `evidence/reference_mapping_diagnostics`.
- Supplementary Table S6: selected historical-query recovery, threshold
  sensitivity and actual consensus provenance in `evidence/known_query_recovery`.
- Supplementary Table S7: cross-cohort QC and real-diagnostic panel sources in
  `evidence/multispecies_input_qc/figures_v2` and
  `evidence/multispecies_real_diagnostics/figures_v2`.
- Supplementary Table S8: all held-out conditional copy-number, localization,
  comparison and per-command resource rows in `evidence/abundance_heldout`.
- Supplementary Table S9: fresh paired baseline/frozen-rule comparison rows in
  `evidence/abundance_multik_collapse_heldout`.
- Supplementary Table S10: all domain-shift baseline, localization, paired
  frozen-rule and panel-source rows in `evidence/abundance_domain_shift_heldout_baseline`
  and `evidence/abundance_domain_shift_multik_heldout`.
- Supplementary Table S11: MorexV3 contig lengths and base counts in
  `evidence/MorexV3_reference_qc/reference_receipt.json`.
- Supplementary Table S12: Morex sequence-native interface resource and parity
  receipts in `evidence/Morex_115Mb_index_interface_v2`.
- Supplementary Table S13: failed/selected localizer development, fresh held-out
  localization, copy-number/comparison, resource and paired classifier rows in
  `evidence/abundance_localizer_development_v1`,
  `evidence/abundance_localizer_development_v2`,
  `evidence/abundance_localizer_heldout_baseline` and
  `evidence/abundance_localizer_multik_heldout`.

E1: `Mo17_alignment_workspace`; E2: `factorial_discovery_s6301_5x`;
E3: `TRASH2_factorial_s6301`; E4: `TRASH_factorial_s6301`;
E5: `factorial_multik_replay`; E6: `factorial_joint_multik`;
E7: ten `*_input_qc` directories and their source manifests;
E8: `multispecies_real_diagnostics`; E9: `reference_mapping_diagnostics`;
E10: `known_query_recovery`; E11: `abundance_heldout`;
E12: paired `abundance_heldout_v2` and `abundance_multik_collapse_heldout`;
E13: paired `abundance_domain_shift_heldout_baseline` and
`abundance_domain_shift_multik_heldout`; E14: `MorexV3_reference_qc`.
E15: `Morex_115Mb_index_interface_v2`; E16: paired localizer development and
held-out directories listed for Supplementary Table S13.
These are authoritative result locations, not replacements for the remaining
final table/figure packaging and journal-specific formatting checks.

## References

1. Benson G. Tandem repeats finder: a program to analyze DNA sequences.
   Nucleic Acids Research 27,573–580 (1999).
   [doi:10.1093/nar/27.2.573](https://doi.org/10.1093/nar/27.2.573).
2. TideHunter. [doi:10.1093/bioinformatics/btz376](https://doi.org/10.1093/bioinformatics/btz376).
3. TRASH. [doi:10.1093/bioinformatics/btad308](https://doi.org/10.1093/bioinformatics/btad308).
4. Satellite Repeat Finder: primary publication and source.
   [Genome Research, doi:10.1101/gr.278005.123](https://doi.org/10.1101/gr.278005.123);
   [source](https://github.com/lh3/srf).
5. A complete assembly of the rice Nipponbare reference genome.
   [doi:10.1016/j.molp.2023.08.003](https://doi.org/10.1016/j.molp.2023.08.003).
6. Minimap2 author documentation and source.
   [manual](https://lh3.github.io/minimap2/minimap2.html);
   [source](https://github.com/lh3/minimap2).

Bibliographic metadata and the full primary dataset references require final
reference-manager curation. No author, funding, conflict or accession-deposition
details are invented in this draft.
