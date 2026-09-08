# TandemX: read-level evidence and sampling uncertainty for plant satellite repeat analysis

**Evidence-backed development manuscript, 7 September 2026.** Author information
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
and 1,485 family conditions. On three additional genomes, a predeclared control-
depth gate reduced aggregate mean absolute relative copy-number error from
0.4088 to 0.3632 and passed all nine frozen simulation gates; 661/1,485 family
conditions improved, 495 were unchanged and 329 worsened. In the 20× stratum
with 2% unit divergence and high
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
from 98.65% to 94.16%. A subsequent seed-robust blend was frozen before three
additional genomes were evaluated. It increased sensitivity from 62.83% to
71.60%, but false-positive rate increased from 2.98% to 4.63% and precision
decreased from 96.93% to 95.87%, failing the predeclared held-out gates.
Exact native optimizations shortened a 1.129-Gb maize replay by 23.61% relative
to its preceding indexed version while preserving seven output files byte for
byte; peak memory increased 0.81%. External comparisons showed task-dependent
trade-offs. After a failed first cascade promotion and a guarded development
revision, a frozen 96-run validation passed all 14 gates: the
TandemX/TideHunter wall-time and direct-child peak-RSS geometric-mean ratios
were 1.9779 and 0.3969, minimum positive array recall/precision were 1.0 and
the negative-control call rate was zero. Some scenario runtime ratios exceeded
4, so this does not establish universal speed, memory or accuracy superiority.
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

The family audit also emits a directional candidate architecture graph instead
of collapsing all related representative lengths. Related 171/342/684-style
pairs retain shorter-to-longer edges, their nearest integer length multiple and
the underlying identity, overlap and shared-word evidence. Near-integer ratios
are labelled putative period multiples; non-integer related pairs remain
unresolved. All alternatives remain in the graph. This makes harmonic and
possible higher-order relationships machine-readable without selecting a rooted
tree or treating a length multiple as a validated higher-order repeat.

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

We then reduced the dominant elastic-alignment storage and call overhead without
changing the recurrence: traceback directions were packed from one byte to two
bits per cell, and all candidate periods for one read reused one native uppercase
buffer. Complete Mo17 discovery replays preserved all six historical core
products at 11.681 Mb and all seven products at 111.506 Mb byte for byte. Elapsed
time decreased from 17.348 to 12.361 s and from 134.684 to 100.039 s,
respectively; peak RSS decreased from 81.703 to 62.484 MiB and from 178.438 to
154.766 MiB. Each comparison used one historical baseline and one replay, so the
magnitudes remain engineering diagnostics rather than repeated performance
estimates (Evidence E19). Existing external-tool diagnostics still place
TideHunter ahead of TandemX in elapsed time on these real-data scales.

The staged native-screen cascade was next evaluated once on predeclared seeds
3101--3103. The complete matrix contained 48 scenario-seed datasets and three
technical repetitions per tool. TandemX completed all 144 runs, its minimum
positive array recall and precision were both 0.985714, no negative-control read
was called, related-family cyclic monomer recall was 1.0, and the direct-child
peak-RSS geometric-mean ratio to TideHunter was 0.455130. Promotion nevertheless
failed. All nine TRF low-complexity-control runs reached the frozen 180-s timeout,
violating the required comparator-completion gate, and TandemX/TideHunter wall-
time ratio was 2.457943 against the predeclared maximum 2.0. Ten of 12 gates
passed, but the mode remains non-default. The failed processes, gates and both
favorable and unfavorable metrics are retained together (Figure 8; Evidence
E20).

We next used fresh development seed 1201 to diagnose the runtime deficit without
reusing those held-out rows. A broad 30%-read-span/95%-shifted-identity gap-free
path reduced the paired runtime ratio from 2.357698 to 1.933332, but worsened
0.1%-indel boundary MAE from 1.364 to 14.221 bp and was rejected. A read-level
audit then added two observable guards: at least 95% valid shifted columns and
at most 2% residual from an integer number of proposed units. The resulting
candidate completed the same 96-run development matrix with a ratio of 1.989148,
minimum base-union F1 0.997597 and maximum boundary MAE 2.921 bp. Commit
`054b935`, the exact rule, seed 2201 and 14 gates were fixed before both
working-branch and main hosted workflows passed.

Seed 2201 was then run once. All 96 TandemX/TideHunter executions completed and
all 14 gates passed. TandemX's minimum positive array recall and precision were
1.0, minimum base-union F1 was 0.997445, maximum positive boundary MAE was
2.35 bp and no negative-control read was called. Relative to TideHunter, the
wall-time geometric-mean ratio was 1.977877 and direct-child peak-RSS ratio was
0.396864; the worst positive recall and precision differences were both zero,
satisfying non-inferiority. The wall-time result is an aggregate across the
frozen distribution: individual scenario ratios reached 4.203, and technical
repetitions are not biological replicates. The validation therefore resolves
the predeclared synthetic promotion gate without establishing real-data or
per-condition dominance (Figure 11; Evidence E25).

A subsequent prioritization audit reused these fixed outputs without rerunning
or retuning the benchmark. At 0.1%, 1% and 4% total indels, both TandemX and
TideHunter had array recall 1.0. TandemX precision was 1.0 in all three
conditions, whereas TideHunter precision was 1.0, 1.0 and 0.958904. TandemX
boundary MAE was 1.229--1.621 bp, but its median elapsed time was 3.34--3.80
times that of TideHunter and its measured direct-child peak RSS was 0.32--0.43
times TideHunter's. We therefore did not introduce a result-driven
seed-and-chain replacement: on this synthetic distribution the remaining clear
deficit was elapsed time rather than read-local accuracy. This decision does
not resolve real-read interval truth or indel profiles outside the generator
(Evidence E27).

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

### Public depth calibration reduces aggregate bias but exposes low-depth and interval failures

We then evaluated the public single-k `quantify` command on the same three
independent 10-Mb development genomes. Each of 27 read conditions was run with
total-bases depth, the planted aggregate error rate, 1,000 simulation-truth-
assisted whole-genome-unique control k-mers, and controls plus the planted error
rate. All 108 commands completed, yielding 1,485 family conditions per method.
Mean absolute relative error was 0.401898 for total-bases depth, 0.359640 for the
oracle error correction and 0.365728 for empirical controls. Mean signed errors
were -0.262533, -0.145827 and -0.113112, respectively. Controls improved the
mean absolute error in every independent genome but were worse at nominal 1x;
969/1,485 paired family conditions improved and 450 worsened.

Applying the same global survival factor to repeat and control counts cancelled
algebraically: controls plus oracle produced exactly the same 1,485 copy-number
estimates as controls alone, while median wall time increased from 1.771 to
4.934 s in the original implementation. Replacing the Python per-base FASTA
survival scan with its exact all-ACGT formula retained byte-identical 5,940-row
metrics and all 108 copy-number products in a full replay. Driver time changed
from 610.507 to 349.217 s; oracle-error and controls-plus-oracle median times
decreased 58.26% and 55.85%, respectively. This single same-machine replay is an
engineering check rather than a publication timing distribution (Evidence E23).
A post-hoc development rule that used controls only when their mean depth was at
least 2 reduced aggregate error to 0.356117 and improved all three development
genome means. We froze that threshold, the same generating process, three new
seeds and nine gates in commit `9eda196`; Ubuntu and macOS workflows passed
before any new data were generated. The three validation datasets then passed an
independent audit of 30 manifests and all 93 declared payloads. All 54 public
commands completed. The frozen rule reduced MARE from 0.408767 to 0.363222,
with reductions of 0.037855, 0.043531 and 0.055248 in seeds 6401--6403, and
passed all nine gates (Figure 10; Evidence E24). Of 1,485 paired family
conditions, 661 improved, 495 were unchanged and 329 worsened. The rule retained
the total-bases estimate for all nine 1× conditions and used controls for all 18
5×/20× conditions. Ungated controls reached a slightly lower validation MARE of
0.359443; this adverse contrast was retained and the threshold was not refit.
The validated rule is now an explicit `--single-copy-min-depth 2` option, while
the controls-only behavior remains unchanged by default because control selection
still used simulation truth. Across every development method and coverage
stratum, truth inclusion in the exported diagnostic 10th--90th spread was far
below 0.95, confirming that these endpoints cannot be interpreted as a sampling
confidence interval (Figure 9; Evidence E22).

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

We next separated classifier development from the consumed localizer validation.
The predeclared seeds 5601–5603 used the same six divergence/fragmentation
scenarios and produced 2,430 paired family conditions. Single k=21 yielded
TP/FN/FP/TN=873/585/10/962. Among 20 log-space single/multi-k blends and
decision thresholds, the pooled selector chose alpha 0.75 and threshold 0.45,
yielding 1042/416/8/964. Sensitivity increased from 0.598765 to 0.714678, FPR
decreased from 0.010288 to 0.008230 and precision increased from 0.988675 to
0.992381. However, all three leave-one-seed-out training folds selected
different parameters. Their combined held-back predictions had FPR 0.021605
and precision 0.980374. The predeclared stability, FPR and precision gates
therefore failed; this is a development failure, not held-out validation.

A second development selector was frozen after this failure. It requires each
development genome separately to preserve baseline FPR and precision, then
maximizes the minimum genome-level sensitivity gain. Its inputs and source
hashes are fixed. It selected alpha 0.5 and threshold 0.5, yielding
TP/FN/FP/TN=1002/456/7/965. Relative to k=21, the worst genome-level sensitivity
gain was 0.076132, the maximum FPR delta was zero and the minimum precision delta
was 0.001558. Full-development deltas were +0.088477 sensitivity, -0.003086 FPR
and +0.004387 precision. All six development gates passed. The localizer and
classifier parameters and evidence hashes were then embedded in a held-out
configuration while seeds 5701–5703 remained untouched. This post-v1
refinement cannot itself establish generalization.

Commit 300e48d and its hosted Ubuntu/macOS checks passed before seeds 5701–5703
were used once. Across 2,430 paired family conditions, single k=21 produced
TP/FN/FP/TN=916/542/29/943 and the frozen blend produced 1044/414/45/927
(Figure 6; Evidence E17). Sensitivity increased from 0.628258 to 0.716049, but
false-positive rate increased from 0.029835 to 0.046296 and precision decreased
from 0.969312 to 0.958678. Seed 5703 had the largest adverse changes: FPR rose
by 0.040123 and precision fell by 0.025985. The overall and worst-seed FPR and
precision gates therefore failed. All 16 additional false positives occurred at
nominal 1× coverage; at 5× and 20× both methods had zero false positives while
the blend improved sensitivity by 0.094650 and 0.082305. This stratification is
a post-hoc diagnosis on consumed held-out data, not evidence that a revised rule
has passed independent validation.

The same 5701–5703 baseline retained the frozen localizer behavior at aggregate
level. Full-assembly mean base recall was 0.973858, positive-assembly mean base
precision was 0.999631 and no predicted bases occurred in 54 absent-family rows.
The 5%-divergent three-segment stratum had mean recall 0.937682, below the 0.95
aggregate gate used for the whole matrix. Successful localization in aggregate
therefore did not rescue classifier calibration at low read depth.

Because all additional false positives in that failed test occurred at nominal
1×, we treated 5701–5703 as consumed development evidence and specified a
transparent depth gate before generating another split. The rule retains single
k=21 with threshold 0.6 below estimated haploid depth 2 and uses the alpha-0.5,
threshold-0.5 blend otherwise. Across the six consumed development genomes it
increased sensitivity by 0.060700 with no FPR increase and a 0.001883 precision
increase. The rule, five development-artifact hashes, unchanged localizer and
acceptance gates were frozen in commit 62892a6; hosted Ubuntu and macOS checks
passed before seeds 5801–5803 were generated once.

The fresh held-out matrix contained 2,430 paired family conditions (Figure 7;
Evidence E18). Single k=21 produced TP/FN/FP/TN=867/591/16/956, whereas the
depth-gated rule produced 953/505/16/956. Sensitivity increased from 0.594650 to
0.653635, FPR remained 0.016461 and precision increased from 0.981880 to
0.983488. Sensitivity improved in every seed (minimum delta 0.030864), the
maximum seed-level FPR delta was zero and the minimum seed-level precision delta
was zero. All predeclared overall and seed-level gates passed. The 810
estimated-depth-below-2 rows were unchanged; sensitivity gains were 0.072016 at
nominal 5× and 0.104938 at 20×. No standard-depth multi-k fallback occurred.
The same run gave full-assembly mean localization recall 0.980059,
positive-assembly mean precision 0.999617 and zero predicted bases in 54 absent-
family rows. These results establish rule transfer within the frozen known-
catalogue IID-substitution matrix, not biological collapse accuracy or general
superiority over external software.

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

We separately audited four retrospective old-to-new assembly candidates using
exact ENA run metadata. Twelve PacBio genomic-WGS runs were retained across
Arabidopsis, maize, rice and soybean. The initial audit correctly found that the
rice HiFi run and AGIS1.0 assembly use different BioSamples and that the seven
selected Mo17 CCS runs use BioSamples different from the T2T assembly. It did
not, however, enroll the stronger Ey15-2 comparison in the Arabidopsis source
paper. Rabanal et al. explicitly compare CLR and HiFi assemblies of the same
Ey15-2 sample (9994/CS76399), and the archived HiFi run `ERR8666125` is assigned
to BioSample `SAMEA13018399` [7]. This meets our predeclared Tier A source rule.

Before inspecting TandemX old/new localization differences, we froze a primary
comparison of Bionano-scaffolded `9994.CLR_Canu` against
`9994.HiFi_Hifiasm`, with the authors' final HiFi-Hifiasm plus CLR-Canu assembly
reserved for sensitivity. Primary source eligibility requires at least 15 kb
localized in the newer assembly and deliberately does not depend on new/read
agreement. The analysis preserves every family fate and fixes 5-, 15- and
50-kb denominator sensitivities. Acquisition and execution remain incomplete,
so no donor-matched performance result is reported here yet. The newer assembly
also shares HiFi evidence with TandemX; it is a donor-matched high-quality
reference proxy rather than absolute independent copy truth (Evidence E21).

We screened cassava TME204 as a second donor-matched comparison because its CLR
and HiFi reads were generated from the same DNA sample and the published
CLR-Falcon/Falcon-Unzip assembly was less haplotype-resolved than the HiFi-
hifiasm assembly [10]. However, the inspected 86-file GigaDB inventory and
linked Mendeley record expose the final HiFi haplotypes but not the historical
CLR assembly sequence. We therefore assigned
`metadata_blocked_no_public_old_assembly`, downloaded none of the approximately
56.14-GB compressed read set and did not replace the missing historical
assembly with a new reassembly. This retained failure also identifies a design
requirement for any future TME204 analysis: its two haplotypes need a frozen
diploid scoring model before they can be used as reference evidence (Evidence
E28).

We also screened the newer Heinz 1706 SL5.0-to-SL-T2T contrast. Both assemblies,
the exact SL5.0 HiFi run and the official 831,451,202-byte SL-T2T FASTA are
public. However, the SL-T2T study explicitly allows mixed-seed heterozygosity
and sample differences between the ONT and HiFi sequencing [11]. We therefore
retained `not_source_eligible_donor_mismatch_risk`, downloaded none of the
28,768,190,557-byte HiFi archive and did not score a donor-matched result. The
pair remains eligible only for a separately frozen same-cultivar descriptive
analysis (Evidence E29).

B73-Ab10 provided a stronger paper-level same-DNA candidate: the v2 study
states that its HiFi library reused the HMW DNA from the CLR-based v1 assembly
[12]. We did not promote it, because its newer assembly is reported to retain
N-gaps predominantly within tandem-repeat arrays, with only about 28% of the
approximately 30.67-Mb Ab10 knob assembled and slightly smaller assembled
knob180/TR-1 totals than v1. In addition, the archival BioSample identifiers
differ and the public 257,275,342,048-byte sequence deposit exposes split
subreads rather than a ready CCS-read product. This candidate therefore failed
the reference-truth gate as `reference_ineligible_new_assembly_incomplete`; no
large download or accuracy run was made (Evidence E30).

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

We then ran pinned TideCluster 1.21.2 on deterministic nested 10- and 100-Mb
sets of 1-Mb MorexV3 windows. Both external TideHunter and clustering stages
completed. The 10-Mb set produced 87 final intervals in 28 operational families
covering 239,939 union bp (2.39939%); the 100-Mb set produced 1,372 intervals in
130 families covering 2,356,743 bp (2.356743%). Positive calls occurred in
10/10 and 99/100 windows. TideHunter stage time/RSS changed from 25.03 s/561,556
kB to 331.52 s/5,152,444 kB, whereas clustering changed from 45.92 s/7,809,052
kB to 60.99 s/7,770,936 kB. Thus clustering already required about 7.4 GiB at
both sampled sizes, and the 1-Gb run was withheld pending a safer resource plan.
TideCluster overlap resolution also made final intervals differ from raw
TideHunter intervals: at 100 Mb, 1,200 final intervals retained exact
intermediate coordinates, 41 were clipped and 131 merged multiple same-family
intervals. A family-consensus membership map and complete overlap coverage were
therefore required to retain sequence provenance; copy number remained
unavailable for 202 merged/resolved intervals (Figure S5; Evidence E26).

These nested windows quantify execution behavior and expose a normalization
requirement. They are neither whole-chromosome runs nor independent accuracy
tests, because no curated MorexV3 family/array denominator exists and windowing
removes long-range chromosome context.

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
Separate classifier-development data then showed that a pooled blend could
improve aggregate sensitivity, FPR and precision while failing seed-level
stability. This prompted a minimax development selector that treats every seed
as a guardrail. Its independent test failed despite a sensitivity gain: all
additional false positives occurred at 1×, and the held-out FPR and precision
regressed overall and in the worst seed. This shows that seed-level robustness
within three development genomes was insufficient to guarantee transfer.
The guarded cascade subsequently passed its separate frozen synthetic
validation, correcting the earlier aggregate runtime-gate failure while using
about 40% of TideHunter's measured direct-child RSS. The 1.977877 runtime ratio
was close to the predeclared limit and several scenarios remained more than
four-fold slower. It therefore supports the cascade routing rule on the tested
distribution rather than a blanket superiority claim. In another read
simulation all three methods recovered every founder; TRF retained slightly
higher base precision, and existing real-read diagnostics still place
TideHunter ahead in elapsed time. TRASH2 was a strong assembly baseline after
its primary unit output was interpreted correctly. Maintaining these results
is necessary for a defensible comparison and identifies where algorithmic work
is still required.

The uncertainty analysis also distinguishes usable inference from a plausible
looking numerical interval. Correlated k values cannot supply independent
replication, and a high conditional coverage rate can conceal substantial
missingness. A practical estimator must report interval availability, model
diagnostics and the sampling unit together. The validated depth guard is
conservative at low depth, while real-genome control specificity and model bias
are not yet calibrated. In its validation split, controls alone had slightly
lower aggregate error than the gate; the gate's support is therefore robustness
against total-bases normalization rather than superiority to every internal
alternative. The observed organellar contribution adds a denominator issue that
idealized haploid simulations do not reproduce. Consequently, the experimental
estimator is exposed as an opt-in public mode and has not replaced the controls
default.

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
PyO3/Rust extension. Frozen validation commit `054b935` passed 531 local Python
tests, compileall, Rust formatting and release clippy with warnings denied.
Hosted runs `34095006477` and `34095171697` passed Linux/macOS Python tests,
executable Rust checks and distributable-wheel builds before seed 2201 was used.
The validation source snapshot included and hashed two untracked local native
extensions, which triggered its conservative revision warning; scoped tracked
source matched commit `054b935` before execution.
Per-run manifests and compact evidence archives preserve
the exact source, input and output hashes used for each result and take
precedence over a manuscript-level version label.
Public inputs and reuse paths are described in the repository documentation.

After sequence clustering, TandemX streams an all-pair or exactly gated
related-pair catalogue audit. For relationships classified as possible
higher-order or partial, it directs an architecture edge from the shorter to the
longer representative. The nearest integer multiple is accepted as a candidate
period multiple only when it is at least two and differs from the observed
length ratio by no more than 0.05. The output retains local identity, overlap,
shared canonical-word fraction and orientation. These are heuristic catalogue
relationships; the method does not infer unit order, choose a biological root,
or distinguish true HORs from harmonic calls and partial representatives.

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

Classifier development used seeds 5601–5603 and reserved 5701–5703. For each
paired family condition, candidate read copy number was the log-space blend
`exp((1-alpha) log(C21) + alpha log(Cmulti))`, with fallback to k=21 when a
multi-k estimate was unavailable or either input was nonpositive. Alpha values
0, 0.25, 0.5, 0.75 and 1 crossed decision thresholds 0.45, 0.5, 0.55 and 0.6;
alpha 0 with threshold 0.6 was the exact baseline anchor. Version 1 selected
maximum sensitivity subject to pooled FPR and precision constraints and tested
parameter stability in leave-one-seed-out folds. After that gate failed,
version 2 retained the same candidates but required the constraints within each
seed and ranked eligible candidates by minimum seed-level sensitivity gain,
then full sensitivity, worst seed FPR delta, worst seed precision delta, alpha
and threshold. Version 2 is explicitly development refinement and must be
frozen before any reserved seed is evaluated. The selected alpha 0.5 and
threshold 0.5, localizer hashes and both development gates were embedded in the
held-out configuration. The runner re-hashed and recomputed both gates before
creating any output, and commit 300e48d passed hosted Ubuntu/macOS checks before
seeds 5701–5703 were used. Held-out acceptance required an overall sensitivity
gain of at least 0.05, a gain of at least 0.03 in every seed, and no FPR or
precision regression overall or in the worst seed. No parameter was selected or
refitted on held-out rows. Coverage-specific results reported after the failed
gate are diagnostic and make these seeds development data for any future rule.

Version 3 used that diagnostic only to define a new development rule: retain the
single-k21 estimate and threshold 0.6 when the single-k21 estimated haploid depth
was below 2; otherwise apply the already evaluated alpha-0.5 blend and threshold
0.5. Development comprised all six consumed seeds 5601–5703, while 5801–5803
were newly reserved. Nine full/cohort/seed development deltas, five development-
artifact hashes and the unchanged localizer were embedded in
`abundance_classifier_depth_gated_heldout_v1.json`. The runner checked those
artifacts and seed disjointness before output creation. The multi-k stage emitted
only raw paired single/multi-k rows. The final evaluator required one-to-one
method pairs at every family-condition key, applied the frozen depth rule, and
recorded that held-out fitting and selection were absent. The same six held-out
criteria used above were retained. Commit 62892a6 passed hosted Ubuntu/macOS
checks before seeds 5801–5803 were executed once.

### Comparator execution and scoring

Read comparators used identical FASTA input, one thread, period scope 30–1000 bp
and minimum repeat span 100 bp. Truth fragments, array matches, family recovery
and interval unions were evaluated separately. Array matching was one-to-one,
requiring interval IoU ≥ 0.5 and period error ≤ max(2 bp, rounded 2% truth period).
Cyclic sequence recovery used the documented exact edit-identity threshold 0.9,
with reverse-complement treatment and native duplicate handling retained.
Failure, timeout or malformed output was unavailable evidence, not zero recall.

For cascade speed refinement, seed 1201 was used only for development. The
screen audit exposed read length, proposed interval, shifted identity, valid
A/C/G/T pair fraction, composition-adjusted identity and residual from an
integer number of units; simulation truth scored candidate rules but was not an
input to the final routing decision. A gap-free proposal was accepted only when
its interval spanned at least 30% of the read, shifted identity was at least
0.95, valid-pair fraction was at least 0.95, unit-span residual was at most 0.02
and composition-adjusted identity was at least 0.7. Other reads retained the
elastic path. The validation configuration froze seed 2201, three repetitions,
16 scenarios and 14 process, determinism, array, base-union, boundary,
negative-control, related-family, TideHunter non-inferiority, runtime and RSS
gates. Source/config commit and hosted-CI success preceded the one-time run.
Seeds 3201--3203 remain reserved and seed 2201 cannot be reused for selection.

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

For the public single-k calibration ablation, control words were sampled outside
planted array intervals, excluded catalogue diagnostic/low-complexity words and
were retained only when their canonical whole-genome occurrence was exactly one.
This use of simulation truth does not represent a blind real-data control panel.
FASTQ survival is the mean product of calibrated Phred-derived correctness
probabilities over valid k-mer windows; FASTA survival uses the supplied
independent-error approximation `(1-error)^k`. Corrected diagnostic depth is
normalized by explicit haploid depth, mean control depth including zero-count
controls, or total read bases/genome size in that precedence order. Median,
MAD and zero fraction are retained for audit. The 10th--90th diagnostic-word
spread is scored only as a diagnostic range.

The depth-gated candidate was fixed at a mean control depth of 2 after the
development run. A separate configuration froze seeds 6401--6403, the same
factorial generator, its input/source hashes and nine gates before generation.
The evaluator ran total-bases and controls modes for every condition, then chose
one complete condition branch without using family truth. The public
`--single-copy-min-depth` implementation was added after validation; it performs
the same observed-depth decision in one command and records the fallback path.
It does not select or certify control k-mers.

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

**Figure 6. Frozen seed-robust classifier development and held-out failure.**
Six panels show the frozen evaluation sequence, development-to-held-out metric
deltas, held-out aggregate metrics, seed-specific deltas, coverage-specific
trade-offs and localization recall by divergence and planted segment count. The
blend improves sensitivity, but increases held-out false-positive rate and
lowers precision; all additional false positives occur at nominal 1×. Inputs,
panel values, hashes and the complete legend are in
`evidence/abundance_classifier_validation_v1/figures_v2` and
`evidence/abundance_classifier_validation_v1/figure_legend.md`. Coverage
stratification is post-hoc diagnosis, and the known-catalogue IID simulation is
not biological collapse truth.

**Figure 7. Fresh held-out validation of the frozen depth-gated classifier.**
Six panels show the failure-to-refinement sequence, development and held-out
metric deltas, held-out confusion counts, seed-specific deltas, coverage-
specific deltas and localization recall by planted divergence and segment count.
The rule passes all predeclared gates on untouched seeds 5801–5803: sensitivity
increases by 0.058985, FPR is unchanged and precision increases by 0.001608.
Inputs, panel values, hashes and the complete legend are in
`evidence/abundance_classifier_depth_gated_validation_v1/figures_v2` and
`evidence/abundance_classifier_depth_gated_validation_v1/figure_legend.md`.
This known-catalogue IID substitution simulation is not biological collapse
truth or a general comparison with external tools.

**Figure 8. Frozen held-out audit of cascade discovery.** Six panels report
mean array recall and raw-call precision across 13 positive scenarios, negative-
read call rates for three controls, paired TandemX/TideHunter wall-time and
direct-child peak-RSS ratios, and all 12 predeclared gates. TandemX passed its
accuracy, determinism, false-call, related-family and RSS requirements, but
promotion failed because nine TRF low-complexity-control executions timed out
and the TandemX/TideHunter wall-time geometric mean was 2.457943 (>2.0).
`evidence/cascade_native_screen_heldout_v1/figures_v2` contains the accepted
editable SVG, PDF/PNG, panel source, complete legend and hashes. Failed TRF rows
are NA rather than fabricated zero-accuracy observations.

**Figure 9. Development audit of public copy-number calibration.** A, aggregate
mean absolute relative error with independent-genome points, including the
post-hoc depth-gated candidate. B, residual signed error. C, error by nominal
coverage, exposing the empirical-control regression at 1x. D, fractions of
paired family conditions that improve, tie or worsen relative to total-bases
normalization. E, original single-execution wall time by coverage; oracle input
is unavailable in blind data. F, truth inclusion within the exported 10th--90th
diagnostic-k-mer spread, with nominal 0.95 shown only to demonstrate lack of
calibration. The accepted editable figure, source rows, full legend and hashes
are in `evidence/quantify_calibration_development_v1/figures_v1`. All controls
were selected with simulation truth, and the candidate is not held-out evidence.

**Figure 10. Frozen validation of depth-gated copy-number normalization.** A,
aggregate MARE with independent-genome points. B, paired genome means. C,
coverage-specific errors for total-bases, controls-only and depth-gated methods.
D, all improved/equal/worse family pairs. E, condition-level use of both frozen
branches. F, runtime and direct-child RSS for all 54 public commands, with group
medians highlighted. The rule passed all nine predeclared gates, but controls
alone retained slightly lower aggregate MARE and all 329 family regressions are
shown. `evidence/quantify_depth_gated_validation_v1/figures_v3` contains the
accepted editable SVG/PDF/PNG, uniquely paired panel source, full legend and
hashes. v1 is a rejected legend-overlap layout; v2 fixed the layout but did not
uniquely key runtime/RSS source rows. Controls used simulation truth, catalogues
were supplied, within-genome conditions are dependent and resource rows are
single executions.

**Figure 11. Frozen validation of the guarded cascade speed path.** A-B,
one-to-one array recall and raw-call precision for TandemX and TideHunter across
13 positive scenarios. C, negative-read call rates for three controls. D-E,
paired median wall-time and direct-child peak-RSS ratios across all 16
scenario-seed groups, with geometric means and predeclared limits. F, all 14
frozen gates. All 96 commands and all gates passed; runtime and RSS ratios were
1.977877 and 0.396864, minimum TandemX base-union F1 was 0.997445 and maximum
boundary MAE was 2.35 bp. Condition-level runtime ratios above 4 are retained.
`evidence/cascade_gap_free_validation_v1/figures_v2` contains the accepted
editable SVG/PDF/PNG, panel source, complete legend and hashes. This is one
synthetic validation seed with technical timing repetitions, not biological
replication or a real-data superiority claim.

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

**Figure S5. TideCluster MorexV3 reference-window scaling and provenance.** A,
internal GNU-time wall time for the TideHunter and clustering stages. B,
container-stage maximum RSS. C, descriptive final interval and operational
family counts. D, union repeat coverage of sampled bases. E, exact, clipped and
merged coordinate-provenance fractions. F, selected representative-period
distributions. Inputs are deterministic nested 1-Mb windows; calls lack
independent accuracy truth and whole-chromosome context. The inspected editable
SVG/PDF/PNG, panel source and hashes are in
`evidence/tidecluster_morex_reference_scaling_v1/figures_v1`.

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
- Supplementary Table S14: classifier development v1, seed-robust development
  v2, complete paired held-out rows, per-seed metrics, localization rows and
  command resource receipts in `evidence/abundance_classifier_validation_v1`.
- Supplementary Table S15: six-genome depth-gated development, fresh 5801–5803
  baseline/localization, raw multi-k pairs, selected classifier rows, per-seed
  metrics and 1,062 command resource receipts in
  `evidence/abundance_classifier_depth_gated_validation_v1`.
- Supplementary Table S16: exact-output Mo17 discovery optimization resources,
  parity status and compact receipts in
  `evidence/discovery_packed_trace_batch_v1`.
- Supplementary Table S17: complete cascade held-out raw/summary rows, paired
  TideHunter ratios, failed runs, gate receipt and Figure 8 panel source in
  `evidence/cascade_native_screen_heldout_v1`.
- Supplementary Table S18: exact four-project ENA queries and the 12 selected
  PacBio genomic-WGS run rows in
  `evidence/retrospective_collapse_source_audit`.
- Supplementary Table S19: all 5,940 public quantify calibration rows, 108
  execution resources, 36 strata, control receipts, per-execution artifact
  hashes, development decision and Figure 9 source data in
  `evidence/quantify_calibration_development_v1`.
- Supplementary Table S20: exact-output hashes for 5,940 metric rows, 108
  copy-number products and three control panels, plus paired method resources and
  replay receipts in `evidence/quantify_calibration_fast_fasta_replay_v1`.
- Supplementary Table S21: all 4,455 frozen-validation metric rows, 1,485 paired
  family outcomes, 54 command resources and artifact hashes, three compact input-
  dataset manifest sets, nine gate observations and Figure 10 panel source in
  `evidence/quantify_depth_gated_validation_v1`.
- Supplementary Table S22: all 96 cascade validation executions, 32
  scenario-tool summaries, 16 paired TandemX/TideHunter rows, all 14 frozen
  gate observations and Figure 11 panel source in
  `evidence/cascade_gap_free_validation_v1`.
- Supplementary Table S23: nested MorexV3 TideCluster stage resources, call
  summaries, interval-provenance classes, copy-number availability and Figure
  S5 source rows in `evidence/tidecluster_morex_reference_scaling_v1`.
- Supplementary Table S24: no-rerun seed-2201 indel-condition comparison,
  including accuracy, boundary, elapsed-time and memory rows with exact source
  hashes in `evidence/indel_detector_gap_audit_v1`.
- Supplementary Table S25: TME204 same-DNA source relationship, raw-read
  accessions and volumes, complete GigaDB inventory result, missing historical
  assembly gate and no-execution fate in
  `evidence/tme204_donor_matched_source_audit_v1`.
- Supplementary Table S26: Heinz 1706 SL5.0/SL-T2T source relationship, exact
  HiFi run and file volumes, official download endpoints, explicit donor-risk
  statement and no-execution fate in
  `evidence/tomato_heinz1706_source_audit_v1`.
- Supplementary Table S27: B73-Ab10 v1/v2 direct same-HMW-DNA statement,
  conflicting archival BioSamples, exact public subread volumes, assembly file
  identity and the target-array reference-ineligibility decision in
  `evidence/b73_ab10_donor_matched_source_audit_v1`.

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
held-out directories listed for Supplementary Table S13; E17:
`abundance_classifier_validation_v1`; E18:
`abundance_classifier_depth_gated_validation_v1`; E19:
`discovery_packed_trace_batch_v1`; E20:
`cascade_native_screen_heldout_v1`; E21:
`retrospective_collapse_source_audit`; E22:
`quantify_calibration_development_v1`; E23:
`quantify_calibration_fast_fasta_replay_v1`; E24:
`quantify_depth_gated_validation_v1`; E25:
`cascade_gap_free_validation_v1`; E26:
`tidecluster_morex_reference_scaling_v1`; E27:
`indel_detector_gap_audit_v1`; E28:
`tme204_donor_matched_source_audit_v1`; E29:
`tomato_heinz1706_source_audit_v1`; E30:
`b73_ab10_donor_matched_source_audit_v1`.
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
7. Pushing the limits of HiFi assemblies reveals centromere diversity between
   two Arabidopsis thaliana genomes.
   [Nucleic Acids Research, doi:10.1093/nar/gkac1115](https://doi.org/10.1093/nar/gkac1115).
8. A complete telomere-to-telomere assembly of the maize genome.
   [Nature Genetics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10335936/).
9. Near-gapless genome assemblies of Williams 82 and Lee cultivars for
   accelerating global soybean research.
   [Plant Genome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807316/).
10. The haplotype-resolved chromosome pairs of a heterozygous diploid African
    cassava cultivar reveal novel pan-genome and allele-specific transcriptome
    features.
    [GigaScience, doi:10.1093/gigascience/giac028](https://doi.org/10.1093/gigascience/giac028).
11. A telomere-to-telomere reference genome assembly of tomato cultivar Heinz
    1706.
    [Plant Communications, doi:10.1016/j.xplc.2025.101618](https://doi.org/10.1016/j.xplc.2025.101618).
12. Conflicting Kinesin-14s in a single chromosomal drive haplotype.
    [Genetics, doi:10.1093/genetics/iyaf091](https://doi.org/10.1093/genetics/iyaf091).

Bibliographic metadata and the full primary dataset references require final
reference-manager curation. No author, funding, conflict or accession-deposition
details are invented in this draft.
