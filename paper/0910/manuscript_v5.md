# TandemX reveals family-specific tandem-repeat under-representation in high-quality plant genome assemblies

**Running title:** Read-based assessment of plant tandem-repeat representation

## Abstract

Long tandem-repeat arrays remain among the least reliably represented components
of plant genomes, even as chromosome continuity approaches telomere-to-telomere
resolution. We asked whether abundant repeat families can remain
under-represented in high-quality assemblies and developed TandemX to compare
family abundance in raw long reads with representation of the same families in
an assembly. In *Arabidopsis thaliana* Ey15-2, all eight repeat families that
gained substantial representation from a historical CLR-Canu assembly to a HiFi
assembly were independently prioritized from HiFi reads; the other 11 eligible
families were not. In *Macadamia jansenii*, the same analysis identified both
families preferentially recovered in the newer assembly and did not label the
other 41. Independent Illumina data supported residual deficits for two Ey15-2
families, and Illumina plus direction-level Nanopore evidence supported one
Macadamia family. Three additional preselected candidates remained unresolved,
revealing how k-mer length and sequence composition can change family-abundance
estimates. Controlled simulations defined the analytical scope of the workflow:
TandemX recovered all 54 planted family conditions in a unified comparison,
whereas SRF recovered 45 at k=151 and 53 at k=101. Ordinary competitive mapping
with the TandemX catalogue produced lower abundance error, and shared-fragment
backgrounds exposed lower catalogue specificity for TandemX. Together, these
results show that raw-read abundance can reveal family-specific repeat deficits
that contiguity statistics do not capture, while explicit cross-k and
cross-platform tests distinguish supported candidates from unstable estimates.

## Introduction

Tandemly repeated DNA is a prominent and rapidly evolving component of plant
genomes. Large arrays occur at centromeres, pericentromeric regions, nucleolar
organizers and other chromosome domains, where their abundance and sequence
organization contribute to genome structure and chromosome variation [1–4].
Their monomer lengths, sequence divergence and copy numbers can differ sharply
among families, chromosomes, accessions and species. Describing a plant genome
therefore requires more than locating a repeat family: it also requires asking
whether the assembly represents its genomic abundance.

Long-read sequencing and telomere-to-telomere assembly have resolved many regions
that were absent from earlier references [3,4,6–8]. Tandem arrays nevertheless
remain difficult because highly similar copies provide few unique anchors and
the arrays can exceed individual sequencing molecules. A chromosome-scale or
gapless assembly can contain a repeat family while still representing fewer
copies than occur in the sequenced material. N50, gap count and genome-wide
k-mer completeness summarize important properties, but none directly answers
this family-level abundance question.

Available methods address complementary parts of the problem. TRF and AniAnn's
annotate tandem arrays in sequence [9,19]; TideHunter detects repeat intervals and
consensus units in long reads [10]; TRASH and TideCluster annotate and cluster
assembly-derived repeats [11,12]; and SRF reconstructs abundant satellite units
from reads or assemblies and estimates mapped-base abundance [13]. Merqury and
KAT compare read and assembly k-mer spectra to assess genome-wide quality and
completeness [16,17], whereas TandemTools uses long reads to assess assembled
extra-long tandem repeats [18]. What remains difficult is to follow a de novo
repeat family from raw-read discovery through abundance estimation and assembly
localization, then report the family-specific discrepancy together with the
evidence needed to interpret it.

We developed TandemX for this purpose. It constructs sequence-supported
operational families from HiFi reads, estimates their abundance with diagnostic
k-mers, localizes the same family catalogue in an assembly and reports the
read--assembly abundance deficit. The linked output distinguishes TandemX from
tools that reconstruct repeat units, annotate arrays or score whole-assembly
completeness as separate tasks (Table 1). Cross-k and orthogonal-platform evidence
is used to separate stable family-specific deficits from estimates that depend on
sequence composition or measurement mode.

We first define the workflow and its analytical performance with controlled,
endpoint-matched comparisons. We then establish its range across ten public HiFi
libraries from eight plant species. The biological tests focus on historical and
newer assemblies of *Arabidopsis thaliana* Ey15-2 and *Macadamia jansenii*, asking
whether read-prioritized families were preferentially recovered by improved
assemblies and whether any deficits persisted in the newer assemblies. Independent
Illumina data, supplemented by Nanopore evidence in Macadamia, support three
residual family-specific deficits. These results establish repeat-family
representation as a dimension of plant assembly quality that remains partly
independent of chromosome continuity.

## Results

### TandemX connects raw-read repeat families to assembly representation

TandemX addresses a distinction that is invisible to contiguity statistics: an
assembly can contain a tandem-repeat family while representing fewer copies than
are supported by the sequenced material (Fig. 1A). The workflow begins with HiFi
reads, detects periodic intervals, derives candidate monomers and groups them by
cyclic, strand-aware sequence similarity (Fig. 1B). Each operational family is
represented by a consensus sequence and its supporting read intervals. Related
representatives and possible period multiples are retained as candidate
relationships, without assigning ancestry or higher-order-repeat structure.

Diagnostic canonical k-mers that distinguish a family from other catalogue
representatives provide a read-derived abundance estimate. The same family
representative is then localized in an assembly, and overlapping intervals are
unioned to measure assembled repeat bases. For read-derived abundance *R* and
assembly representation *A*, TandemX reports *D* = max[*R* - *A*, 0) and the
ratio *A*/*R* (Fig. 1C). Candidate deficits can then be tested across k-mer
lengths or with independent sequencing platforms. A family with discordant
measurements remains unresolved rather than contributing to a genome-wide score.

This linked family-level output is the main distinction from existing tools
(Table 1). Repeat detection, repeat-unit reconstruction, assembly annotation,
long-read mapping and genome-wide k-mer completeness are all established tasks.
TandemX connects them around one operational family catalogue so that a reader
can move directly from a read-supported family to its estimated abundance,
assembly coordinates, representation ratio and evidence state (Fig. 1D).

**Table 1. Native analytical scope of tools relevant to tandem-repeat assembly
assessment.** The table describes standard inputs and outputs rather than every
possible downstream combination.

| Tool or class | Native input | Principal relevant output | Family-resolved read abundance | Linked family-level read--assembly audit |
| --- | --- | --- | --- | --- |
| TideHunter | Long reads | Per-read tandem intervals and consensus units | No catalogue-level estimate | No |
| SRF | Reads or assemblies | Reconstructed satellite units/HORs and mapped-base abundance | Yes, for reconstructed units | No linked assembly audit |
| TRASH/TRASH2, TideCluster, AniAnn's | Assembly sequence | Tandem-array annotation, clustering or repeat structure | No read-derived estimate | No |
| TandemTools | Assembly and long reads | Mapping, polishing and quality assessment of assembled extra-long tandem repeats | Not a primary output | Assembly assessment without a de novo abundance catalogue |
| Merqury/KAT | High-accuracy reads and assembly | K-mer spectra and whole-assembly quality summaries | No repeat-family catalogue | Genome-wide rather than family-resolved |
| TandemX | HiFi reads; assembly optional | Operational family catalogue, abundance, localization, deficit and evidence state | Yes | Yes |

### Controlled benchmarks establish analytical performance and method complementarity

We evaluated repeat discovery, abundance estimation and assembly localization
against planted truth before combining them in a unified comparison (Fig. 2A).
In a factorial 10-Mb simulation containing 55 planted families, TandemX recovered
all founders from clean and high-error 5x long-read datasets. A broader validation
spanning period, divergence, indels, array number and negative controls retained
high array and family recovery and made no calls in the three purely negative
classes. Sequence-based clustering also recovered all three related planted
families in a dedicated validation.

We next measured how discovery changed from 0.5x to 30x in nested read sets from
three independently generated 10-Mb validation genomes (Supplementary Fig. S1;
Supplementary Tables S3 and S4). Median operational-family count increased from
34 at 0.5x to 55 at 5x and 58 at 30x, whereas the median number of unique
candidate monomers continued to increase from 36 to 119. All 55 planted
families, including the 20-copy low-abundance tier, were recovered in every
genome at 5x, and no additional planted family was recovered at higher depth.
Under the pre-specified catalogue criterion, 5x--10x and 10x--20x both had a
new-family fraction below 0.05 and a family-level Jaccard above 0.95 in all
three genomes; 20x was therefore the first operational saturation depth. The
20x--30x transition also passed in all genomes. This result separates early
recovery of the planted families from the more conservative stabilization of
the complete operational catalogue.

Quantification was scored with the planted catalogue supplied so that abundance
error remained separate from discovery error (Fig. 2B). Across three independent
10-Mb validation genomes and 1,485 family conditions, depth-gated normalization
reduced mean absolute relative error from 0.409 to 0.363. The largest benefit
occurred in the deliberately difficult 20x, 2%-divergence, high-error condition,
where multi-k estimation reduced mean error from 56.6% at k=21 to 13.5%.
Diagnostic support was insufficient in some low-coverage conditions, so the
workflow reports availability and dispersion with each estimate.

Assembly localization remained accurate when nearby anchors were bridged across
no more than one monomer length (Fig. 2C). On new simulated genomes, mean recall
was 0.981 and mean precision was 0.9996; recall in the weakest 5%-divergent,
three-segment stratum was 0.949. The downstream binary comparison was more
sensitive to abundance error: in the held-out test, sensitivity was 0.654, false-
positive rate was 0.016 and precision was 0.983.

Task-matched comparisons showed substantial overlap with existing tools rather
than one universal winner (Fig. 2D). TandemX, TRF and TideHunter recovered all 55
founders in the read simulation, although one-to-one high-error array precision
was 1.000, 0.699 and 0.920, respectively. Base-union precision exceeded 0.996 for
all three, showing that overlapping and harmonic calls affected interval counts
more than covered bases. TideHunter was generally faster. On assemblies, TRASH2
recovered all planted families, and TideCluster achieved high base-union accuracy
in completed runs.

A unified experiment compared TandemX, SRF and ordinary competitive mapping on
18 newly generated read datasets spanning six conditions (Fig. 2E;
Supplementary Table S2). TandemX recovered all 54 planted family conditions,
whereas SRF recovered 45 at k=151 and 53 at k=101. Positive-family abundance
accuracy followed a different ranking: competitive mapping with the TandemX
catalogue had the lowest error in every condition, and SRF was more accurate than
TandemX under substitutions and 2% unit divergence. In shared-fragment
backgrounds, TandemX annotated a mean 30,887 negative-read bases per dataset,
compared with 0 for SRF k=151, 227 for SRF k=101 and 11,166 for competitive
mapping. The combined benchmarks therefore define TandemX as a connected
family-level audit with strong recovery, while abundance estimation and catalogue
specificity remain endpoint-dependent.

### TandemX operates across diverse plant repeat landscapes

We next asked whether the workflow could be applied across the range of genome
sizes and repeat contents encountered in plant genomics. Ten complete HiFi
libraries from eight reported species passed checksum, FASTQ traversal and
record/base reconciliation, together comprising 14,937,608 reads and 262.7 Gb
of sequence (Fig. 3A; Supplementary Table S1). The cohort included Arabidopsis,
rice, maize, barley, rye, wheat, oat and wild soybean, with individual libraries
ranging from 5.62 to 77.09 Gb and median read lengths from 13.46 to 21.51 kb.

TandemX produced repeat calls in each tested plant dataset, including wheat,
barley, rye and oat libraries from large, repeat-rich genomes (Fig. 3B-D).
Catalogue size, monomer length, supporting-read fraction and repeat-base fraction
varied among libraries and input scales. Nested subsamples retained the major
input-specific patterns while increasing the number of lower-support calls at
larger scale. These profiles demonstrate that the same workflow can generate
family catalogues across contrasting plant datasets without loading complete
libraries into memory.

The cross-species profiles are used here to establish operating range and to
display the diversity of candidate repeat families. Several libraries represent
pooled plants or technical batches, and no matched family-level truth exists
across the cohort. The biological assembly tests therefore use the two systems
with the strongest available historical/newer assembly relationship rather than
treating the ten libraries as replicate estimates of species differences.

### Improved assemblies preferentially recover repeat families prioritized from raw reads

The first biological test asked whether raw-read abundance could identify repeat
families later recovered by improved assemblies (Fig. 4A). In *A. thaliana*
Ey15-2, we compared the paper-supplied Bionano-scaffolded CLR-Canu assembly with
the newer HiFi-Hifiasm assembly from the same accession and sample designation
[4]. Nineteen read-derived families had at least 15 kb in the newer assembly and
entered the primary analysis. Eight showed an old/new assembly ratio below 0.6.
The HiFi-read comparison independently prioritized all eight against the
historical assembly and did not label the other 11 (Fig. 4B).

Representation changed most clearly at the family level. Across the eligible
Ey15-2 families, the total HiFi-versus-old abundance deficit was 840,658 bp and
the new-versus-old assembly gain was 973,959 bp. The family-level Pearson
correlation was 0.589, and the mean absolute difference was 45,551 bp, indicating
agreement in which families gained rather than exact agreement in gain magnitude
(Fig. 4C). At a permissive 5-kb eligibility threshold the larger family set
contained one missed and two additional calls; six families remained at 50 kb.

We repeated the analysis in the phylogenetically distant *M. jansenii*, using an
earlier CLR Falcon-Unzip primary assembly and a newer HiFi IPA primary assembly
(Fig. 4D). The newer assembly was more contiguous, with 284 rather than 762
primary contigs, although it was slightly shorter overall. Of 1,227 read-derived
families, 43 had at least 15 kb represented in the newer assembly. Two families
showed strong historical under-representation; TandemX prioritized both and did
not label the other 41. The result was unchanged when read depth was normalized
by the paper-reported 28-Gb yield instead of the public FASTQ total.

The Macadamia HiFi-versus-old abundance deficit was 2.12 Mb, whereas the
new-versus-old gain among eligible families was 0.159 Mb. Family-level
correlation was 0.535 and the mean absolute difference was 50.3 kb (Fig. 4E).
Thus both plant systems supported the same qualitative result: families with a
large raw-read signal relative to a historical assembly were the families that
preferentially gained representation in an improved assembly. The discordant
magnitudes led us to test remaining deficits with independent reads.

### Independent reads support residual deficits in three high-quality assembly families

We selected six newer-assembly candidates before examining orthogonal data:
TXF000002, TXF000154 and TXF001517 in Ey15-2, and TXF000496, TXF000563 and
TXF000695 in Macadamia (Fig. 5A; Table 2). Family definitions, the 15-kb
eligibility criterion and the 0.6 assembly/read ratio were held constant.
Illumina abundance was evaluated at k=21 and k=31; Macadamia also had
direction-level PromethION evidence.

Ey15-2 TXF000002 had an Illumina abundance estimate of 519,994-555,544 bp,
compared with 49,758 bp localized in the newer assembly. TXF000154 had an
Illumina estimate of 33,895-44,496 bp, compared with 18,297 bp in the assembly.
Both remained below the assembly/read threshold at both k values (Fig. 5B).
None of the other 16 eligible Ey15-2 context families crossed the orthogonal
deficit rule at either k.

Macadamia TXF000496 provided the strongest cross-platform example (Fig. 5C).
Its Illumina estimate was 1.586-1.604 Mb, whereas 55,562 bp was localized in the
newer assembly. Calibrated ONT occupancy corresponded to 0.620-0.785 Mb across
the genome-size sensitivity and supported the same direction. None of the other
40 eligible Macadamia context families crossed the Illumina rule at both k
values. The two Ey15-2 families and TXF000496 therefore form three independently
supported cases in which a newer HiFi assembly retains substantially less family
sequence than is indicated by a second sequencing library or platform.

The remaining candidates show why sequence-level diagnostics must accompany the
headline result (Fig. 5D-F). Ey15-2 TXF001517 changed interpretation between
k=21 and k=31. Macadamia TXF000563 crossed the 0.6 boundary at k=21 but not k=31,
and ONT did not support a deficit. TXF000695 showed the strongest instability:
its 467 diagnostic words produced median Illumina depths of 77,091 at k=21 and
101 at k=31, yielding abundance estimates of 1,104,859 and 1,649 bp. HiFi and
ONT mapping instead supported the low-abundance direction.

The TXF000695 distribution was bimodal rather than driven by one extreme word.
At k=21, 259 of 467 words occupied the high-depth component; at k=31, only 202
did so and the median shifted to the low-depth component. Shorter words may be
sampling a conserved segment shared with other genomic sequence, although the
current data cannot distinguish that explanation from mosaic family structure
or other sources of homology. TXF001517, TXF000563 and TXF000695 therefore
remain unresolved. Together with the three supported families, they show that
independent platforms can identify reproducible residual deficits while cross-k
disagreement exposes repeat families that require sequence or cytological
follow-up.

**Table 2. Family-level evidence from historical/newer assemblies and orthogonal
reads.** Counts use the primary 15-kb newer-assembly eligibility criterion.

| System | Eligible families | Historical families preferentially recovered in newer assembly | Correctly prioritized from raw reads | Orthogonal candidates | Supported residual deficits | Unresolved candidates |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *A. thaliana* Ey15-2 | 19 | 8 | 8 | 3 | 2 | 1 |
| *M. jansenii* | 43 | 2 | 2 | 3 | 1 | 2 |

## Discussion

The central biological result is that repeat-family representation can remain
incomplete after large gains in assembly continuity. In Ey15-2 and Macadamia,
raw-read abundance identified exactly the primary families that were
preferentially recovered when CLR-era assemblies were replaced by HiFi
assemblies. Independent reads then supported residual deficits for two Ey15-2
families and one Macadamia family. Assembly continuity and family-level repeat
representation are therefore related but distinct dimensions of genome quality.

The historical/newer comparisons provide a useful retrospective test because
they ask whether a prediction made from reads points to sequence recovered by a
subsequent assembly. The agreement was exact at the primary family threshold in
both species: eight of eight Ey15-2 families and two of two Macadamia families
were prioritized, with no calls among the remaining 52 eligible families. The
continuous estimates were less concordant with assembly gains. TandemX is thus
most strongly supported as a method for identifying which families merit review,
rather than for reconstructing the precise number of bases absent from an older
assembly.

Orthogonal data extend the result beyond a comparison between assemblies. The
Ey15-2 Illumina library changed extraction, library construction and sequencing
technology while retaining the source tissue pool. Macadamia Illumina reads
provided a separate diagnostic-k-mer measurement, and PromethION reads changed
the measurement mode to long-read occupancy. Concordant evidence at two k-mer
lengths supported TXF000002 and TXF000154; concordant Illumina and ONT direction
supported TXF000496. These three cases show that substantial family-specific
deficits can persist in newer HiFi assemblies.

TandemX makes this question accessible by preserving the identity of a repeat
family across discovery, abundance estimation and assembly localization. SRF
reconstructs repeat units and estimates their abundance, TRASH and related tools
annotate assembled arrays, TandemTools examines read support for assembled
extra-long repeats, and Merqury/KAT assess genome-wide k-mer completeness. The
distinct output of TandemX is a linked family record that combines the raw-read
signal, assembly representation and evidence state. This organization turns a
general observation that repeats challenge assembly into a list of specific
families that can be inspected, reassembled or tested experimentally.

The controlled comparisons also define where this output is most useful.
TandemX showed strong family recovery and accurate assembly localization, but it
was not uniformly best at every endpoint. TideHunter was faster in the tested
read comparisons; TRASH2 and TideCluster performed strongly on assembly
annotation; SRF produced more specific catalogues in the shared-fragment
background; and competitive mapping gave the lowest abundance error when the
TandemX catalogue was supplied. These results support a task-specific view of
repeat analysis. TandemX contributes the connected audit and family
prioritization, while established tools remain preferable for several component
tasks.

The unresolved candidates reveal a second biological and analytical issue:
different parts of a repeat consensus may have very different genomic
specificity. TXF000695 contained two diagnostic-word depth components, and the
median changed components when k increased from 21 to 31. This pattern is
consistent with a repeat containing locally conserved sequence shared with a
larger genomic background, although mosaic family structure and other forms of
homology remain possible. Cross-k stability, observed-word fraction and depth
dispersion provide practical indicators for distinguishing robust family signals
from such composition-sensitive estimates.

The present evidence has clear boundaries. The old/new analysis covers two
species, the orthogonal analysis was restricted to six preselected candidates,
and no physical copy-number measurement is available. Macadamia public metadata
support a close sample relationship but not identical DNA extraction for every
library. Simulations simplify ploidy, population heterogeneity and long-range
array organization, while ONT occupancy supplies direction rather than exact
copy number. The supported cases should therefore be interpreted as
family-specific under-representation rather than a prevalence estimate for all
plant assemblies.

These results nominate a direct next biological test. Candidate families can be
mapped to chromosome domains, examined with targeted long-read or assembly
analysis and validated cytologically with FISH. A repeat discovered from reads,
poorly represented in an assembly and independently localized on chromosomes
would connect the abundance deficit to a visible genomic structure. The current
study establishes the analytical and cross-platform basis for that experiment;
the discovery and validation of additional repeat families remain future work.

## Methods

### Study design and terminology

The study was organized around three levels of evidence. Controlled simulations
tested discovery, abundance estimation, localization and binary comparison
against planted truth. Historical/newer assembly pairs tested whether
read-predicted families preferentially gained representation as assembly data
and methods improved. Orthogonal Illumina and ONT data then tested candidate
deficits in the newer assemblies without assuming that either assembly was
absolute copy-number truth.

An operational repeat family is a group of candidate monomers connected to a
fixed representative by the documented cyclic, strand-aware sequence criterion.
It is an analytical unit, not an evolutionary clade. `Read-derived abundance`
denotes the repeat bases implied by diagnostic-k-mer depth after the stated
normalization. `Assembly representation` denotes unioned bases localized for the
same family. `Read--assembly abundance deficit` is the positive difference
between these quantities. `Estimated under-representation` is the same
continuous quantity interpreted within the evidence limits. A binary candidate
state was defined by an assembly/read ratio below 0.6 in the pre-specified experiments.
`Residual under-representation` required concordant orthogonal support;
discordant k values, platforms or sensitivities were `unresolved`.

### TandemX discovery and family construction

TandemX scans long reads for local periodicity over the stated period range,
refines candidate array intervals by alignment and derives a consensus monomer
for each supported interval. Candidate monomers are canonicalized across
rotation and reverse complement before sequence-based clustering. Families use
fixed representatives rather than transitive merging, and supporting read
intervals, member assignments and ambiguous relationships are retained.
Candidate relationships among representatives are reported separately: a
shorter-to-longer edge can be labelled as a possible period multiple when its
length ratio is close to an integer, but no unit order, biological root or
higher-order-repeat status is inferred.

Discovery validation used independently generated genomes spanning monomer
periods, array copy numbers, GC composition, unit divergence, substitutions,
indels, multiple arrays per read and negative controls. Array recovery required
one-to-one interval matching with intersection-over-union of at least 0.5 and
period error no greater than 2 bp or 2% of the true period, whichever was
larger. Sequence-family recovery used cyclic edit identity with reverse-
complement equivalence. Exact simulation configurations, split identities and
all adverse development results are provided in Supplementary Methods and
Supplementary Results.

### Coverage-saturation validation

Three independently generated 10-Mb genomes used validation seeds 6501, 6502
and 6503. Each contained 54 factorial families spanning periods of 61, 171 and
421 bp, copy counts of 20, 80 and 200, GC fractions of 0.3, 0.5 and 0.7, and
unit divergence of 0 or 0.02, together with one 171-bp, 6,000-copy family.
Reads followed the complete Mo17 HiFi length distribution and one fixed IID
error model with substitution, insertion and deletion probabilities of 0.001
each. For each genome, the 0.5x, 1x, 2x, 5x, 10x, 20x and 30x datasets were
nested prefixes of one simulated read stream. Depths were paired within a seed;
the three genome seeds were the independent simulation units.

Discovery used periods of 30--1,000 bp, a minimum repeat span of 100 bp,
minimum support of one read, elastic interval refinement, sequence clustering
at 0.95 identity and the related-family audit. Candidates and operational
families were counted after cyclic and reverse-complement canonicalization.
Planted-family recall and adjacent-depth catalogue matching used independent
one-to-one cyclic global edit matching at identity at least 0.90. The abundance
tiers were fixed by planted copy count: low, 20 copies; medium, 80 copies; and
high, at least 200 copies. For adjacent depths with *M* matched families and *L*
and *H* lower- and higher-depth families, catalogue Jaccard was *M*/(*L* + *H*
- *M*) and new-family fraction was (*H* - *M*)/*H*. A transition passed when
new-family fraction was below 0.05 and Jaccard was above 0.95. Operational
saturation required two consecutive transitions to pass in all three seeds;
the earliest upper depth was reported. This rule, all seeds, all depths and the
matching thresholds were committed before validation outputs were examined.

### Diagnostic-k-mer abundance estimation

For each family and k, TandemX constructs canonical circular k-mers from the
representative, removes low-complexity words and retains words not assigned to
another catalogue family. Counts are divided by their multiplicity in the
representative. The primary single-k estimator summarizes diagnostic-word depth
and converts it to family copy number using either explicit haploid depth,
empirical single-copy controls or total read bases divided by haploid genome
size, in that order of precedence. Estimated copy number is multiplied by
monomer length to obtain read-derived family abundance. Median depth, median
absolute deviation, zero or unobserved fractions and warnings are retained.

The conditional multi-k analysis used k=15, 21, 27 and 31. A log-linear trend
across supported k values estimated the intercept at zero word length; a
single-k value was retained when multi-k support was insufficient. Joint
read-level first and cross moments propagated the shared-read sampling variance
through this fit. These intervals condition on the supplied genome size and
catalogue and do not include family error, normalization error or systematic
extrapolation bias. The final depth gate retained the conservative single-k
rule below estimated haploid depth 2 and used the fixed multi-k blend at higher
depth. No model was refitted in the old/new or orthogonal plant analyses.

### Assembly localization and comparison

Catalogue representatives were localized to assemblies by exact diagnostic
anchors followed by bounded extension and bridging. Overlapping intervals were
unioned within family before their bases were summed. The final simulated
localizer allowed anchor gaps of no more than one monomer length and applied the
fixed identity criterion. Localization truth was scored as base-union recall
and precision, separately from family recovery and binary comparison.

For a family with read-derived abundance *R* and assembly representation *A*,
TandemX reports max[*R* − *A*, 0] and the ratio *A*/*R* when *R* is positive. The fixed
binary candidate threshold was *A*/*R* < 0.6. This threshold is an analytical rule
evaluated in the reported benchmark, not a universal biological boundary.
Insufficient read evidence, localization warnings and denominator exclusions
remain explicit family fates.

### Comparator evaluation

Read-level comparisons used identical FASTA inputs and matched period and
minimum-array scopes for TandemX, TRF and TideHunter. Raw interval precision,
one-to-one array recovery, base-union recovery and cyclic family recovery were
scored separately. Assembly-side evaluations used TRASH/TRASH2 and TideCluster
with their native output relationships retained. Approximate region summaries,
unit-level intervals, primary and secondary consensuses, and merged intervals
were not treated as interchangeable. A failed, timed-out, malformed or
unstarted process yielded unavailable evidence rather than zero accuracy.
Complete settings, resource measurements and failure histories are in
Supplementary Methods.

The unified SRF comparison was specified before generating seeds
2026091001–2026091003. Each of six conditions contained 100 approximately 5-kb
reads, three unrelated families and 12 units per planted array, with 70% positive
reads except the 10% low-abundance condition. Conditions were clean 171-bp units,
1% substitutions, 0.1% insertions plus 0.1% deletions, 2% unit substitutions,
low abundance, and a 120-bp-unit shared-fragment background. Negative reads in
the last condition contained a 100-bp fragment from one family separated by
independent 30–80-bp random gaps. Observed post-mutation coordinates defined
positive truth; neutral input read IDs contained no labels.

TandemX used the existing cascade detector, sequence clustering, periods
30–1000 bp, minimum span 100 bp and support of one read. Its Rust diagnostic
counter used k=21 and depth=1; no simulated error rate was supplied to the
FASTA quantifier. SRF used its native KMC–SRF–mapping–filtering–abundance
workflow, with count cutoff 20 and both k=151 and k=101 reported. A clean KMC
rebuild passed independent canonical-count checks at k=17,101,151. Competitive
mapping used 10-kb tandemized native TandemX representatives and minimap2
map-hifi, retaining primary and secondary alignments with block length >=100 bp
and identity >=90%. Cross-native-family overlaps were withheld from abundance
before truth correspondence but retained in global interval precision.

Independent global edit correspondence allowed circular rotations, reverse
complements and pure integer-repeat units up to 64 times at >=90% identity.
Multiple matching families remained ambiguous. Native SRF retained mapped bp
was verified against positive-keep BED intervals. Abundance errors included all
positive-truth families; unmatched native amounts and negative-read interval
attribution were reported separately. All 72 method/configuration cells
completed, including three valid empty SRF catalogues. Native stages ran serially;
wall time, CPU and maximum child RSS were recorded per dataset. These three
simulation seeds were not technical timing repeats or a large-input scaling study.

### Plant long-read cohort

Complete public long-read files were accepted only after source byte count and
checksum verification, gzip and FASTQ traversal through end of file, record and
base reconciliation, and an exact duplicate-identifier audit. Ten HiFi
libraries from eight reported species were retained. Deterministic whole-file
sampling produced nested input sizes for descriptive catalogue and comparator
analyses. Library, BioSample and source-paper provenance were kept separate;
technical batches and pooled plants were not counted as independent biological
replicates. Reference and organellar content were recorded where available,
because total-library bases do not necessarily equal nuclear sequencing depth.

### Ey15-2 and Macadamia historical/newer assembly analyses

The Ey15-2 analysis used complete public HiFi reads and the paper-supplied
Bionano-scaffolded CLR-Canu and HiFi-Hifiasm assemblies. One HiFi-derived
catalogue was quantified and independently localized to each assembly. The
primary denominator required at least 15,000 localized bases in the newer
assembly; 5,000- and 50,000-bp denominators were retained as sensitivities. The
old/new assembly reference-proxy label used the same 0.6 ratio as the read/old
comparison. A combined final assembly and explicit 107× normalization were
sensitivity analyses only.

The Macadamia analysis used two complete HiFi runs, the purged CLR Falcon-Unzip
primary assembly and the HiFi IPA primary assembly. The update study's
paper-level same-sample statement supported enrollment, while differing archival
BioSamples precluded a same-extraction claim. Discovery used a deterministic
whole-library sample and abundance used both complete read files. The primary
normalization used public read bases and the reported 780-Mb genome size; a
separate sensitivity used the paper-reported 28-Gb yield. Eligibility,
localization and binary rules matched Ey15-2.

For explanatory context only, newer assemblies were aligned to their historical
counterparts and family intervals were summarized by primary aligned-query
coverage. Ey15-2 allowed same-chromosome scopes; Macadamia did not because old
and new contig identifiers were not shared. These alignments did not define
family eligibility or truth.

### Orthogonal Illumina and ONT validation

The orthogonal analysis retained the unchanged HiFi catalogues, newer-assembly
localization and 15-kb denominators. Ey15-2 used PCR-free Illumina run
ERR8666067. Macadamia used Illumina run SRR11191912 and PromethION run
SRR11191910. Input objects were accepted only after official checksum,
`vdb-validate` and converted record/base reconciliation. The Ey15-2 library
contained 157,774,340 reads and 23,666,151,000 bases; Macadamia Illumina and ONT
contained 225,016,144 reads/33,669,394,677 bases and 2,841,932
reads/23,186,565,438 bases, respectively.

Illumina abundance was estimated independently at k=21 and k=31 from
family-exclusive diagnostic k-mers. The denominator was the mean depth of
20,000 deterministic controls present exactly once in both enrolled assemblies
and absent from all catalogue monomers. Counts below the minimum reported by the
counter were conservatively treated as the midpoint of zero and one. Control
median, dispersion and zero fraction were retained as diagnostics. Original
base qualities were not used in the sequence-only count.

For Macadamia ONT validation, the 43 eligible representatives were tandemized
and reads were mapped competitively. Primary family assignments required at
least 500 aligned bases and 0.75 identity; reads assigned to multiple families
were excluded. Occupancy was corrected by a single mapping-efficiency factor
fixed from the corresponding HiFi mapping before candidate ONT results were
examined. ONT direction was evaluated at haploid genome sizes of 616, 653, 738
and 780 Mb. The mapping estimate was used only to determine whether read
occupancy supported more or less family sequence than the assembly.

Residual under-representation required an assembly/orthogonal-read ratio below
0.6, agreement with the HiFi direction and stability across k and applicable
genome-size sensitivities. The strongest Macadamia interpretation additionally
required concordant ONT direction. A quantification-bias interpretation
required the orthogonal estimate to be compatible with the newer assembly while
the fixed HiFi analysis indicated a deficit and exceeded orthogonal abundance by
at least 1.5-fold. All other outcomes were unresolved. Binary interpretation was
reported separately from deficit magnitude.

### Cross-k diagnostic assessment of TXF000695

The complete diagnostic-k-mer count distributions were retained for the
preselected Macadamia families. For TXF000695, the 467 family k-mers at each k
were joined to the exact sequence counts; words absent from the counter output
were assigned the same fixed 0/1 midpoint used by the abundance analysis. We
reported diagnostic-set size, fraction observed at least twice, median, median
absolute deviation, MAD/median and the cross-k fold differences in median depth
and estimated abundance. Counts were also summarized by fixed depth ranges to
determine whether the cross-k reversal reflected one outlier or a redistribution
of the diagnostic set. The labels `unstable_across_k`,
`depth_median_boundary_sensitive` and `reliability=low` are descriptive QC
annotations for this manuscript. They were not used to change family selection,
the estimator, the 0.6 threshold or the pre-specified interpretation rules.

### Statistical analysis and reproducibility

The family was the unit of the old/new and orthogonal comparisons. Binary
sensitivity, false-positive rate and precision were kept separate from
continuous abundance error. Correlations between read-versus-old deficits and
new-versus-old assembly gains were descriptive and did not treat the latter as
truth. Simulation summaries distinguished independent source genomes from
dependent coverage/error conditions and technical timing repetitions.

Analyses were specified before the relevant validation data were examined where
stated. Independent implementations recomputed the Ey15-2 and Macadamia family
tables and the successful TideCluster accuracy endpoints. Commands, software
versions, validation splits, input identities and resource profiles accompany
the source data.

## Data and code availability

The TandemX source code and committed source data are available at
[github.com/CongyangY/TandemX](https://github.com/CongyangY/TandemX). Public read
and assembly accessions are given in Methods and the source-data manifests.
Large FASTQ/FASTA files and intermediate mapping or k-mer databases can be
regenerated from the recorded public sources and are not versioned in Git.

## Figure legends

**Figure 1. TandemX links raw-read tandem-repeat families to assembly
representation.** (A) Conceptual distinction between chromosome continuity and
family-level repeat representation. (B) HiFi reads are converted into periodic
intervals, candidate monomers and a sequence-supported operational family
catalogue. (C) Diagnostic-k-mer depth estimates read-derived family abundance
(*R*), while localization of the same representative measures assembly
representation (*A*) and the read--assembly abundance deficit. (D) The family
report connects abundance, assembly coordinates, the *A*/*R* ratio and
cross-k or orthogonal evidence. Candidate relationships among representatives
are shown separately from the operational family assignments.

**Figure 2. Controlled validation defines the analytical scope of TandemX.**
(A) Simulation design and the separate truth endpoints for discovery,
quantification, localization and binary comparison. (B) Family and interval
recovery by TandemX, TRF and TideHunter across clean and high-error reads. (C)
Abundance error across coverage, read-error and sequence-divergence conditions.
(D) Base-union recall and precision for assembly localization, including
divergent and interrupted arrays. (E) Unified TandemX, SRF and competitive-
mapping comparison, showing family recovery, positive-family MARE and
negative-read attribution as distinct endpoints. Runtime and peak-memory
measurements are reported in Supplementary Figure S5.

**Figure 3. TandemX operates across diverse plant long-read datasets.** (A)
Species, accessions and total HiFi bases for ten complete libraries from eight
plant species. (B) Median and N50 read lengths. (C) Tandem-repeat interval call
density and positive-read fraction in matched-size subsamples from wheat,
barley, rye and oat. (D) Repeat-base fraction across nested input sizes,
including the 1.17-Gb barley workload. These comparisons describe the observed
data and operating range; they are not cross-species accuracy estimates.

**Figure 4. Repeat families prioritized from raw reads are preferentially
recovered by improved assemblies.** (A) Historical/newer assembly comparison
for Ey15-2 and Macadamia. (B) Ey15-2 family representation in the historical and
newer assemblies, with the eight primary proxy-positive families highlighted.
(C) Ey15-2 read--historical deficits versus newer--historical gains and
eligibility-threshold sensitivities. (D) Corresponding Macadamia family
representation, highlighting the two primary proxy-positive families. (E)
Macadamia deficit-versus-gain comparison and genome-size normalization
sensitivity. (F) Summary of primary family-level prioritization in the two
species.

**Figure 5. Independent reads support three residual family-specific deficits
and resolve the confidence boundary of three others.** (A) Evidence design for
the six candidates selected from the newer assemblies. (B) HiFi, Illumina
k=21/k=31 and assembly estimates for Ey15-2 TXF000002 and TXF000154. (C)
Illumina and direction-level ONT evidence for Macadamia TXF000496. (D) Cross-k
and ONT results for the unresolved TXF001517 and TXF000563 candidates. (E)
Complete k-mer depth distributions for TXF000695 at k=21 and k=31, showing the
median transition between low- and high-depth components. (F) Family-level
classification of the six candidates as orthogonally supported or unresolved.

## Supplementary figure legends

**Supplementary Figure S1. TandemX discovery saturation across three simulated
validation genomes.** (A) Unique candidate monomers and operational families
from 0.5x to 30x. (B) New operational families per added 1x and the fraction
new relative to the higher-depth catalogue. (C) One-to-one family-level Jaccard
between adjacent nested depths. (D) Planted-family recall for low (20-copy),
medium (80-copy) and high (at least 200-copy) families. Lines show medians and
ribbons show ranges across three independently generated genomes. Dashed lines
mark the pre-specified new-family fraction of 0.05 and Jaccard of 0.95;
the vertical dotted line marks the first depth satisfying the two-consecutive-
transition, all-seed saturation rule.

## Supplementary table legends

**Supplementary Table S1. Public plant HiFi libraries used to establish the
operating range of TandemX.** The table reports species and material labels, run
accessions, complete-file read and base counts, read-length summaries, GC and N
fractions and input-validation results for ten libraries from eight reported
plant species. Source: `paper/tables/input_cohort.tsv`.

**Supplementary Table S2. Unified abundance comparison on newly generated
planted-read data.** Values are mean MARE across three independent simulated
datasets per condition. All positive-truth families, including undetected
families, enter each mean; unmatched native abundance and negative-read
attribution are reported separately.

**Supplementary Table S3. Discovery yield and abundance-stratified recovery
across nested read depths.** Values are medians and ranges across three
independently generated validation genomes. Candidate and family counts were
canonicalized across rotation and strand; recall used one-to-one cyclic
sequence matching to the 55 planted families.

**Supplementary Table S4. Marginal family discovery and catalogue stability
between adjacent read depths.** Values are medians and ranges across three
validation genomes. The table reports both new operational families and newly
recovered planted families. Passing-seed counts apply the pre-specified
new-family-fraction and Jaccard criteria.

## References

1. Garrido-Ramos MA. Satellite DNA: an evolving topic. *Genes* 8, 230 (2017).
   [doi:10.3390/genes8090230](https://doi.org/10.3390/genes8090230).
2. Melters DP et al. Comparative analysis of tandem repeats from hundreds of
   species reveals unique insights into centromere evolution. *Genome Biology*
   14, R10 (2013).
   [doi:10.1186/gb-2013-14-1-r10](https://doi.org/10.1186/gb-2013-14-1-r10).
3. Naish M et al. The genetic and epigenetic landscape of the Arabidopsis
   centromeres. *Science* 374, eabi7489 (2021).
   [doi:10.1126/science.abi7489](https://doi.org/10.1126/science.abi7489).
4. Rabanal FA et al. Pushing the limits of HiFi assemblies reveals centromere
   diversity between two *Arabidopsis thaliana* genomes. *Nucleic Acids
   Research* 50, 12309–12327 (2022).
   [doi:10.1093/nar/gkac1115](https://doi.org/10.1093/nar/gkac1115).
5. The Arabidopsis Genome Initiative. Analysis of the genome sequence of the
   flowering plant *Arabidopsis thaliana*. *Nature* 408, 796–815 (2000).
   [doi:10.1038/35048692](https://doi.org/10.1038/35048692).
6. Wang B et al. High-quality *Arabidopsis thaliana* genome assembly with
   Nanopore and HiFi long reads. *Genomics, Proteomics & Bioinformatics* 20,
   4–13 (2022).
   [doi:10.1016/j.gpb.2021.08.003](https://doi.org/10.1016/j.gpb.2021.08.003).
7. Navrátilová P et al. Prospects of telomere-to-telomere assembly in barley:
   analysis of sequence gaps in the MorexV3 reference genome. *Plant
   Biotechnology Journal* 20, 1373–1386 (2022).
   [doi:10.1111/pbi.13816](https://doi.org/10.1111/pbi.13816).
8. Shang L et al. A complete assembly of the rice Nipponbare reference genome.
   *Molecular Plant* 16, 1232–1236 (2023).
   [doi:10.1016/j.molp.2023.08.003](https://doi.org/10.1016/j.molp.2023.08.003).
9. Benson G. Tandem repeats finder: a program to analyze DNA sequences.
   *Nucleic Acids Research* 27, 573–580 (1999).
   [doi:10.1093/nar/27.2.573](https://doi.org/10.1093/nar/27.2.573).
10. Gao Y, Liu B, Wang Y, Xing Y. TideHunter: efficient and sensitive tandem
    repeat detection from noisy long reads using seed-and-chain. *Bioinformatics*
    35, i200–i207 (2019).
    [doi:10.1093/bioinformatics/btz376](https://doi.org/10.1093/bioinformatics/btz376).
11. Wlodzimierz P, Hong M, Henderson IR. TRASH: Tandem Repeat Annotation and
    Structural Hierarchy. *Bioinformatics* 39, btad308 (2023).
    [doi:10.1093/bioinformatics/btad308](https://doi.org/10.1093/bioinformatics/btad308).
12. Novák P et al. TideCluster: tandem-repeat annotation and clustering
    software. Source and versioned software record available from
    [GitHub](https://github.com/kavonrtep/TideCluster) (accessed 10 September 2026).
13. Zhang Y, Chu J, Cheng H, Li H. De novo reconstruction of satellite repeat
    units from sequence data. *Genome Research* 33, 1994–2001 (2023).
    [doi:10.1101/gr.278005.123](https://doi.org/10.1101/gr.278005.123).
14. Murigneux V et al. Comparison of long-read methods for sequencing and
    assembly of a plant genome. *GigaScience* 9, giaa146 (2020).
    [doi:10.1093/gigascience/giaa146](https://doi.org/10.1093/gigascience/giaa146).
15. Murigneux V et al. Improvements in the sequencing and assembly of plant
    genomes. *GigaByte* (2021).
    [doi:10.46471/gigabyte.24](https://doi.org/10.46471/gigabyte.24).

16. Rhie A, Walenz BP, Koren S, Phillippy AM. Merqury: reference-free quality,
    completeness, and phasing assessment for genome assemblies. *Genome Biology*
    21, 245 (2020). [doi:10.1186/s13059-020-02134-9](https://doi.org/10.1186/s13059-020-02134-9).
17. Mapleson D et al. KAT: a K-mer analysis toolkit to quality control NGS
    datasets and genome assemblies. *Bioinformatics* 33, 574–576 (2017).
    [doi:10.1093/bioinformatics/btw663](https://doi.org/10.1093/bioinformatics/btw663).
18. Mikheenko A et al. TandemTools: mapping long reads and assessing/improving
    assembly quality in extra-long tandem repeats. *Bioinformatics* 36, i75–i83
    (2020). [doi:10.1093/bioinformatics/btaa440](https://doi.org/10.1093/bioinformatics/btaa440).
19. Sweeten A, Schatz MC, Phillippy AM. AniAnn’s: alignment-free annotation of tandem repeat
    arrays using fast average nucleotide identity estimates. *Bioinformatics*
    42, btag581 (2026). [doi:10.1093/bioinformatics/btag581](https://doi.org/10.1093/bioinformatics/btag581).
