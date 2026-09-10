# TandemX reveals residual satellite-repeat under-representation in high-quality plant genome assemblies

**Running title:** Read-based assessment of plant satellite-repeat representation

**Manuscript status.** First complete structural rewrite, 10 September 2026. The
original `paper/manuscript.md` remains the unchanged technical master and evidence
index. This version reorganizes the same frozen analyses as a genomics research
article; it does not replace the underlying evidence archives or promote
development results beyond their validated scope.

## Abstract

Tandem and satellite repeats are major components of plant genomes and often
form long arrays at centromeres, ribosomal DNA loci and other structurally
important regions. Their length, homogeneity and copy number make them difficult
to assemble, and high contiguity alone cannot establish that every repeat family
is represented at the abundance present in the sequenced material. We developed
TandemX to obtain an assembly-independent view of repeat abundance. TandemX
discovers candidate repeat units de novo from long reads, organizes
sequence-supported repeat families, estimates family abundance from diagnostic
k-mers and localizes the same families in an assembly, thereby identifying a
read--assembly abundance deficit consistent with possible under-representation.
Controlled simulations and comparisons with TideHunter, TRF, TRASH and
TideCluster showed that discovery, quantification and localization were reliable
within the tested conditions, while also revealing method- and endpoint-specific
trade-offs; these experiments do not establish universal superiority over any
comparator. We then analysed historical and newer assemblies of *Arabidopsis
thaliana* Ey15-2 and *Macadamia jansenii*. Families predicted from HiFi reads to
be under-represented in the historical assemblies preferentially gained sequence
in the newer high-quality assemblies, but the read estimates also identified
candidate deficits that remained after assembly improvement. Independent
PCR-free Illumina evidence in Ey15-2 and Illumina plus direction-level Oxford
Nanopore evidence in Macadamia supported family-specific residual
under-representation for Ey15-2 families TXF000002 and TXF000154 and Macadamia
family TXF000496. Three other preselected families remained unresolved because
k-mer lengths or sequencing platforms disagreed. Thus, high-quality HiFi
assemblies can retain family-specific satellite-repeat under-representation, and
raw-read-derived abundance provides a complementary measure of assembly
completeness. TandemX reports this evidence with explicit uncertainty; its
continuous deficits are estimated under-representation, not physical missing-base
truth.

## Introduction

Tandemly repeated DNA is a prominent and rapidly evolving component of plant
genomes. Satellite arrays can extend for hundreds of kilobases or megabases and
are frequently concentrated at centromeres, pericentromeric regions, nucleolar
organizers and other chromosome domains with important structural and
cytological properties [1–4]. These repeats are not a single sequence class:
their monomer lengths, sequence divergence, higher-order organization and
genomic copy numbers vary among families, chromosomes, accessions and species.
Accurate representation of repeat families is therefore important both for
describing genome composition and for interpreting centromere evolution,
chromosome structure and repeat-associated variation.

The same properties that make tandem arrays biologically distinctive make them
difficult to reconstruct. Long arrays contain many highly similar copies,
provide few unique anchors and can exceed the length of individual sequencing
molecules. Earlier plant references consequently omitted or compressed many
centromeric satellites and ribosomal DNA arrays [3,5]. PacBio HiFi reads, Oxford
Nanopore ultra-long reads and telomere-to-telomere assembly strategies have
transformed this landscape, resolving repeat-rich regions that were previously
inaccessible and revealing chromosome-specific satellite organization in plants
[3,4,6–8]. Nevertheless, assembly continuity and repeat representation are
different properties. A chromosome-scale or even gapless sequence may contain
an array while still under-representing its copy number, and agreement among
assembly metrics cannot by itself determine how much of a repetitive family was
present in the sequenced material.

This distinction creates a measurement problem. When an assembly contains fewer
copies of a repeat family than expected, the assembly alone does not reveal
whether the difference reflects collapse, biological variation between source
materials, platform or library bias, ambiguity in family definition, or error in
abundance estimation. Comparisons between historical and newer assemblies are
informative when source material is closely matched: repeat families that gain
representation with improved data and assembly methods provide a high-quality
assembly reference proxy for earlier under-representation. A newer assembly is
not, however, absolute copy-number truth. An independent measurement must come
from the reads themselves, and the most informative validation changes library
preparation or sequencing platform while preserving the biological material as
closely as public provenance permits.

Existing methods address important but partly separate parts of this problem.
TRF detects approximate tandem repeats in assembled sequences [9]; TideHunter
detects tandem repeats in noisy long reads [10]; TRASH and TideCluster provide
assembly-oriented repeat annotation and clustering [11,12]; and Satellite Repeat
Finder reconstructs satellite repeat units from sequencing data [13]. These
tools differ in input, repeat representation and output granularity, so no single
comparator supplies a complete test of de novo family discovery, abundance
estimation and assembly representation. What has been missing is a unified,
read-first framework that discovers repeat families without an assembly, measures
their abundance in the read population and then asks how strongly each family is
represented in a candidate assembly.

We developed TandemX to connect these evidence layers. Starting from HiFi reads,
TandemX detects tandem arrays, derives candidate monomers, groups them into
sequence-supported operational families and estimates family abundance using
diagnostic k-mers. When an assembly is available, the same catalogue is localized
to obtain family-level assembly representation. The difference between
read-derived and assembly-derived abundance is reported as a read--assembly
abundance deficit and, when sufficiently large and supported, as possible
under-representation. Family relationships, sparse evidence and conflicting
estimates remain explicit rather than being forced into a single resolved model.

Here we first establish the analytical scope of TandemX with controlled
simulations and task-matched comparisons. We then use public plant long-read
datasets to characterize repeat-family landscapes and focus on two old/new
assembly systems: *A. thaliana* Ey15-2 and *M. jansenii*. Finally, we test the
remaining deficits using independent Illumina evidence and, for Macadamia,
direction-level Nanopore evidence. The results show that assembly improvement
recovers many repeat-rich sequences but does not guarantee complete
representation of every satellite family. They also identify a clear confidence
boundary: abundance estimates can become unstable when diagnostic k-mers sample
different components of a heterogeneous or background-homologous repeat. TandemX
therefore provides an assembly-completeness view that complements contiguity
metrics rather than replacing physical or cytological validation.

## Results

### TandemX links read-derived repeat families with assembly representation

TandemX was designed around a distinction that is not captured by conventional
assembly-contiguity statistics: a repeat family may be present in an assembly
but represented by fewer copies than occur in the sequenced material. The
workflow begins with long reads rather than an assembly (Fig. 1). Periodic
signals identify candidate tandem arrays within individual reads, local
alignment refines array boundaries, and repeated segments are converted to
candidate monomers. Candidates are compared by cyclic, strand-aware sequence
similarity and grouped into operational repeat families. The resulting catalogue
is therefore supported by sequence similarity rather than by shared monomer
length or identifier patterns alone.

Family construction deliberately stops short of evolutionary interpretation.
A representative sequence and its supporting read intervals define each
operational family, while related representatives are retained in a separate
evidence graph. Near-integer relationships between representative lengths can
flag possible period multiples, and other related pairs remain unresolved.
These links can expose harmonic calls, partial units or candidate higher-order
organization, but they do not define ancestry, a rooted hierarchy or a validated
higher-order repeat. This separation prevents a convenient catalogue structure
from being mistaken for a biological model.

For each family, TandemX selects diagnostic canonical k-mers that are exclusive
among catalogue representatives and corrects their counts for multiplicity
within the monomer. Read depth over these diagnostic words provides an estimate
of the genomic abundance represented by the family. When an assembly is
available, the same representative is localized independently and overlapping
intervals are unioned to estimate assembled repeat bases. TandemX then reports
the positive difference between read-derived abundance *R* and assembly-
localized bases *A* (*D* = max[*R* − *A*, 0]). We refer to this
quantity as the read--assembly abundance deficit or estimated
under-representation. It is not a direct estimate of physical missing sequence,
because read sampling, genome-size normalization, family specificity and
assembly localization all contribute uncertainty.

The comparison output separates continuous magnitude from categorical
interpretation. In the controlled and retrospective analyses below, a frozen
assembly/read ratio threshold of 0.6 defined a candidate under-representation
state, but the evidence supporting that state was evaluated independently of the
size of the deficit. Families with inadequate diagnostic support, inconsistent
k-mer lengths or discordant platforms retain explicit warning or `unresolved`
states. Every catalogue family also retains its analysis fate, so a family
excluded by the biological denominator is not conflated with a technical
failure.

HiFi reads were used to build the TandemX catalogue and its initial abundance
prediction. Illumina and ONT data were introduced only later as orthogonal
tests. They are not presented here as additional TandemX input modes: Illumina
provides an independent diagnostic-k-mer abundance estimate, whereas ONT is used
for direction-level long-read support. This design makes the main biological
question explicit—whether a HiFi-derived family deficit persists when examined
with evidence that did not generate the newer assembly.

The software now presents these outputs through an offline family report. A
single run command accepts HiFi reads, an optional assembly and an output
directory; detailed settings remain available through a configuration file.
The report links family abundance, assembly representation, confidence labels
and warnings to the underlying catalogue and candidate relationship graph.
Editable vector exports and source tables support inspection and reuse. These
interface additions do not change the frozen analytical thresholds or provide
independent validation of repeat-family architecture.

### TandemX recovers tandem-repeat families and abundance across controlled benchmarks

We evaluated discovery, abundance estimation and assembly localization as
separate tasks before combining them (Fig. 2). In a factorial 10-Mb simulation
containing 55 planted families, TandemX recovered all founders from both clean
and high-error 5× long-read datasets. On a broader validation set spanning
repeat period, sequence divergence, indel rate, array number and negative
controls, the final discovery path retained high array and family recovery and
made no calls in the three negative-control classes. Closely related monomers
were evaluated by cyclic sequence identity rather than exact representative
length; this recovered three of three related planted families in the final
sequence-clustering validation. These tests establish recovery within the
simulated sequence models, not complete recall for biological satellites.

Quantification was evaluated with the planted founder catalogue supplied, so
copy-number error could be separated from discovery error. Across three
independent 10-Mb validation genomes and 1,485 family conditions, a frozen
depth-gated normalization reduced mean absolute relative error from 0.409 to
0.363. The change improved 661 conditions, left 495 unchanged and worsened 329,
showing that a lower aggregate error did not imply uniform improvement. A
multi-k model was most useful in the deliberately difficult 20× condition with
2% unit divergence and high read error, where mean absolute relative error fell
from 56.6% at k=21 to 13.5%. The benefit was conditional: diagnostic support was
often unavailable at low coverage, and some families became less accurate.
Intervals derived from joint read-level moments covered the planted value in
401 of 424 supported 20× conditions, but only 424 of 495 conditions yielded an
interval. We therefore report both conditional coverage and missingness rather
than describing the interval as generally calibrated.

Assembly localization was challenged with divergent monomers and interrupted
arrays. An exact-anchor localizer initially lost most array bases at 3–5% unit
divergence. Allowing nearby anchors to bridge no more than one monomer length
restored full-assembly mean recall to 0.981 on fresh simulated genomes while
maintaining mean positive-assembly precision of 0.9996; the weakest 5%-divergent,
three-segment stratum had recall of 0.949. The downstream binary comparison
remained more sensitive to abundance error than localization alone. In the final
depth-gated held-out test, sensitivity increased from 0.595 to 0.654 without
increasing the 0.016 false-positive rate, and precision was 0.983. Earlier
domain-shift and classifier experiments that increased false-positive rates are
retained in Supplementary Results and define why the final method reports
uncertainty rather than a universal collapse call.

Comparator results depended on the matched task and endpoint. In the 55-family
read simulation, TandemX, TRF and TideHunter all recovered every planted founder.
For high-error reads, one-to-one eligible-array precision was 1.000, 0.699 and
0.920, respectively, whereas base-union precision exceeded 0.996 for every
method; overlapping or harmonic calls therefore affected interval counts more
than the underlying covered bases. TideHunter was generally faster in the
tested long-read comparisons, while memory and runtime rankings varied with
input scale. On the assembly side, TRASH2 recovered all 55 planted families and
localized their native unit intervals with high base precision and recall once
its output hierarchy was interpreted correctly. TideCluster also achieved high
base-union accuracy in successful planted-truth runs, but its clustering stage
approached the available memory ceiling in the tested container. Failed or
unstarted comparator cells remain unavailable rather than being assigned zero
accuracy.

Together, the controlled analyses show that TandemX can recover repeat families,
estimate abundance and localize arrays across the tested conditions. They also
show why no global superiority claim is warranted: alternative tools can match
family recovery, exceed TandemX in speed, or provide strong assembly annotation,
and the ranking changes with the metric being measured. The development
failures, parameter-selection history, resource profiles and frozen validation
records are provided in Supplementary Results and the evidence archive rather
than used as the organizing narrative of the main text.

### Read-derived abundance identifies repeat families under-represented in historical plant assemblies

We first established that the workflow could be applied to plant long-read data
spanning a broad range of genome sizes and repeat landscapes. Ten complete HiFi
libraries from eight reported species passed source-checksum and streaming
FASTQ validation, together comprising 14,937,608 reads and 262.7 Gb of sequence
(Fig. 3; Table 1). These data include Arabidopsis, rice, maize, barley, rye,
wheat, oat and wild soybean. Their catalogues differed widely in family number,
monomer length, read support and estimated abundance. Because several datasets
represent pooled plants or technical libraries and because matched repeat truth
is unavailable, these contrasts are descriptive family landscapes rather than
cross-species accuracy estimates or biological replicates.

The central biological test used old and new assemblies generated from closely
matched plant material. In *A. thaliana* Ey15-2, the source study compared
Bionano-scaffolded CLR-Canu and HiFi-Hifiasm assemblies from the same accession
and sample designation [4]. We constructed one catalogue from the complete public
HiFi library and froze family eligibility and comparison thresholds before
examining old/new localization. Nineteen families had at least 15 kb represented
in the newer assembly and entered the primary reference-proxy analysis. Eight
had an old/new assembly ratio below 0.6, indicating that representation had
increased substantially in the HiFi assembly. TandemX independently identified
all eight from the HiFi-read abundance relative to the historical assembly and
did not label the other 11 families. At the more permissive 5-kb denominator,
the larger and more difficult set exposed one missed and two additional calls,
whereas only six families remained at 50 kb (Fig. 4).

The direction of the Ey15-2 result was stronger than its quantitative agreement.
The total HiFi-versus-old abundance deficit was 840,658 bp, compared with a
973,959-bp new-versus-old assembly gain, but family-level magnitudes differed:
Pearson correlation was 0.589 and the mean absolute difference was 45,551 bp.
Assembly-to-assembly alignments provided consistent context—the eight proxy-
positive families had lower coverage in the historical assembly than the other
families—but neither the newer assembly nor the alignment is independent
copy-number truth. The analysis therefore supports the ability to prioritize
families that preferentially gain representation as assemblies improve, not the
ability to infer exact missing bases.

We repeated the design in *M. jansenii*, a phylogenetically distant species with
an earlier CLR Falcon-Unzip primary assembly and a newer HiFi IPA assembly. The
update study described the HiFi material as the same sample used for the earlier
comparison, although different archival BioSamples mean that identical DNA
extraction cannot be established [14,15]. The newer assembly was more contiguous
(284 versus 762 primary contigs), but slightly shorter overall. From 1,227
read-derived families, 43 met the frozen 15-kb newer-assembly eligibility
criterion. Two were strongly under-represented in the historical assembly and
41 were in other states; the HiFi read comparison identified both and did not
label the other 41. This agreement was preserved when read depth was normalized
by the reported 28-Gb yield rather than the public FASTQ total.

As in Ey15-2, the Macadamia binary result did not imply accurate magnitude. The
HiFi-versus-old abundance deficit was 2.12 Mb under the primary normalization,
whereas the new-versus-old assembly gain among the eligible families was
0.159 Mb; the mean absolute family difference was 50.3 kb and Pearson
correlation was 0.535. Primary alignment coverage was also lower for the two
proxy-positive families, but old and new contig identifiers were not directly
comparable. Thus, two independent plant systems showed the same qualitative
pattern: families predicted from raw reads to be depleted in a historical
assembly preferentially gained representation in a newer assembly, while the
continuous deficit remained an estimate rather than physical missing-sequence
truth.

Other old/new assembly candidates were screened but not converted into
favourable results when their evidence design was inadequate. Cassava lacked the
exact historical assembly sequence, tomato retained donor-mismatch risk, and a
maize B73-Ab10 update remained incomplete specifically within the tandem arrays
needed as a reference proxy. Their exclusion criteria and all unexecuted or
failed acquisition paths are retained in Supplementary Results. Ey15-2 and
Macadamia were carried forward because they provided the strongest available
old/new source relationships, not because their newer assemblies were assumed
to be complete.

### High-quality HiFi assemblies retain family-specific satellite-repeat under-representation

The old/new analysis established that read-derived abundance could identify
families under-represented in historical assemblies, but it did not answer
whether candidate deficits in the newer assemblies reflected residual
under-representation or an upward bias in the HiFi estimator. We therefore
froze six newer-assembly deficit candidates before examining other sequencing
platforms: TXF000002, TXF000154 and TXF001517 in Ey15-2, and TXF000496,
TXF000563 and TXF000695 in Macadamia (Fig. 5). All other eligible families were
retained as context rows. Family definitions, the 15-kb denominator, diagnostic
k-mer rules and the 0.6 ratio threshold were unchanged.

For Ey15-2, PCR-free Illumina reads provided evidence from an independent DNA
extraction and library preparation from the same ground-tissue pool used for
high-molecular-weight DNA. Family abundance was estimated at k=21 and k=31
using family-exclusive diagnostic k-mers and empirical single-copy controls.
TXF000002 had an Illumina abundance estimate of 519,994–555,544 bp, compared
with 49,758 bp localized in the newer assembly. TXF000154 had an Illumina
estimate of 33,895–44,496 bp, compared with 18,297 bp in the assembly. Both
families remained below the frozen assembly/read ratio threshold at both k
values and therefore had stable orthogonal support for residual
under-representation. The HiFi estimate for TXF000154 was nevertheless
1.52–1.99-fold above the Illumina estimate, illustrating that agreement in
direction does not imply agreement in abundance magnitude.

Macadamia provided both Illumina and PromethION evidence from the accession-1005
study. TXF000496 showed the clearest cross-platform result. Its Illumina
estimate was 1.586–1.604 Mb, whereas only 55,562 bp was localized in the newer
assembly. Direction-level ONT occupancy, evaluated across a 616–780-Mb
genome-size sensitivity, corresponded to 0.620–0.785 Mb and supported the same
assembly-deficit direction. ONT was used only to test whether substantially more
family sequence was present in long reads than in the assembly; it was not
treated as precise physical copy-number truth.

The independent data therefore support family-specific residual
under-representation for Ey15-2 TXF000002 and TXF000154 and Macadamia
TXF000496. None of the 16 Ey15-2 or 40 Macadamia eligible context families
crossed the orthogonal deficit rule at either k, providing a check against a
global shift in the Illumina normalization. The selected outcome—three supported
families and three unresolved families—is not a 50% accuracy estimate. It is a
family-level evidence statement about six candidates chosen before orthogonal
inspection. Likewise, the difference between independent read abundance and
newer-assembly localized bases is a read--assembly abundance deficit, not a
measurement of physically missing DNA.

### Discordant repeat families define the confidence boundary of abundance estimation

The other three candidates did not satisfy a stable cross-k and cross-platform
interpretation (Fig. 6). Ey15-2 TXF001517 changed from an inconclusive k=21
estimate to a k=31 estimate compatible with quantification bias, leaving the
cross-k result unresolved. Macadamia TXF000563 lay close to the decision
boundary: k=21 placed the assembly/read ratio just below 0.6, k=31 placed it
above the threshold, and ONT did not support a residual deficit. It therefore
also remained unresolved rather than being assigned according to the more
favourable k value.

TXF000695 exposed a more fundamental abundance boundary. Its 467-bp
representative yielded 467 diagnostic k-mers at both k=21 and k=31. Thus, the
instability was not caused by a smaller nominal diagnostic set at the longer k.
Of these words, 381 at k=21 and 299 at k=31 were observed at least twice. The
median Illumina depth nevertheless changed from 77,091 to 101, a 763-fold
difference, and the corresponding abundance estimate changed from 1,104,859 to
1,649 bp, a 670-fold difference. The result also did not arise from one extreme
word: the median absolute deviation was 18,812 at k=21 and 101 at k=31.

Inspection of the complete diagnostic-word distribution showed two sharply
separated depth components. At k=21, 259 of 467 words occupied the high-depth
component between approximately 50,000 and 100,000 counts, placing the median
inside that component. At k=31, only 202 words remained in the high-depth
component and 265 fell below 10,000 counts, including 168 absent or singleton
words; the median consequently fell to the low-depth component. This median
cliff explains the numerical discontinuity without establishing its biological
cause. A plausible interpretation is that shorter diagnostic words primarily
sample a locally conserved, high-abundance segment shared with related genomic
sequence rather than the abundance of the complete 467-bp repeat unit. That
interpretation remains a hypothesis because the current catalogue and read data
do not uniquely distinguish mosaic family structure, background homology and
other sources of k-mer sharing.

Independent mapping supported the low-abundance direction: raw HiFi mapping
occupancy corresponded to approximately 960 bp, and calibrated ONT occupancy
gave 1,344–1,702 bp across the genome-size sensitivity, close to the k=31
estimate and far below both the k=21 and frozen HiFi k-mer estimates. The newer
assembly contained 20,424 bp. These observations argue against calling
TXF000695 a stable residual-deficit family, but they also do not validate a
general quantification-bias mechanism. We therefore retain `unresolved` as the
only supported conclusion.

No candidate satisfied the predeclared stable cross-k quantification-bias rule.
The orthogonal gate therefore did not trigger a change to the frozen abundance
estimator.

This example motivates reporting-level quality-control annotations rather than
a new tuned classifier. For each family, cross-k fold difference, fraction of
diagnostic words observed, and MAD/median can accompany heuristic warnings such
as `unstable_across_k`, `depth_median_boundary_sensitive` and
`reliability=low`. For TXF000695, the observed fractions were 0.816 at k=21 and
0.640 at k=31, while MAD/median changed from 0.244 to 1.000. These annotations
describe why the estimate is fragile; they are not independently validated
decision rules, do not modify the core estimator and must not be interpreted as
proof that one k value is correct.

### TandemX provides an assembly-completeness view complementary to contiguity metrics

The combined analyses support a role for read-derived family abundance that is
distinct from conventional assembly evaluation. In Ey15-2 and Macadamia,
families depleted from historical assemblies preferentially gained
representation in newer assemblies, demonstrating that the read--assembly
comparison tracks biologically relevant improvements in repeat-rich sequence.
Orthogonal reads then showed that three family-specific deficits persisted in
the newer HiFi assemblies. These findings do not imply that the assemblies are
globally incomplete or that contiguity metrics are uninformative. Rather, they
show that a high N50, small contig number or nominally gapless chromosome cannot
by itself establish complete representation of every tandem-repeat family.

TandemX provides an external view by asking whether the abundance implied by raw
reads is compatible with the amount localized in an assembly. Agreement across
k values and independent platforms can prioritize families for assembly review,
centromere analysis or cytological follow-up. Disagreement is equally
informative: it identifies repeat families for which sequence homology,
diagnostic-word composition or platform effects prevent a stable conclusion.
The appropriate output is therefore a family-level evidence state with explicit
uncertainty, not a single genome-wide completeness score.

### A bounded recovery test does not establish a missing-sequence repair

We tested whether the two Ey15-2 families with stable orthogonal abundance
support could yield a useful sequence recovery from the historical assembly
and the original HiFi reads. The frozen historical localization contained six
intervals for TXF000002 and none for TXF000154. Of 48 predefined flank windows,
20 met the assembly-relative uniqueness criterion. Full-read recruitment
produced one unpolished spanning-read candidate at a TXF000002 interval,
supported by 26 concordant dual-anchor reads. The other five intervals lacked
a usable anchor pair or sufficient unambiguous spanning-read support;
TXF000154 retained an explicit no-localized-locus outcome.

The candidate was 5,681 bp between the selected anchor interiors, compared
with 5,679 bp in the historical assembly. Its span included non-repeat
sequence. After candidate sequences and tables were locked, the same anchor
pair mapped uniquely around a 5,679-bp interval in the newer assembly. The
frozen full-catalogue localizer identified 679 repeat bp in both assembly
intervals and 644 bp in the unpolished candidate. Thus this candidate did not
increase measured repeat representation or improve agreement with the newer
reference proxy. After excluding the candidate source read, 97 of 120
flank-recruited reads aligned across at least 90% of each old, candidate and
newer span at at least 98% identity. This supports the observed span but does
not distinguish a recovered repeat expansion. These results do not establish a recovery method and do not
justify further algorithm expansion for this release. They also do not show
that the families are intrinsically unresolvable: the test examined specified
historical loci, flank offsets and direct read spans, rather than every
possible locus or assembly path.

## Discussion

Long-read sequencing and telomere-to-telomere assembly have changed the study of
plant repetitive DNA, but they have not made repeat completeness a solved
problem. The central result of this study is that family-specific satellite-
repeat under-representation can persist in highly contiguous HiFi assemblies.
This conclusion is deliberately narrower than stating that the assemblies are
inaccurate overall. Ey15-2 and Macadamia both contain repeat families for which
the newer assembly represents much more sequence than its historical
predecessor, yet independent reads still support more family sequence than is
localized in the newer assembly. Contiguity and repeat-family representation
must therefore be evaluated as related but distinct dimensions of assembly
quality.

The main contribution of TandemX is not a universal speed or sensitivity
advantage over existing repeat software. Its contribution is to connect three
measurements that are usually separated: de novo family discovery in raw long
reads, family abundance in the read population and representation of the same
family in an assembly. That connection makes it possible to ask an
assembly-completeness question without defining the assembly as truth. TRF,
TideHunter, TRASH, TideCluster and SRF remain useful, task-specific alternatives
or complements. Our controlled comparisons show that several recover the same
planted families, TideHunter is often faster, and assembly-oriented tools can
provide highly accurate interval annotation. TandemX instead emphasizes a
family-level read--assembly comparison with explicit uncertainty and a direct
path to orthogonal validation.

The old/new assembly analyses support this framing. In Ey15-2, all eight
families that gained substantial representation from CLR-Canu to HiFi-Hifiasm
were independently prioritized by their HiFi-read abundance, while the other 11
primary families were not. Macadamia provided the same qualitative result in a
distant plant lineage: both families that gained representation from the CLR to
HiFi assembly were read-predicted, and the other 41 primary families were not.
The strength of this agreement lies in its direction and family ranking. The
continuous read deficit did not closely reproduce the size of the assembly
gain, especially in Macadamia. A newer assembly is also influenced by the same
HiFi data used for the initial TandemX prediction. Old/new concordance is
therefore a high-quality assembly reference proxy, not an independent estimate
of genomic copy number.

The platform-independent tests address that dependence. Ey15-2 PCR-free
Illumina data changed DNA extraction, library construction and sequencing
technology while preserving the source tissue pool. Macadamia Illumina data
provided an independent diagnostic-k-mer estimate, and Macadamia ONT reads
added long-read directional evidence. Agreement at two Illumina k values
supported residual under-representation for Ey15-2 TXF000002 and TXF000154;
agreement between Illumina and ONT supported Macadamia TXF000496. These are
three family-specific observations, not an estimate of a genome-wide residual-
collapse rate. ONT occupancy is especially useful as a change in measurement
mode, but its platform error, mapping efficiency and genome-size dependence
preclude treating it as exact copy-number truth.

The discordant families are part of the principal result rather than failed
examples to be hidden. TXF001517 and TXF000563 show that conclusions near the
0.6 boundary can depend on k and platform. TXF000695 demonstrates a sharper
failure mode: two equally large diagnostic sets sampled a bimodal depth
distribution on different sides of a median boundary and produced abundance
estimates separated by approximately three orders of magnitude. Mapping-based
evidence favoured the lower-abundance direction, but the present data cannot
determine whether the high-depth k=21 component arose from a mosaic repeat,
local sequence conservation shared across families, background homology or a
combination of these effects. Retaining this family as `unresolved` is therefore
more informative than forcing either residual under-representation or confirmed
quantification bias. Cross-k stability, observed-word fraction and dispersion
should be reported as heuristic reliability indicators, not promoted to a new
classifier without independent validation.

This family-level view has several applications. It can prioritize repeat
families for targeted reassembly or manual inspection when chromosome-scale
metrics appear satisfactory; identify centromeric or ribosomal arrays that
merit cytological validation; provide a quantitative context for satellite
evolution across accessions or species; and help select abundant, specific
candidate sequences for FISH probe development. Cohort comparisons are also
possible when extraction, ploidy, genome-size normalization and family
correspondence are controlled. In each application, TandemX is best used as a
screen for where assembly and read evidence disagree, followed by an independent
test appropriate to the biological question.

Several limitations define the current boundary. First, no independent
physical copy-number measurement is available for the reported families;
read-derived abundance and read--assembly deficit must not be converted to
missing-base truth. Second, the old/new systems represent two species and small
family denominators, and public metadata do not prove identical DNA extraction
for all Macadamia libraries. Third, the orthogonal set was deliberately selected
from six pre-existing candidates and cannot be interpreted as population
accuracy or prevalence. Fourth, diagnostic k-mers can be affected by family
definition, sequence divergence, low complexity, background homology and
platform/library bias. Fifth, the controlled simulations simplify ploidy,
population heterogeneity, long-range array organization and empirical error
structure, and real-read interval truth remains unavailable. Sixth, ONT mapping
was designed for direction-level support only, whereas the Illumina estimates
depend on empirical controls and k-mer specificity. Finally, TandemX estimates
family abundance and localizes supported arrays; it does not claim to reconstruct
every megabase-scale array or resolve every higher-order organization. Within
these limits, the data establish the main biological point: raw-read-derived
family abundance provides information about repetitive-sequence completeness
that assembly contiguity alone cannot supply.

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
state was defined by an assembly/read ratio below 0.6 in the frozen experiments.
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
rule below estimated haploid depth 2 and used the frozen multi-k blend at higher
depth. No model was refitted in the old/new or orthogonal plant analyses.

### Assembly localization and comparison

Catalogue representatives were localized to assemblies by exact diagnostic
anchors followed by bounded extension and bridging. Overlapping intervals were
unioned within family before their bases were summed. The final simulated
localizer allowed anchor gaps of no more than one monomer length and applied the
frozen identity criterion. Localization truth was scored as base-union recall
and precision, separately from family recovery and binary comparison.

For a family with read-derived abundance *R* and assembly representation *A*,
TandemX reports max[*R* − *A*, 0] and the ratio *A*/*R* when *R* is positive. The frozen
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
the frozen HiFi result indicated a deficit and exceeded orthogonal abundance by
at least 1.5-fold. All other outcomes were unresolved. Binary interpretation was
reported separately from deficit magnitude.

### Cross-k diagnostic assessment of TXF000695

The complete diagnostic-k-mer count distributions were retained for the
preselected Macadamia families. For TXF000695, the 467 family k-mers at each k
were joined to the exact sequence counts; words absent from the counter output
were assigned the same frozen 0/1 midpoint used by the abundance analysis. We
reported diagnostic-set size, fraction observed at least twice, median, median
absolute deviation, MAD/median and the cross-k fold differences in median depth
and estimated abundance. Counts were also summarized by fixed depth ranges to
determine whether the cross-k reversal reflected one outlier or a redistribution
of the diagnostic set. The labels `unstable_across_k`,
`depth_median_boundary_sensitive` and `reliability=low` are descriptive QC
annotations for this manuscript. They were not used to change family selection,
the estimator, the 0.6 threshold or the frozen interpretation rules.

### Bounded targeted recovery and post-lock proxy assessment

The recovery proof of concept used only the historical Ey15-2 assembly, its
frozen family localization, the complete frozen catalogue and original HiFi
library ERR8666125 during candidate generation. TXF000002 and TXF000154 were
selected from the previously completed orthogonal abundance gate. Their
family-level newer-assembly abundance results were already known; the design
therefore separates recovery inputs from later proxy inspection and is not a
wholly unobserved validation set. We examined 2-kb flanks at offsets of 0, 2, 5
and 10 kb from each historical interval. Eligible flanks required an expected
historical alignment with at least 90% coverage and 98% identity, without a
competing alignment meeting 90% coverage and 95% identity. The nearest eligible
flank on each side was selected using historical evidence alone.

Minimap2 recruitment retained repeat-supporting alignments of at least 500 bp
and 85% identity, and flank alignments covering at least 90% of the flank at
98% identity. Candidate read spans required distinct read IDs, unambiguous
placements, consistent strand and order, and at least three dual-anchor reads
whose length range did not exceed max(100 bp, 1% of their median). The observed
read span nearest the median supplied an unpolished candidate; it was not a
consensus or a resolved assembly. The initial 30-kb recruitment-template run
was operator-stopped for resource cost. A separately recorded 1-kb minimum
template continuation tested family presence while retaining the same
biological filters. A streaming evidence collector preserved the complete
alignment set after a smaller collector limit was reached. Both technical
outcomes remain archived separately from biological results.

Candidate hashes were sealed before querying newer-assembly sequence. The
selected historical flanks were then mapped to the newer assembly, and the
original, candidate and proxy spans were evaluated using the frozen complete
catalogue, k = 21 and the existing iid-base localization setting at 0.9. This
localizer is a consistent representation measure, not independent physical
copy truth. Candidate spans and localized repeat bases were recorded
separately. For read remapping, we extracted the 120 flank-recruited reads other than the
candidate source read, verified the complete FASTQ digest and mapped all three
comparison spans to this subset. Supporting read IDs were counted once at
90% query coverage and 98% alignment identity. This same-library alignment
check is distinct from orthogonal abundance validation. Original assemblies
were not modified.

### User workflow and report provenance

The run interface detects supported FASTA/FASTQ and gzip inputs, validates
configuration values and reuses completed stages only when their input and
configuration fingerprints match. The documented k = 21 default is an operating
setting, not a data-adaptive optimum. If the user supplies an assembly but no
genome-size denominator, its streamed sequence length is explicitly labelled a
provisional normalization denominator. Without an assembly, genome size or
sufficient depth information, discovery remains available but absolute abundance
is withheld. Family summaries retain separate discovery, abundance and
comparison confidence fields; missing evidence is distinct from zero
representation. Report products include source hashes, editable SVG and PDF,
PNG previews, source tables and export receipts. Candidate relationship graphs
preserve the existing pairwise evidence and do not infer monomer order.

### Statistical analysis and reproducibility

The family was the unit of the old/new and orthogonal comparisons. Binary
sensitivity, false-positive rate and precision were kept separate from
continuous abundance error. Correlations between read-versus-old deficits and
new-versus-old assembly gains were descriptive and did not treat the latter as
truth. Simulation summaries distinguished independent source genomes from
dependent coverage/error conditions and technical timing repetitions.

Analyses were frozen or specified before the relevant validation data were
examined where stated. Independent implementations recomputed the Ey15-2 and
Macadamia family tables and the successful TideCluster accuracy endpoints.
Exact commands, software identities, validation splits, run fates, file hashes,
resource profiles and rejected analyses are retained in the technical master
and evidence archives. They provide reproducibility support but are not treated
as biological results.

## Data and code availability

The TandemX source code and committed evidence are available at
[github.com/CongyangY/TandemX](https://github.com/CongyangY/TandemX). The
technical master manuscript, complete evidence index, source-data tables,
accepted and rejected figure versions, and machine-readable run fates are
retained under `paper/` in the repository. Public read and assembly accessions
are given in Methods and the evidence manifests. Large raw FASTQ/FASTA files,
mapping outputs and scratch databases are regenerated from the recorded public
sources and are not versioned in Git. A versioned software release, Bioconda
recipe and Zenodo deposition remain in preparation; no DOI or final release
identifier is claimed in this draft.

## Table legend

**Table 1. Public plant HiFi libraries used for input validation and descriptive
real-data analyses.** The table reports species and material labels, run
accessions, complete-file read and base counts, read-length summaries, GC and N
fractions, duplicate-read-ID checks, source-file identities and provenance
warnings for ten libraries from eight reported plant species. These are
validated input files, not eight biological accuracy experiments. Source:
`paper/tables/input_cohort.tsv`.

## Figure legends and construction plan

### Figure 1. Biological problem and TandemX workflow

**Biological problem and TandemX workflow.** **A,** A long, homogeneous plant
satellite array can be present but under-represented in an otherwise highly
contiguous assembly. **B,** HiFi reads are scanned for tandem arrays and
sequence-supported monomers are organized into an operational family catalogue.
**C,** Diagnostic-k-mer depth estimates read-derived family abundance, whereas
the same representative is localized independently in an assembly. **D,** The
comparison produces represented, possible under-representation and unresolved
states; the continuous read--assembly abundance deficit is not physical missing
bp. **E,** HiFi data generate the TandemX prediction, after which independent
Illumina k-mer abundance and ONT direction can provide orthogonal support. The
latter are validation measurements rather than additional TandemX input modes.
This figure requires a new editable-vector schematic and a compact definitions
source table.

### Figure 2. Integrated analytical validation across controlled benchmarks

**Integrated analytical validation across controlled benchmarks.** **A,** De
novo read-level array and family recovery across repeat periods, divergence and
read-error conditions. **B,** Abundance error across coverage and error strata,
including both improved and worsened family conditions. **C,** Assembly
localization recall and precision, with the divergent and interrupted-array
boundary shown explicitly. **D,** Binary under-representation sensitivity,
false-positive rate and precision, displayed separately from continuous
abundance error. **E,** Task-matched read comparisons with TideHunter and TRF
and assembly comparisons with TRASH/TRASH2 and TideCluster. **F,** Summary of
accuracy, runtime and memory trade-offs. The figure combines only the final
validation and principal adverse boundary; development iterations remain in
Supplementary Figures. Source panels are assembled from the current controlled-
benchmark evidence tables without selecting only favourable strata.

### Figure 3. Repeat-family landscapes across diverse plant long-read datasets

**Repeat-family landscapes across diverse plant long-read datasets.** **A,**
Species, material, genome-size denominator and validated sequence yield for ten
libraries from eight reported plant species. **B,** Number and abundance
distribution of operational repeat families at matched nested sampling scales.
**C,** Monomer-length and sequence-complexity landscape. **D,** Descriptive
fraction of read-derived family abundance represented in available assemblies.
**E,** Representative family profiles from small and large plant genomes.
**F,** Provenance and denominator map identifying pooled material, technical
batches and unresolved donor/reference matches. Called-base fraction is not an
accuracy endpoint, and the libraries are not treated as interchangeable
biological replicates. A consolidated source table must be generated from the
existing cohort-QC and real-diagnostic tables before this figure is rendered.

### Figure 4. Old/new assemblies provide reference-proxy evidence for historical under-representation

**Old/new assemblies provide reference-proxy evidence for historical
under-representation.** **A–B,** Ey15-2 and Macadamia family representation in
historical and newer assemblies, highlighting families independently predicted
from HiFi reads to be depleted in the historical assembly. **C,** HiFi-derived
abundance versus historical-assembly representation. **D,** New-versus-old
assembly gain compared with HiFi-versus-old abundance deficit; binary agreement
and magnitude disagreement are displayed separately. **E,** Sensitivity to the
5-, 15- and 50-kb newer-assembly eligibility denominators. **F,** Assembly-to-
assembly aligned-query coverage as explanatory context. The newer assemblies
are high-quality reference proxies, not absolute or independent copy-number
truth. The accepted Ey15-2 panels and Macadamia family tables will be rebuilt as
one two-species composite rather than retaining the old Ey15-only main figure.

### Figure 5. Orthogonal reads support residual under-representation in three repeat families

**Orthogonal reads support residual under-representation in three repeat
families.** **A,** Frozen HiFi abundance and newer-assembly representation for
six candidates selected before orthogonal inspection. **B,** Independent
Illumina k=21 and k=31 estimates. **C,** Direction-level Macadamia ONT support
across the 616–780-Mb genome-size sensitivity. **D,** Stable support for Ey15-2
TXF000002 and TXF000154. **E,** Cross-platform support for Macadamia TXF000496.
**F,** Final evidence matrix showing three supported and three unresolved
families. The split is not an accuracy rate, and all continuous differences are
read--assembly abundance deficits rather than physical missing-base truth. The
current accepted orthogonal-validation figure supplies the source data but will
be expanded to six focused panels.

### Figure 6. Discordant families and the confidence boundary of read-derived abundance

**Discordant families define the confidence boundary of read-derived
abundance.** **A,** Cross-k and cross-platform outcomes for TXF001517,
TXF000563 and TXF000695. **B,** Complete TXF000695 diagnostic-k-mer depth
distributions at k=21 and k=31, showing the high- and low-depth components and
the shift of the median between them. **C,** Diagnostic-set size and observed
fraction: both k values contain 467 diagnostic words, of which 381 and 299 were
observed at least twice. **D,** Median depth differs by approximately 763-fold
and estimated abundance by 670-fold. **E,** ONT and HiFi mapping occupancy are
closer to the low-abundance direction, while the mechanism remains unresolved.
**F,** Reporting-level QC summary showing cross-k fold, observed fraction,
MAD/median and the heuristic labels `unstable_across_k`,
`depth_median_boundary_sensitive` and `reliability=low`. These warnings are not
a validated classifier and do not alter the estimator or family interpretation.
The derived QC values and their raw-input mapping are retained in
`paper/0910/source_data/figure6_txf000695_qc.tsv` and its accompanying notes.

## Supplementary migration map

The technical master, all source data and every adverse or incomplete result are
retained. Commit identifiers, CI run identifiers, exact seeds, hashes,
byte-identity checks, failed transfers, operator-stopped runs, build histories,
development versions and evidence E1–E42 are moved out of the narrative Results
and remain available through `paper/manuscript.md`, Supplementary Methods,
Supplementary Results, Source Data and `paper/evidence/`. No technical record is
deleted or relabelled as successful evidence.

| Technical-master material | Destination in the new article package |
| --- | --- |
| Original Figures 1–2 (multi-k and joint-read uncertainty) | Condensed in new Figure 2; full panels and all unavailable intervals retained as Supplementary Figures S1–S2 |
| Original Figures 3–7 (conditional comparison, domain shift, localizer and classifier sequence) | Final validation and principal boundary condensed in new Figure 2; complete development, failed held-out and final depth-gated results retained as Supplementary Figures S3–S7 |
| Original Figures 8–11 (cascade, calibration and performance gates) | Essential task-matched trade-offs condensed in new Figure 2; complete gates, runtime distributions and adverse cells retained as Supplementary Figures S8–S11 |
| Original Figure 12 (Ey15-2 proxy analysis) and Macadamia proxy evidence | Rebuilt together as new Figure 4; full species-specific panels retained in Supplementary Figures |
| Original Figure 13 (orthogonal validation) | Source for new Figure 5; accepted and rejected layouts remain archived |
| New TXF000695 depth-distribution analysis | New Figure 6; exact diagnostic-word table and summary become Source Data |
| Current Supplementary Figures S1–S6 | Renumbered after the above transfer without dropping input-QC, real-diagnostic, TideCluster or failed-run panels |
| Original Supplementary Tables S1–S39 and Evidence E1–E42 | Retained in full and renumbered only after main-text figures and tables are finalized |

The main text will cite the biological or methodological conclusion; the
supplement will carry the exact frozen configuration, complete unfavourable
rows, run fate, reproduction record and resource context. This separation
changes article hierarchy, not evidence content.

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
    [GitHub](https://github.com/kavonrtep/TideCluster). Final bibliographic
    metadata require reference-manager verification.
13. Zhang Y, Chu J, Cheng H, Li H. De novo reconstruction of satellite repeat
    units from sequence data. *Genome Research* 33, 1994–2001 (2023).
    [doi:10.1101/gr.278005.123](https://doi.org/10.1101/gr.278005.123).
14. Murigneux V et al. Comparison of long-read methods for sequencing and
    assembly of a plant genome. *GigaScience* 9, giaa146 (2020).
    [doi:10.1093/gigascience/giaa146](https://doi.org/10.1093/gigascience/giaa146).
15. Murigneux V et al. Improvements in the sequencing and assembly of plant
    genomes. *GigaByte* (2021).
    [doi:10.46471/gigabyte.24](https://doi.org/10.46471/gigabyte.24).

Reference 12 and complete dataset-specific metadata remain to be verified in a
reference manager before submission. No author, funding or competing-interest
information is inferred in this draft.
