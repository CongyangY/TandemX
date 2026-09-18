# Frozen read-span development endpoints, 2026-09-18

The existing Ey15-2 and Macadamia nine-case sets are **development** data,
frozen in
`benchmarks/controlled_collapse/read_span_validation_v1_20260918/protocol.json`
before any new independent source is selected. Their prior
5% threshold, 1,024-bp natural flank, >=99% source alignment, >=1-kb
alignment flank, <=10-bp flank/CIGAR disagreement and >=3 distinct eligible
read/record requirements are unchanged. All edits from a source share one
array and the same original reads; there are two source clusters, not 18
independent biological validations.

The descriptive continuous-endpoint archive is
`benchmarks/controlled_collapse/read_span_validation_v1_20260918/development_continuous_endpoints.json`.
It computes signed read-minus-
assembly span deficit, absolute bp error versus the **injected deletion**,
and relative error only for nonzero deletions. This tests read-span arithmetic
against deliberately edited assembly contexts, not physical family copy
number or biological collapse magnitude.

| Development source | Eligible reads/records | Injected deletions marked discordant | Intact false discordance | Smallest tested deletion | Mean absolute bp error across nine edits |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ey15-2 | 7 HiFi reads | 8/8 | 0/1 | 811 bp | 0 bp |
| Macadamia | 7 SRA records; ZMW unverified | 8/8 | 0/1 | 788 bp | 1 bp |

The near-exact magnitude errors follow from unchanged read spans and exact
assembly deletions on one nominated interval. They are **not** evidence that
TandemX reconstructs unknown genome-wide missing bases. A formal minimum
detectable deletion below 788/811 bp was not measured. Localization error
is unavailable because source loci and flank coordinates were supplied to
the scorer. The seven selected reads per locus do not support an independent
read-length effect analysis. No preregistered 5×/10×/20×/30× whole-source
downsample has been run; coverage response is not evaluated. The full
Macadamia source is approximately 28.9× under the archived 780-Mb genome
size assumption, so a 30× subset cannot be obtained without upsampling;
the Ey15 enrolled read subset is far below a full 30× source. Future coverage
studies must report nominal whole-library coverage and observed eligible
local-spanner count separately, and retain zero-support cases in their
denominator.

M1 absolute family abundance remains NO-GO. The read-span signed deficit is
a locus-specific development statistic and must not be substituted for the
production family-wide copy-number estimate or treated as independent
physical copy truth.
