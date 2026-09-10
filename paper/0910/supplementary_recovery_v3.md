# Supplementary Recovery Test

This supplement preserves the complete bounded recovery result and methods moved
out of the main manuscript. It is a negative, bounded test; it does not establish
a recovery method or physical missing-sequence truth.

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
