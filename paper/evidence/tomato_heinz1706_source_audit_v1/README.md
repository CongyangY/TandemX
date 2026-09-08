# Tomato Heinz 1706 donor-match source audit v1

This is a negative donor-eligibility result, not a TandemX accuracy result.
The same named cultivar has an attractive old/new assembly pair: the HiFi+Hi-C
SL5.0 assembly is public, and the later ONT-based SL-T2T assembly adds sequence
concentrated in difficult repetitive regions. The official tomato database also
exposes the 831,451,202-byte `SL-T2T.fa` file directly.

The source paper nevertheless does not support an identical donor. It states
that small SL5.0/SL-T2T/Verkko differences can reflect heterozygosity from mixed
seeds or sample differences between the ONT and HiFi sequencing. Cultivar-name
agreement is therefore insufficient for the preregistered Tier A rule. The
candidate is retained as `not_source_eligible_donor_mismatch_risk`, not promoted
to a favourable benchmark.

The old assembly download is 239,182,278 bytes compressed with published MD5
`4ea3fc4a3ef92083bed44273ed194b5b`; the exact Heinz 1706 HiFi run
`SRR15243707` contains 34,431,605,563 reported bases in 2,072,446 reads and is
28,768,190,557 bytes compressed. None of these large files were downloaded.
The audit queried only metadata, checksums and HTTP response headers.

Current fate: `not_source_eligible_donor_mismatch_risk`. The pair may be used
later as same-cultivar descriptive context under a separately frozen Tier B
design, but it must not be scored as donor-matched assembly-collapse truth.

