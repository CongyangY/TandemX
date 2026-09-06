# TRASH1 completed default factorial assembly evaluation

The unchanged10-Mb independent s6301 genome contains55 founder families and
2,201,400 planted repeat bases, including one1.026-Mb array. Pinned TRASH commit
f7d53a0 ran de novo without templates in the same immutable ARM64 R4.4.3 image
as TRASH2; one CPU and8 GiB allocation. Native GNU time:2,299.44 s,513.2266 MiB
peak RSS. Cgroup memory peak624,246,784 bytes includes cache. This older wrapper
has no complete cgroup CPU counters. It is an assembly test, not comparable to
macOS HiFi-read timings. All original outputs and failure receipts are retained.

The49 native regions recover43/55 families using primary consensus alone and
**48/55** using all primary and secondary native consensus sequences at the
same exact cyclic edit-identity threshold0.9. The latter bank has98 in-scope
sequences and84 distinct canonical consensuses. Both endpoints remain visible;
the native alternative sequence is not silently omitted. Scope is30–1000 bp.

Native period/region matching gives26/55 array matches under IoU0.5 plus
period tolerance max(2 bp,2% truth period), for both coordinate alternatives.
Conditional boundary MAE is350.12 bp under R-extraction coordinates. These
period-sensitive matches must not be reported as pure base-detection recall.
All11,257 monomer rows exactly match the input reference with the independently
audited−1 bp adjustment; none match offsets−2,0,+1,+2. At−1, unit-base precision
is0.9993707471, recall0.9927064595 and F10.9960274559, from2,186,720 union bases.
The unadjusted policy is also reported; original coordinates are unchanged.

The first evaluator failed because one native period is171.5, a valid fractional
estimate unsupported by its original integer parser. `evaluation_v1` records
that failure. `evaluation_v2` preserves fractional periods and completed;
`evaluation_v3` adds the explicit primary/secondary catalogue sensitivity while
retaining earlier region and unit endpoints. Neither failure nor alternative
endpoint is presented as an independent biological replicate. One IID genome
seed does not establish real-plant accuracy or final resource superiority.
