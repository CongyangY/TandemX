# Divergence-aware localizer development v1

This compact archive retains the failed first development attempt on seeds
5301--5303. It crossed 1%, 3% and 5% independent founder-to-unit substitutions,
one or three same-family array segments separated by 500 bp, and five assembly
retention levels. The supplied founder catalogue isolates localization; no reads,
copy-number estimates or assembly/read classifications were regenerated.

All 90 localization commands completed. Full-assembly mean base recall was
0.903608 and failed the predeclared minimum of 0.95. Positive-assembly precision
was 0.999435, and none of 54 absent-family rows had predicted repeat bases.
At 5% divergence, exact anchors split long arrays into as many as 119 predicted
fragments. This failure motivated bounded anchor bridging and remains part of
the evidence rather than being replaced by the selected development result.

`archive_manifest.json` verifies 24 compact files, including all localization
rows, 90 receipt hashes/resource measurements and the frozen source snapshot.
The simulation uses independent length-preserving substitutions and a known
catalogue. The IID proxy is not alignment identity or biological validation.
