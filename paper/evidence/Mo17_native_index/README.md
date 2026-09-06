# Native-index ablation on actual Mo17 candidates

Committed source `4ae8ca8`; same frozen core/native alignment binary in both
children, same85,663 serialized candidates from1.129 Gb and28,586 families.
All family/member values match exactly (hash in each worker receipt).

| Gate | Clustering seconds | Child seconds | Peak RSS MiB |
| --- | ---: | ---: | ---: |
| Python packed postings | 218.0149 | 220.0886 | 250.2188 |
| Native compact postings | 164.3125 | 166.1930 | 259.6250 |

Clustering is24.63% shorter, while child peak RSS is3.76% higher. This is a
single fixed-order diagnostic with concurrent jobs, not final isolated timing,
external superiority or full live-candidate pipeline parity. Both trade-offs
are retained. Candidate export rounding remains part of this replay's scope.
