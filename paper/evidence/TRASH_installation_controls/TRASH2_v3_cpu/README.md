# Whole-container CPU and repeat-run control

Same100-kb author input and immutable image as the prior successful TRASH2
control. The run completed in30.37 seconds. Whole-cgroup CPU was29.761074 seconds
(28.664724 user,1.096349 system), compared with3.31+0.54 seconds from GNU time.
This validates the need to include R worker processes. Native cgroup counters,
throttling and the exact measurement wrapper are retained.

The two executions are **not byte-identical**. Their seven region geometries
agree, but monomer counts are360 and363 and representatives differ. Unit coverage
Jaccard is0.996439; this measures output agreement, not biological accuracy.
Both results remain preserved. Upstream sets a top-level seed but also samples
units inside parallel workers; that is a plausible source of variability, not
a proven cause. No external source was patched and no result was selected for
more favorable accuracy. The input is a human author example, not a plant cohort
member or new accuracy benchmark.
