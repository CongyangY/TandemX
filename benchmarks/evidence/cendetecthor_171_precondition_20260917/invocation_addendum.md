# Invocation correction before scoring

The first default replay stopped during Snakefile parsing, before any rule,
because `expMonSize=no` on the Snakemake command line is a truthy string; the
source calls `int(config["expMonSize"])` and raises `ValueError`. The upstream
YAML contains `expMonSize: no`, parsed as false. We retain the first manifest,
log and receipt, and remove only the redundant command-line override for both
predeclared conditions. No input sequence, truth, native threshold or scientific
setting changes. The corrected runs use new workdirs and new manifests.
