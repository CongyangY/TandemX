# TRASH comparator Linux environment

This is optional comparator infrastructure, not a TandemX dependency. The
original TRASH primarily supports Linux; its macOS branch would use Windows
paths. TRASH2 additionally requires mafft/nhmmer. A dedicated Linux container
avoids changing the host R installation or the `tandemx-dev` environment.

The R4.4.3 multi-architecture base is pinned by image-index digest, and both
author repositories are pinned by commit. The observed native ARM64 image
manifest is `sha256:0792a280d2a8b0cf6e08b7bf5338ade551af71fa9af3f520c820df253aa42295`.
The host's unrelated existing BRAKER image is amd64 and is not used for native
performance comparisons. Compilation uses two jobs; test containers must set
explicit CPU/memory limits and mount only selected inputs/outputs.

```bash
docker build --platform linux/arm64 -t tandemx-comparators:trash-r443-20260906 \
  benchmarks/containers/trash
docker run --rm --network none --cpus 1 --memory 8g \
  tandemx-comparators:trash-r443-20260906 Rscript --version
python -m benchmarks.scripts.run_trash_container \
  --fasta /absolute/path/control.fa --tool trash \
  --image tandemx-comparators:trash-r443-20260906 --outdir /absolute/path/new_run
python -m benchmarks.scripts.audit_trash_coordinates \
  --fasta /absolute/path/control.fa --native /absolute/path/new_run/native \
  --tool trash --outdir /absolute/path/new_coordinate_audit
pytest -q tests/unit/test_trash_container.py tests/unit/test_trash_coordinates.py
```

No TRASH installation script is run. Dependencies are installed at build time,
and their source archives, SHA-256 values, R/package versions and Debian package
versions are retained under `/opt`. CRAN transitive versions are not yet a
frozen rebuild lock; retain the completed image ID and provenance before testing.
The container supplies TRASH1's expected mafft path through an executable symlink;
its R algorithms are unchanged. Use `--def` for TRASH1 to find the container R
library. TRASH2 must run from `/opt/TRASH_2/src` and must not install anything
during measured inference; an offline preflight is required.

Image construction alone is not a completed comparator benchmark. Test native
array/monomer workflows and inspect outputs before larger assembly runs. Bundled
HOR executables may have a different architecture and are outside this initial
array benchmark. Compare time/memory only with tools run in the same environment;
container timings must not be ranked directly against macOS diagnostics.

The runner resolves the image tag to an immutable ID, validates the input FASTA,
uses one CPU/8 GiB and no network, preserves output/logs plus package provenance,
and checks the native exit marker and output headers. GNU time measures the
command inside Linux; its RSS is a waited-process maximum, not a sum of concurrent
descendants. Separate cgroup memory peak includes cache. The host Docker client
measurement is labelled separately. The input hash is rechecked after execution.
Use `trash2` for its distinct interface; TRASH1 uses seed6101, TRASH2's pinned
source sets seed0. Neither receives repeat templates or HOR options here.

The two initial T7 bind-mount creation requests timed out before inference and
were retained as failed. Their delayed, never-started containers were explicitly
removed by their recorded names. Subsequent runs completed on the identical
author-supplied100-kb human example. That input validates installation/adapters,
not plant performance or an additional cohort species. The coordinate audit
is explicitly restricted to1-Mb controls and evaluates five offsets without
altering native outputs. On this pinned TRASH1 source, all355 monomer sequences
match strand-adjusted reference extraction only at offset-1 from reported
1-based coordinates; TRASH2 has360 matches at offset0. Array-coordinate and
larger plant/simulated accuracy evaluation remain separate work.
