# TideCluster comparator container

This image is a task-matched assembly comparator for TandemX. It installs
TideCluster 1.21.2 from the author's conda channel and verifies the upstream
TideHunter 1.4.3 and `kitehor` runtime pins during the build. The `linux/amd64`
base is pinned by platform-manifest digest because TideHunter 1.4.3 is not
available and its bundled abPOA code does not compile on macOS ARM64.

The Bioconda TideHunter binary itself uses instructions unavailable through
Docker's amd64 emulation on Apple Silicon. The build stage uses the pinned
linux/amd64 BRAKER3 image because it already supplies Git, GCC and zlib headers;
only `g++` and `make` are added from its Ubuntu repositories. The image then
verifies and rebuilds upstream TideHunter commit
`d4a7feb2b9edfefb267a30db5d6f0685fb79ea32` with its pinned abPOA submodule and
the project's documented `sse2=1` target, then replaces only that executable.
The conda MMseqs2 16.747c6 binary has the same emulation incompatibility. It is
replaced with the official release's SSE2 build after checking SHA-256
`2f18cc42a31e64cd33207c5069ef9cf5d54bfd50ec6e550b4d65ee534a4ba402`
and exact version commit `747c64cc8db3b4803a0f1194a3f75b3ba9f81bcb`.

Build and record the resulting image digest:

```bash
docker build --platform linux/amd64 \
  -f benchmarks/containers/tidecluster/Dockerfile \
  -t tandemx/tidecluster:1.21.2 .
docker image inspect tandemx/tidecluster:1.21.2 \
  --format '{{index .RepoDigests 0}} {{.Id}}'
```

Run inputs from a mounted directory and count container startup/emulation in
end-to-end resource measurements:

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" -w /work tandemx/tidecluster:1.21.2 \
  tidehunter -c 1 -pr output/tc -f assembly.fa
```

The image also contains GNU `time` in the comparator environment. Use it inside
the container when a benchmark needs child-process maximum RSS; keep that field
separate from host-side Docker startup and process-tree measurements.

Installation/build failures are retained as reproducibility outcomes. They are
never converted to zero accuracy. The read-local TideHunter comparison and the
assembly-level TideCluster comparison remain separate benchmark tasks.
