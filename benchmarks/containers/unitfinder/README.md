# unitFinder comparator container

This container freezes the official unitFinder checkout at commit
`e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28` and supplies the dependencies
required by its primary single-chromosome command: TRF, MUMmer4, Biopython,
NumPy and Matplotlib. Scikit-learn and Clustal Omega are included for auxiliary
repository and documented cross-chromosome scripts. The conda solver result and
tracked-source hashes are retained in `/opt/provenance`.

The image does not repair or rewrite upstream code. The source audit records one
syntax-invalid auxiliary script, Python 2.7 bytecode files, 85 `os.system`
calls, no dependency manifest and no automated test tree. A successful image
build and `-h` import preflight establish only an executable isolated starting
point.

The primary de novo and supplied-reference modes must be evaluated separately.
The README's chromosome concatenation, manual script modifications and Clustal
Omega steps are not silently folded into a single command; if used, their work
and resource costs must be recorded as separate stages. unitFinder is a plant
centromeric-satellite workflow, not a general read-first repeat detector.

Build without modifying the host Python or conda environments:

```bash
docker build --platform linux/amd64 \
  -t tandemx/unitfinder:e80bff38 benchmarks/containers/unitfinder
```

After the image digest and provenance files pass inspection, run each input in
its own empty working directory so upstream fixed temporary names cannot collide.
Generate the deterministic 308,375-bp installation/interface smoke with:

```bash
python benchmarks/scripts/generate_unitfinder_smoke.py --outdir /new/output/path
```

Its ten 105--150-copy decoys and one 600-copy target exercise unitFinder's
source-defined copy-number outlier gate. The smoke is not an accuracy benchmark;
the target, manifest, input hashes and any failed native output must be retained.
