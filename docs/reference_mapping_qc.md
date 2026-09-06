# Reference-concordance QC on whole-library samples

This is an additional biological-context check after complete FASTQ validation
and fixed-seed sampling. It measures concordance to a declared reference,
including organellar contigs when explicitly supplied. It does not label
unmapped reads contaminants or establish true satellite copy number.

The runner uses the locally recorded minimap2 executable/version/hash and the
author's `map-hifi` preset for HiFi reads. `-c --cs=short` produces base-level
alignments, avoiding approximate PAF match/block fields. Secondary output is
retained with the explicit cap `-N5`; it is not an exhaustive list of repeat
placements. `-t1 -K50M -I3G` sets threads, query batch target and index batch
size. References exceeding3 Gb are rejected in this development QC runner,
so MAPQ is never silently interpreted across a split index. The existing
reference, not an unverified replacement, must match its exact SHA-256.
See the [author manual](https://lh3.github.io/minimap2/minimap2.html).

```bash
python -m benchmarks.scripts.run_reference_mapping \
  --reference /path/to/Col-CEN_v1.2.fasta.gz \
  --reference-sha256 b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8 \
  --sampling-receipt /path/to/ERR6210723_seed6101_v1/sampling_receipt.json \
  --sample-id sample_001 --minimap2 /path/to/minimap2 \
  --organelle-contig ChrM --organelle-contig ChrC \
  --outdir /path/to/new_mapping_qc
pytest -q tests/unit/test_reference_mapping_qc.py
```

All records are counted before mapping. A16-MiB SQLite page cache, disabled
mmap and disk temporary sorting keep identifiers and alignments off the Python
heap. One FASTQ record and one PAF row are processed at a time. Full reference
QC streams sequence lines. Minimap2's own reference index remains in memory;
the controller's bounded memory does not imply constant native index memory.
Source/helper snapshots, exact input hashes and native resources are retained.

PAF names and both sequence lengths must agree with the checked inputs.
Coordinates, strand, counts, unique typed tags and full CIGAR consumption are
validated, including late rows. Native `P/I` records are treated as primary
alignments and `S/i` as secondary. Per-read spans are unioned separately for all
alignments, primary alignments and primary MAPQ20–254 alignments; MAPQ255 is
missing, not high confidence. Unmapped reads remain in every denominator.
Primary query spans aligned to declared organelles and other reference contigs
are also unioned separately. Their overlap is reported explicitly, so a read
with split/overlapping placements cannot silently inflate a compartment total.
These are reference-concordant spans, not proven cellular origins.

Reference coverage uses primary CIGAR `M`, `=` and `X` target blocks. Deletions
and skipped target spans are excluded. Union covered bases and summed aligned
target bases are both reported. The sum divided by reference length is
**alignment-base depth**; overlapping primary/split alignments from one read
can still contribute twice to the sum. It is not independent unique-copy depth.
Secondary placements do not inflate primary reference-depth summaries.

Composition tables group reads into5-percentage-point GC bins and symmetric
max(G+A,C+T)/length bins. The latter is strand invariant but does not measure GA
tract length. N remains in the length denominator. Upper bounds are exclusive
except100%; counts and mapped spans permit explicit denominator checks. These
are composition/concordance associations. Establishing empirical sequencing
bias requires comparison with genomic composition, depth and independent data.

An empty PAF after a successful native exit is zero observed mapping; timeout,
native failure, malformed rows, altered input hashes and incomplete query
counts are failed QC, never zero accuracy. Low MAPQ is retained because repeat
reads need not have unique placements. Reference collapse, known Col-CEN NOR
limitations, genotype differences and material pooling remain unresolved by
this check. Non-organelle reference contigs are labelled `other_reference`,
not assumed nuclear from their names. See the output definitions in
[file_formats.md](file_formats.md#reference-mapping-qc).
