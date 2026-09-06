# Complete input acquisition and file QC

Run inside `tandemx-dev`; store full data/SQLite indexes at the chosen data root.
No filtering is applied. Biological QC: [cohort_and_qc.md](cohort_and_qc.md).

For other versioned GenBank assemblies, `fetch_ncbi_reference.py` resolves the
exact GCA version from the official NCBI directory listing, requires the expected
BioProject, archives source metadata/checksums and applies an explicit compressed
size budget before transfer. It verifies the complete FASTA totals and retains
missing optional FCS metadata. A resume must match the existing plan and archived
metadata hashes. See the [Nipponbare example and reference limits](nipponbare_reference.md).

```bash
python -m benchmarks.scripts.fetch_ena_complete \
  --metadata ena_metadata.tsv --accession ERR6210723 \
  --outdir /path/to/data/ERR6210723_complete --max-download-gb 12
python -m benchmarks.scripts.qc_complete_fastq \
  --fastq /path/to/data/ERR6210723_complete/ERR6210723.fastq.gz \
  --outdir /path/to/data/ERR6210723_complete/qc \
  --expected-reads 933904 --expected-bases 14646601458
pytest -q tests/unit/test_complete_data_qc.py
```

The archived ENA filereport must contain `run_accession`, `fastq_ftp`,
`fastq_md5`, `fastq_bytes`, with exactly one row for the requested run.
Semicolon-separated file lists must agree in cardinality. Only official ENA
HTTPS FASTQ paths are accepted. The decimal GB budget applies to the complete
run, not a truncated prefix. Up to five attempts use a 60-second socket timeout
and validated byte-range continuation; rerun the same command to resume.
Ignored resume ranges and mismatched existing plans fail explicitly.

`download_plan.json` retains the source row, size/MD5, URL and source hashes.
Full size/MD5 must match before `.partial` is renamed. Existing files are
reverified. `download_receipt.json` adds SHA-256 and per-file transfer state;
its completion refers to transfer/checksum only. `transfer.log` retains errors.

QC parses all four-line FASTQ records, validates sequence/quality length,
ACGTN/Phred+33, optional plus IDs, exact archive ID uniqueness and expected
totals, and consumes all gzip members/trailers. Lines at the five-million-byte
limit fail explicitly. Histograms retain distributions; a disk-backed SQLite
B-tree checks IDs with a 16 MiB cache. It never retains all reads or IDs in RAM.
Use a new QC directory; failures leave `complete=false` and an error receipt.

`qc.json` includes input/script hashes, observed read/base totals, length
min/max/median/N50, GC/N fractions, and mean reported error probability with
its Phred transform. Read quality uses error probabilities, not arithmetic
mean Phred. These scores are not empirically measured accuracy or automatic CN
correction factors. `read_ids.sqlite` is retained; archive-renamed IDs cannot
establish independence of source molecules. Reference FASTA QC is separate;
nuclear/organellar denominators and known assembly issues remain explicit.

## Reproducible nested whole-file sampling

```bash
python -m benchmarks.scripts.sample_complete_fastq \
  --fastq /path/to/data/ERR6210723_complete/ERR6210723.fastq.gz \
  --qc-receipt /path/to/data/ERR6210723_complete/qc/qc.json \
  --outdir /path/to/new/samples --namespace ERR6210723 --seed 6101 \
  --fractions 0.001 0.01 0.1 --genome-size 131559676
pytest -q tests/unit/test_complete_sampling.py
```

Sampling requires successful full-file QC and scans every record through EOF,
recomputing compressed-file SHA-256 and observed totals in the same pass. A
mismatch fails and leaves `.partial` files. The shared bounded parser also
hashes QC input during parsing, avoiding a second raw-file scan. Caller and
parser hashes are recorded.

A seeded, library-namespace-specific BLAKE2b 128-bit primary-read-ID digest is
compared with integer fraction thresholds. Reads have equal inclusion probability
regardless of sequence/length/quality/repeat status. QC has already checked exact
ID uniqueness. Input order cannot change membership. Higher fractions include
lower ones for the same seed/namespace; output order follows the source.
Fractions are Bernoulli probabilities, not exact base/read counts. Whole reads
are retained. Up to 12 fractions share a scan. Use a stable accession namespace.

Deterministic gzip output preserves reads/qualities with LF line endings. Each
subset has an ID/hash/length TSV, length and joint length/GC/quality histograms.
`sampling_plan.json` records source/QC/scripts and assumptions; the receipt
records expected and observed read/base totals, hashes, N50, GC/N, reported error
probability and nominal bases/genome size. This ratio is not measured nuclear
depth. Empty subsets have an explicit status, zero counts and missing distribution
statistics; they cannot be presented as tested positive data. Seeds assess
sampling stability, not extra plants. Missing libraries and library-selection
bias remain unresolved by uniform sampling of included files.

The versioned Mo17 reference has a dedicated reproducible acquisition command:

```bash
python -m benchmarks.scripts.fetch_mo17_reference \
  --outdir /path/to/data/references/Mo17_T2T_GCA_022117705.1
python -m benchmarks.scripts.plot_complete_qc \
  --qcdir /path/to/full/qc --samplesdir /path/to/samples \
  --outdir /path/to/new/figure --label 'Mo17 CCS: SRR15447419'
pytest -q tests/unit/test_reference_qc.py
```

The reference command pins GCA_022117705.1, verifies its PRJNA751841 project,
archives NCBI metadata/checksums, limits compressed FASTA to 1 GB, resumes an
existing download plan and verifies source MD5 plus all contigs/base totals.
FASTA is consumed line by line with bounded lines; complete chromosomes are not
kept in RAM. Duplicate/empty IDs, malformed sequence and truncation fail. IUPAC
ambiguities are counted explicitly. Reference identity and zero Ns alone do not
prove every satellite copy correct. Mo17 library and assembly share cultivar and
study context; exact donor/extraction identity remains unresolved.

MorexV3 can be acquired directly from its original IPK publication:

```bash
python -m benchmarks.scripts.fetch_morex_reference \
  --outdir /path/to/data/references/MorexV3_IPK_2021_3
pytest -q tests/unit/test_morex_reference.py
```

The [original FASTA record](https://doi.ipk-gatersleben.de/DOI/b2f47dfb-47ff-4114-89ae-bad8dcc515a1/b6e6a2e5-2746-4522-8465-019c8f56df7f/1)
publishes SHA-256 `54c98a04d13ff97350f5f3a5bfa45ac395ad640df8bb1f7598eca4e7edb437c1`.
The script archives the source metadata and pins that file/checksum under DOI
10.5447/ipk/2021/3 (CC BY 4.0). A 4.5-GB transfer cap bounds the original,
uncompressed FASTA. Chunked transfers are supported; resume requires the exact
HTTP range, and access-denied responses are not retried. No failed partial is
renamed or silently discarded. Full streaming FASTA QC rechecks the pinned SHA,
contig/ambiguity counts and structural validity. `reference_receipt.json` marks
completion only after these checks. Matching cultivar/study does not establish
identical donor DNA or precise satellite copy truth.
