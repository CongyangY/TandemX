# Complete input acquisition and file QC

Run inside `tandemx-dev`; store full data/SQLite indexes at the chosen data root.
No filtering is applied. Biological QC: [cohort_and_qc.md](cohort_and_qc.md).

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
