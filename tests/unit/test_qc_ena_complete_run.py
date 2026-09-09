from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

from benchmarks.scripts.qc_ena_complete_run import qc_run


def test_qc_run_matches_hash_and_aggregate_metadata(tmp_path: Path) -> None:
    fastq_dir = tmp_path / "reads"
    fastq_dir.mkdir()
    paths = [fastq_dir / "ERR1234567_1.fastq.gz", fastq_dir / "ERR1234567_2.fastq.gz"]
    for index, path in enumerate(paths, start=1):
        with gzip.open(path, "wt") as handle:
            handle.write(f"@r/{index}\nACGT\n+\nIIII\n")
    md5s = [hashlib.md5(path.read_bytes()).hexdigest() for path in paths]
    sizes = [path.stat().st_size for path in paths]
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\tread_count\tbase_count\tlibrary_layout\n"
        "ERR1234567\t"
        "ftp.sra.ebi.ac.uk/a/ERR1234567_1.fastq.gz;ftp.sra.ebi.ac.uk/a/ERR1234567_2.fastq.gz\t"
        f"{md5s[0]};{md5s[1]}\t{sizes[0]};{sizes[1]}\t1\t8\tPAIRED\n"
    )
    result = qc_run(
        metadata,
        "ERR1234567",
        fastq_dir,
        tmp_path / "qc.json",
        threads=2,
    )
    assert result["status"] == "complete"
    assert result["FASTQ_record_count"] == 2
    assert result["ENA_read_count_field"] == 1
    assert result["read_count_interpretation"] == (
        "ENA_read_count_equals_spots;two_FASTQ_records_per_spot"
    )
    assert result["total_bases"] == 8
    assert result["all_ENA_file_MD5_matched"] is True
