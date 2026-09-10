import gzip

from benchmarks.scripts.extract_fasta_regions import extract_regions


def test_extract_regions_streams_wrapped_gzip_fasta(tmp_path):
    fasta = tmp_path / "assembly.fna.gz"
    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(">chr1 description\nACGTAC\nGTACGT\n>chr2\nTTTT\n")
    regions = [
        {"region_id": "middle", "sequence_id": "chr1", "start0": 4, "end0": 10},
        {"region_id": "second", "sequence_id": "chr2", "start0": 1, "end0": 4},
    ]
    assert extract_regions(fasta, regions) == {"middle": "ACGTAC", "second": "TTT"}
