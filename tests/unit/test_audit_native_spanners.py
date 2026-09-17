import gzip

import pytest

from benchmarks.scripts.audit_native_spanners import audit, partition_identity


def test_partition_identity_counts_reverse_complement_and_indel():
    reference = "ACGTACGTACGT"
    query = "ACGTTACGTACGT"
    forward = partition_identity(reference, query, "+", 0, 13, 0,
                                 "4M1I8M", (4, 8))
    assert [(row["matches"], row["columns"]) for row in forward] == [
        (4, 4), (4, 5), (4, 4)]
    reverse = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert partition_identity(reference, reverse, "-", 0, 13, 0,
                              "4M1I8M", (4, 8)) == forward


def test_audit_requires_high_identity_despite_geometric_span(tmp_path):
    context = tmp_path / "context.fa"
    context.write_text(">context\n" + "ACGT" * 5 + "\n")
    reads = tmp_path / "reads.fastq.gz"
    with gzip.open(reads, "wt") as out:
        out.write("@read1\n" + "ACGT" * 5 + "\n+\n" + "I" * 20 + "\n")
        out.write("@read2\n" + "T" * 20 + "\n+\n" + "I" * 20 + "\n")
    paf = tmp_path / "map.paf"
    paf.write_text(
        "read1\t20\t0\t20\t+\tcontext\t20\t0\t20\t20\t20\t60\ttp:A:P\tcg:Z:20M\n"
        "read2\t20\t0\t20\t+\tcontext\t20\t0\t20\t5\t20\t60\ttp:A:P\tcg:Z:20M\n"
    )
    result = audit(paf, context, reads, 6, 14, flank_bp=2, minimum_identity=0.95)
    assert result["geometric_spanner_count"] == 2
    assert result["qualified_read_ids"] == ["read1"]
    assert result["geometric_alignments"][1]["array"]["identity"] < 0.95
    with pytest.raises(ValueError, match="minimum identity"):
        audit(paf, context, reads, 6, 14, flank_bp=2, minimum_identity=1.1)
