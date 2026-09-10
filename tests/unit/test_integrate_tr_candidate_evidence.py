from benchmarks.scripts.integrate_tr_candidate_evidence import (
    best_short_period,
    gene_context,
    interval_overlap,
    max_homopolymer,
    merge_intervals,
    read_fasta,
    shannon_entropy,
)


def test_sequence_complexity_metrics_distinguish_uniform_and_periodic_sequences():
    assert shannon_entropy("ACGT" * 10) == 2.0
    assert shannon_entropy("A" * 40) == 0.0
    assert max_homopolymer("ACCCCGTT") == 4
    period, fraction = best_short_period("ACGT" * 10, maximum=10)
    assert period == 4
    assert fraction == 1.0


def test_merge_intervals_joins_only_same_sequence_within_gap():
    rows = [
        {"seqid": "chr1", "start0": 100, "end0": 150, "score": 800},
        {"seqid": "chr1", "start0": 200, "end0": 250, "score": 900},
        {"seqid": "chr1", "start0": 400, "end0": 430, "score": 700},
        {"seqid": "chr2", "start0": 120, "end0": 170, "score": 600},
    ]
    merged = merge_intervals(rows, maximum_gap=50)
    assert merged == [
        {"seqid": "chr1", "start0": 100, "end0": 250, "score": 900.0, "part_count": 2},
        {"seqid": "chr1", "start0": 400, "end0": 430, "score": 700.0, "part_count": 1},
        {"seqid": "chr2", "start0": 120, "end0": 170, "score": 600.0, "part_count": 1},
    ]


def test_interval_and_gene_context_use_zero_based_half_open_coordinates():
    assert interval_overlap(10, 20, 19, 30)
    assert not interval_overlap(10, 20, 20, 30)
    genes = {
        "chr1": [
            {"start0": 0, "end0": 10, "gene_id": "left"},
            {"start0": 20, "end0": 30, "gene_id": "overlap"},
            {"start0": 40, "end0": 50, "gene_id": "right"},
        ]
    }
    assert gene_context("chr1", 15, 25, genes) == (1, "left", 5, "right", 15)


def test_read_fasta_uses_structured_family_identifier(tmp_path):
    fasta = tmp_path / "monomers.fa"
    fasta.write_text(
        ">family_id=TXF000001;monomer_id=TXM000001;length_bp=8;confidence=high\nACGTACGT\n",
        encoding="utf-8",
    )
    assert read_fasta(fasta) == {"TXF000001": "ACGTACGT"}
