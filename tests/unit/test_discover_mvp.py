from __future__ import annotations

from pathlib import Path

from tandemx.discover.mvp import (
    DiscoverConfig,
    FastaRecord,
    RepeatFamily,
    compare_family_pair,
    discover_toy_repeats,
    find_best_periodic_candidate,
    periodicity_score,
    read_fasta,
)


def test_periodicity_score_detects_exact_period() -> None:
    sequence = "ACGT" * 20
    assert periodicity_score(sequence, 4) == 1.0
    assert periodicity_score(sequence, 5) < 0.5


def test_discover_core_does_not_reference_simulator_truth_files() -> None:
    source = Path("tandemx/discover/mvp.py").read_text(encoding="utf-8")

    assert "truth_monomers.fa" not in source
    assert "truth_arrays.bed" not in source
    assert "truth_copy_number.tsv" not in source


def test_find_best_periodic_candidate() -> None:
    record = FastaRecord(read_id="read1", description="read1;strand=-", sequence="AACCGGTT" * 8)
    candidate = find_best_periodic_candidate(
        record,
        min_period=4,
        max_period=12,
        min_repeat_span=16,
        candidate_index=1,
    )
    assert candidate is not None
    assert candidate.period_bp == 8
    assert candidate.strand == "-"
    assert candidate.confidence == "high"


def test_read_fasta_accepts_fastq_input(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")

    records = list(read_fasta(reads))

    assert len(records) == 1
    assert records[0].read_id == "r1"
    assert records[0].sequence == "ACGT"


def test_discover_toy_repeats_writes_documented_outputs(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">r1;strand=+\nACGTACGTACGTACGT\n>r2;strand=-\nACGTACGTACGTACGT\n", encoding="utf-8")
    outdir = tmp_path / "discover"
    outdir.mkdir()

    candidates, families = discover_toy_repeats(
        DiscoverConfig(
            reads=reads,
            outdir=outdir,
            min_monomer_len=4,
            max_monomer_len=8,
            min_support_reads=1,
            min_repeat_span=8,
        )
    )

    assert len(candidates) == 2
    assert len(families) == 1
    assert (outdir / "candidate_reads.tsv").is_file()
    assert (outdir / "monomers.fa").is_file()
    assert (outdir / "families.tsv").is_file()
    assert (outdir / "family_similarity.tsv").is_file()
    families_text = (outdir / "families.tsv").read_text(encoding="utf-8")
    assert families_text.startswith("family_id\tmonomer_id\tmonomer_length_bp")


def test_family_similarity_flags_possible_related_monomers() -> None:
    family_a = RepeatFamily(
        family_id="TXF000001",
        monomer_id="TXM000001",
        monomer_sequence=("ACGTTCAGGACTAACCGTGA" * 8)[:120],
        monomer_length_bp=120,
        support_read_count=10,
        support_span_bp=1200,
        mean_identity=0.95,
        low_complexity_flag=False,
        confidence="high",
        warning="",
    )
    family_b = RepeatFamily(
        family_id="TXF000002",
        monomer_id="TXM000002",
        monomer_sequence=family_a.monomer_sequence + family_a.monomer_sequence[:80],
        monomer_length_bp=200,
        support_read_count=2,
        support_span_bp=400,
        mean_identity=0.9,
        low_complexity_flag=False,
        confidence="medium",
        warning="",
    )

    similarity = compare_family_pair(family_a, family_b, k=11)

    assert similarity.relationship in {"likely_redundant", "possible_higher_order_or_partial"}
    assert similarity.local_identity >= 0.9
    assert similarity.shared_kmer_fraction >= 0.8
