from __future__ import annotations

from pathlib import Path

from tandemx.discover.mvp import (
    DiscoverConfig,
    FastaRecord,
    discover_toy_repeats,
    find_best_periodic_candidate,
    periodicity_score,
    read_fasta,
)


def test_periodicity_score_detects_exact_period() -> None:
    sequence = "ACGT" * 20
    assert periodicity_score(sequence, 4) == 1.0
    assert periodicity_score(sequence, 5) < 0.5


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
    families_text = (outdir / "families.tsv").read_text(encoding="utf-8")
    assert families_text.startswith("family_id\tmonomer_id\tmonomer_length_bp")
