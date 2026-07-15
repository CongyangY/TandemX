from __future__ import annotations

from pathlib import Path

from tandemx.quantify.mvp import (
    QuantifyConfig,
    canonical_kmer,
    is_low_complexity_kmer,
    monomer_kmer_counts,
    quantify_toy_copy_number,
)


def test_canonical_kmer() -> None:
    assert canonical_kmer("ACGA") == "ACGA"
    assert canonical_kmer("TCGT") == "ACGA"


def test_low_complexity_kmer_filter() -> None:
    assert is_low_complexity_kmer("AAAAAAAAAA")
    assert is_low_complexity_kmer("ATATATATAT")
    assert not is_low_complexity_kmer("ACGTTCAGGA")


def test_monomer_kmer_multiplicity() -> None:
    counts = monomer_kmer_counts("ACGTACGT", 4)
    assert counts[canonical_kmer("ACGT")] == 2


def test_quantify_uses_haploid_depth_and_multiplicity_correction(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    outdir = tmp_path / "quantify"
    outdir.mkdir()
    reads.write_text(">r1\nACGTACGTACGTACGT\n", encoding="utf-8")
    monomers.write_text(">family_id=TXF000001;length_bp=8\nACGTACGT\n", encoding="utf-8")

    estimates = quantify_toy_copy_number(
        QuantifyConfig(
            reads=reads,
            monomers=monomers,
            genome_size=16,
            outdir=outdir,
            k=4,
            haploid_depth=1.0,
        )
    )

    assert len(estimates) == 1
    assert estimates[0].diagnostic_kmer_count > 0
    assert estimates[0].estimated_copy_number > 1
    assert estimates[0].warning == "genome_background_uniqueness_not_verified"
    assert estimates[0].confidence == "medium"
    assert (outdir / "copy_number.tsv").is_file()


def test_quantify_warns_when_haploid_depth_is_estimated(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    outdir = tmp_path / "quantify"
    outdir.mkdir()
    reads.write_text(">r1\nACGTACGTACGTACGT\n", encoding="utf-8")
    monomers.write_text(">family_id=TXF000001;length_bp=8\nACGTACGT\n", encoding="utf-8")

    estimates = quantify_toy_copy_number(
        QuantifyConfig(
            reads=reads,
            monomers=monomers,
            genome_size=16,
            outdir=outdir,
            k=4,
            haploid_depth=None,
        )
    )

    assert "haploid_depth_estimated" in estimates[0].warning
    assert estimates[0].confidence == "medium"


def test_python_and_rust_targeted_counting_match(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    reads.write_text(">r1\nACGTTCAGGACACGTTCAGGAC\n", encoding="utf-8")
    monomers.write_text(
        ">family_id=TXF000001;length_bp=11\nACGTTCAGGAC\n",
        encoding="utf-8",
    )

    results = {}
    for backend in ("python", "rust"):
        outdir = tmp_path / backend
        outdir.mkdir()
        results[backend] = quantify_toy_copy_number(
            QuantifyConfig(
                reads=reads,
                monomers=monomers,
                genome_size=22,
                outdir=outdir,
                k=5,
                haploid_depth=1.0,
                kmer_backend=backend,
            )
        )

    assert results["rust"] == results["python"]


def test_quantify_limits_reads_and_bases_consistently(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    reads.write_text(
        ">r1\nACGTTCAGGAC\n>r2\nACGTTCAGGAC\n>r3\nACGTTCAGGAC\n",
        encoding="utf-8",
    )
    monomers.write_text(
        ">family_id=TXF000001;length_bp=11\nACGTTCAGGAC\n",
        encoding="utf-8",
    )

    estimates = quantify_toy_copy_number(
        QuantifyConfig(
            reads=reads,
            monomers=monomers,
            genome_size=22,
            outdir=tmp_path,
            k=5,
            haploid_depth=1.0,
            max_reads=2,
            max_read_bases=22,
        )
    )

    assert estimates[0].median_kmer_depth == 2.0
