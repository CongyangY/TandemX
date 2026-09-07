from __future__ import annotations

from pathlib import Path

import pytest

from tandemx.quantify.mvp import (
    QuantifyConfig,
    canonical_kmer,
    is_low_complexity_kmer,
    estimate_control_depth,
    monomer_kmer_counts,
    quality_window_survival,
    quantify_toy_copy_number,
    read_single_copy_kmers,
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


def test_quality_window_survival_uses_phred_probabilities_and_skips_n() -> None:
    survival, windows, quality_windows, assumed_windows = quality_window_survival(
        "ACGTNACGT", "5555!5555", 3, None
    )
    expected = 4 * (1.0 - 0.01) ** 3
    assert abs(survival - expected) < 1e-12
    assert (windows, quality_windows, assumed_windows) == (4, 4, 0)


def test_quantify_corrects_fastq_kmer_survival(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    monomers = tmp_path / "monomers.fa"
    sequence = "ACGTTCAGGACACGTTCAGGAC"
    reads.write_text(f"@r1\n{sequence}\n+\n{'I' * len(sequence)}\n")
    monomers.write_text(
        ">family_id=TXF000001;length_bp=11\nACGTTCAGGAC\n"
    )
    estimates = quantify_toy_copy_number(QuantifyConfig(
        reads=reads, monomers=monomers, genome_size=len(sequence),
        outdir=tmp_path / "out", k=5, haploid_depth=1.0,
    ))
    estimate = estimates[0]
    assert abs(estimate.kmer_survival_probability - (1.0 - 0.0001) ** 5) < 1e-12
    assert estimate.median_kmer_depth > estimate.raw_median_kmer_depth
    assert estimate.quality_window_count == len(sequence) - 5 + 1
    assert "phred_scores_treated_as_calibrated_probabilities" in estimate.warning


def test_empirical_single_copy_kmers_set_depth_when_explicit_depth_is_absent(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    controls = tmp_path / "controls.tsv"
    reads.write_text(">r1\nACGTTCAGGACACGTTCAGGACTTTGC\n")
    monomers.write_text(
        ">family_id=TXF000001;length_bp=11\nACGTTCAGGAC\n"
    )
    controls.write_text("kmer\texpected_copy_number\nTTTGC\t1\n")
    estimates = quantify_toy_copy_number(QuantifyConfig(
        reads=reads, monomers=monomers, genome_size=27,
        outdir=tmp_path / "out", k=5, haploid_depth=None,
        single_copy_kmers=controls,
    ))
    estimate = estimates[0]
    assert estimate.normalization_method == "empirical_single_copy_kmers_mean"
    assert estimate.single_copy_control_kmer_count == 1
    assert estimate.single_copy_control_mean_depth == 1.0
    assert estimate.single_copy_control_median_depth == 1.0
    assert estimate.haploid_depth == 1.0
    assert estimate.estimated_copy_number >= 2.0


def test_single_copy_controls_reject_reverse_complement_duplicates(tmp_path: Path) -> None:
    controls = tmp_path / "controls.tsv"
    controls.write_text(
        "kmer\texpected_copy_number\nAACGT\t1\nACGTT\t1\n"
    )
    try:
        read_single_copy_kmers(controls, 5)
    except ValueError as error:
        assert "Duplicate canonical control" in str(error)
    else:  # pragma: no cover - explicit failure message is clearer than assert False
        raise AssertionError("reverse-complement controls were not deduplicated")


def test_single_copy_controls_reject_low_complexity_and_missing_values(tmp_path: Path) -> None:
    controls = tmp_path / "controls.tsv"
    controls.write_text("kmer\texpected_copy_number\nAAAAA\t1\n")
    try:
        read_single_copy_kmers(controls, 5)
    except ValueError as error:
        assert "Low-complexity" in str(error)
    else:  # pragma: no cover
        raise AssertionError("low-complexity control was accepted")
    controls.write_text("kmer\texpected_copy_number\nAACGT\n")
    try:
        read_single_copy_kmers(controls, 5)
    except ValueError as error:
        assert "Invalid expected copy number" in str(error)
    else:  # pragma: no cover
        raise AssertionError("missing control copy number was accepted")


def test_quality_window_survival_rejects_length_and_rate_errors() -> None:
    for quality, rate, expected in (("III", None, "lengths"), (None, 1.0, "rate")):
        try:
            quality_window_survival("ACGT", quality, 3, rate)
        except ValueError as error:
            assert expected in str(error)
        else:  # pragma: no cover
            raise AssertionError("invalid quality-survival input was accepted")


def test_quantify_rejects_disabled_supplied_error_correction(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    monomers = tmp_path / "monomers.fa"
    reads.write_text(">r1\nACGTTCAGGAC\n")
    monomers.write_text(">family_id=TXF000001\nACGTTCAGGAC\n")
    with pytest.raises(ValueError, match="cannot be combined"):
        quantify_toy_copy_number(QuantifyConfig(
            reads=reads,
            monomers=monomers,
            genome_size=11,
            outdir=tmp_path / "out",
            k=5,
            haploid_depth=None,
            read_error_rate=0.01,
            quality_correction_enabled=False,
        ))


def test_control_depth_mean_retains_zero_observations_at_low_depth() -> None:
    stats = estimate_control_depth([0.0, 0.0, 1.0, 2.0])
    assert stats.mean == 0.75
    assert stats.median == 0.5
    assert stats.mad == 0.5
    assert stats.zero_fraction == 0.5
