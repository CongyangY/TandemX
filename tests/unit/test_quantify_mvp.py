from __future__ import annotations

from dataclasses import replace
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
    select_haploid_depth,
    validate_quantify_config,
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


def test_quality_window_survival_fast_fasta_path_preserves_ambiguous_breaks() -> None:
    observed = quality_window_survival("AACGTNACGTTA", None, 3, 0.01)
    expected_windows = 3 + 4
    assert observed == (
        expected_windows * (1.0 - 0.01) ** 3,
        expected_windows,
        0,
        expected_windows,
    )


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


def test_single_copy_depth_gate_selects_a_recorded_total_bases_fallback(tmp_path: Path) -> None:
    config = QuantifyConfig(
        reads=tmp_path / "reads.fa",
        monomers=tmp_path / "monomers.fa",
        genome_size=100,
        outdir=tmp_path / "out",
        k=5,
        haploid_depth=None,
        single_copy_kmers=tmp_path / "controls.tsv",
        single_copy_min_depth=2.0,
    )
    stats = estimate_control_depth([0.0, 1.0])
    depth, method = select_haploid_depth(config, stats, total_read_bases=100)
    assert depth == 1.0
    assert method == "total_read_bases_divided_by_genome_size_low_control_depth_fallback"

    ungated = replace(config, single_copy_min_depth=None)
    assert select_haploid_depth(ungated, stats, 100) == (
        0.5,
        "empirical_single_copy_kmers_mean",
    )


def test_single_copy_depth_gate_requires_controls(tmp_path: Path) -> None:
    config = QuantifyConfig(
        reads=tmp_path / "reads.fa",
        monomers=tmp_path / "monomers.fa",
        genome_size=100,
        outdir=tmp_path / "out",
        k=5,
        haploid_depth=None,
        single_copy_min_depth=2.0,
    )
    with pytest.raises(ValueError, match="requires --single-copy-kmers"):
        validate_quantify_config(config)


def test_single_copy_depth_gate_rejects_explicit_haploid_depth(tmp_path: Path) -> None:
    controls = tmp_path / "controls.tsv"
    controls.write_text("kmer\texpected_copy_number\nACGTT\t1\n")
    config = QuantifyConfig(
        reads=tmp_path / "reads.fa",
        monomers=tmp_path / "monomers.fa",
        genome_size=100,
        outdir=tmp_path / "out",
        k=5,
        haploid_depth=1.0,
        single_copy_kmers=controls,
        single_copy_min_depth=2.0,
    )
    with pytest.raises(ValueError, match="cannot be combined with --haploid-depth"):
        validate_quantify_config(config)

@pytest.mark.parametrize("sequence,k,targets", [
    ("ACGTACGTNNacgt", 3, {"ACG", "CGT", "AAA", "TTT"}),
    ("ACGTACGT", 31, {"A" * 31}),
    ("ACGTACGT", 32, {"A" * 32}),
    ("ACGΩTACG", 3, {"ACG", "CGT"}),
])
def test_selected_python_rolling_codes_match_legacy_strings(sequence, k, targets):
    from collections import Counter
    from tandemx.quantify.mvp import (
        iter_kmers, selected_target_code_map, update_selected_python_kmer_counts,
    )

    expected = Counter(word for word in iter_kmers(sequence, k) if word in targets)
    observed = Counter()
    update_selected_python_kmer_counts(
        observed, sequence, k, targets, selected_target_code_map(targets, k)
    )
    assert observed == expected

def test_full_quantify_rolling_codes_match_legacy_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import tandemx.quantify.mvp as quantify

    reads = tmp_path / "reads.fa"
    catalogue = tmp_path / "catalogue.fa"
    reads.write_text(">r1\nACGTTCAGGACACGTTCAGGACNNACGTTCAGGAC\n>r2\nTCCTGAACGTCTCCTGAACGTC\n", encoding="utf-8")
    catalogue.write_text(">family_id=F1\nACGTTCAGGAC\n>family_id=F2\nTCCTGAACGTC\n", encoding="utf-8")
    base = dict(reads=reads, monomers=catalogue, genome_size=1_000, k=5,
                haploid_depth=1.0, kmer_backend="python")
    quantify.quantify_toy_copy_number(QuantifyConfig(outdir=tmp_path / "new", **base))
    new_bytes = (tmp_path / "new" / "copy_number.tsv").read_bytes()

    def legacy(counts, sequence, k, targets, code_targets):
        counts.update(word for word in quantify.iter_kmers(sequence, k) if word in targets)

    monkeypatch.setattr(quantify, "update_selected_python_kmer_counts", legacy)
    quantify.quantify_toy_copy_number(QuantifyConfig(outdir=tmp_path / "legacy", **base))
    assert (tmp_path / "legacy" / "copy_number.tsv").read_bytes() == new_bytes

@pytest.mark.parametrize("quality", [False, True])
def test_selected_count_quality_and_worker_match_legacy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, quality: bool) -> None:
    import tandemx.quantify.mvp as quantify

    reads = tmp_path / ("reads.fastq" if quality else "reads.fa")
    if quality:
        reads.write_text("@r1\nACGTTCAGGACNACGTTCAGGAC\n+\nIIIIIIIIIIIIIIIIIIIIIII\n", encoding="utf-8")
    else:
        reads.write_text(">r1\nACGTTCAGGACNACGTTCAGGAC\n", encoding="utf-8")
    targets = {quantify.canonical_kmer("ACGTT"), quantify.canonical_kmer("TCAGG")}
    if quality:
        observed = quantify.count_selected_read_kmers_quality_and_bases(
            [reads], 5, targets, "python", max_reads=None, max_read_bases=None,
            progress_every=1, logger=__import__("logging").getLogger(__name__), progress=None,
            assumed_error_rate=None,
        )
    else:
        observed = quantify.count_selected_read_kmers_and_bases_one_file(reads, 5, targets, "python")

    def legacy(counts, sequence, word_k, old_targets, _code_targets):
        counts.update(word for word in quantify.iter_kmers(sequence, word_k) if word in old_targets)

    monkeypatch.setattr(quantify, "update_selected_python_kmer_counts", legacy)
    if quality:
        expected = quantify.count_selected_read_kmers_quality_and_bases(
            [reads], 5, targets, "python", max_reads=None, max_read_bases=None,
            progress_every=1, logger=__import__("logging").getLogger(__name__), progress=None,
            assumed_error_rate=None,
        )
    else:
        expected = quantify.count_selected_read_kmers_and_bases_one_file(reads, 5, targets, "python")
    assert observed == expected
