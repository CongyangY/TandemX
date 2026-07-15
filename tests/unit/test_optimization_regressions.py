from __future__ import annotations

from pathlib import Path
import random

import pytest

from tandemx.compare.mvp import union_interval_length
from tandemx.discover.mvp import (
    FastaRecord,
    find_best_periodic_candidate_with_stats,
    find_best_short_periodic_candidate_with_stats,
    orient_monomer,
)
from tandemx.discover.rust_backend import rust_backend_available, scan_read_for_periods
from tandemx.io.sequences import DuplicateIdTracker
from tandemx.locate.mvp import (
    LocateConfig,
    build_family_kmer_index,
    covered_bp,
    locate_toy_arrays,
    locate_record_arrays,
    window_density,
)
from tandemx.pipeline import PipelineConfig, step_fingerprint
from tandemx.probe.mvp import candidate_probe_sequences, long_oligo_tm
from tandemx.quantify.mvp import MonomerRecord, monomer_kmer_counts
from tandemx.utils.kmers import canonical_kmer_code, iter_canonical_kmer_codes


def test_duplicate_tracker_spills_and_remains_exact() -> None:
    with DuplicateIdTracker(memory_limit=2) as tracker:
        assert tracker.add("r1")
        assert tracker.add("r2")
        assert tracker.add("r3")
        assert tracker.spilled
        assert not tracker.add("r1")
        assert not tracker.add("r3")


def test_rolling_kmer_codes_match_direct_encoding_and_reset_on_n() -> None:
    sequence = "ACGTTCNAGGAC"
    observed = dict(iter_canonical_kmer_codes(sequence, 5))
    expected = {
        index: canonical_kmer_code(sequence[index : index + 5])
        for index in range(len(sequence) - 4)
        if "N" not in sequence[index : index + 5]
    }
    assert observed == expected


def test_circular_monomer_kmers_are_phase_invariant_even_when_k_is_longer() -> None:
    monomer = "ACGTTCAG"
    shifted = monomer[3:] + monomer[:3]
    assert monomer_kmer_counts(monomer, 13) == monomer_kmer_counts(shifted, 13)
    assert orient_monomer(monomer) == orient_monomer(shifted)


def test_overlapping_arrays_are_not_double_counted_in_density_or_abundance() -> None:
    intervals = [(10, 40), (20, 60), (50, 80)]
    assert covered_bp(intervals, 0, 100) == 70
    assert union_interval_length(intervals) == 70
    record = FastaRecord("chr1", "chr1", "A" * 100)
    density = window_density(record, intervals, window_size=100, step_size=100)
    assert density[0].score == pytest.approx(0.7)


def test_locate_min_identity_changes_array_acceptance() -> None:
    monomer = MonomerRecord("TXF000001", "ACGTTCAGGACTAACCGTGA")
    sequence = list(monomer.sequence * 6)
    sequence[len(sequence) // 2] = "A" if sequence[len(sequence) // 2] != "A" else "C"
    record = FastaRecord("chr1", "chr1", "".join(sequence))
    index, indexed, shared = build_family_kmer_index([monomer], 11)
    relaxed = locate_record_arrays(record, [monomer], index, indexed, shared, 11, 0.6)
    strict = locate_record_arrays(record, [monomer], index, indexed, shared, 11, 0.95)
    assert relaxed
    assert not strict


def test_chunked_locate_preserves_kmers_across_chunk_boundaries(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import tandemx.io.sequences as sequence_io
    import tandemx.locate.mvp as locate_module

    monomer = "ACGTTCAGGACTAACCGTGA"
    assembly = tmp_path / "assembly.fa"
    catalog = tmp_path / "monomers.fa"
    assembly.write_text(">chr1\n" + "GCTA" * 8 + monomer * 6 + "TGCA" * 8 + "\n", encoding="utf-8")
    catalog.write_text(">family_id=TXF000001;length_bp=20\n" + monomer + "\n", encoding="utf-8")

    default = locate_toy_arrays(
        LocateConfig(assembly, catalog, None, tmp_path / "default", 50, 25, 11, 0.8)
    )
    monkeypatch.setattr(
        locate_module,
        "read_fasta_chunks",
        lambda path: sequence_io.read_fasta_chunks(path, chunk_bases=37),
    )
    chunked = locate_toy_arrays(
        LocateConfig(assembly, catalog, None, tmp_path / "chunked", 50, 25, 11, 0.8)
    )

    assert chunked == default


def test_short_monomers_generate_tandem_context_probes_and_tm_responds_to_formamide() -> None:
    probes = candidate_probe_sequences("ACGTTCAG", min_len=80, max_len=120)
    assert probes
    assert all(len(probe) == 80 for probe in probes)
    assert long_oligo_tm(probes[0], formamide_percent=20) < long_oligo_tm(probes[0])


def test_pipeline_fingerprint_changes_when_input_content_changes(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">r1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        reads=(reads,),
        assembly=None,
        genome_size=100,
        haploid_depth=None,
        outdir=tmp_path / "run",
        max_reads=1,
        max_read_bases=None,
        kmer_backend="python",
        steps=("discover",),
        min_period=2,
        max_period=10,
        top_periods=2,
        threads=1,
        resume=False,
        force=False,
        profile=False,
    )
    before = step_fingerprint(config, "discover")
    reads.write_text(">r1\nTGCA\n", encoding="utf-8")
    after = step_fingerprint(config, "discover")
    assert before != after


@pytest.mark.skipif(not rust_backend_available(), reason="TandemX Rust extension is not installed")
def test_python_and_rust_local_repeat_boundaries_are_identical() -> None:
    monomer = "ACGTTCAGGACTAACCGTGATCGATCGATCG"
    rng = random.Random(17)
    left = "".join(rng.choices("ACGT", k=217))
    right = "".join(rng.choices("ACGT", k=203))
    sequence = left + monomer * 12 + right
    record = FastaRecord("read1", "read1", sequence)
    python_candidate, _ = find_best_periodic_candidate_with_stats(
        record,
        min_period=20,
        max_period=50,
        min_repeat_span=200,
        candidate_index=0,
        kmer_size=11,
        top_periods=5,
        min_seed_occurrences=2,
        min_spacing_support=2,
        max_pairs_per_kmer=100,
    )
    rust = scan_read_for_periods(
        sequence,
        k=11,
        min_period=20,
        max_period=50,
        top_periods=5,
        min_seed_occurrences=2,
        min_spacing_support=2,
        max_pairs_per_kmer=100,
        min_repeat_span=200,
    )
    assert python_candidate is not None
    assert rust.best_period == python_candidate.period_bp
    assert rust.repeat_start == python_candidate.read_start
    assert rust.repeat_end == python_candidate.read_end
    assert rust.periodicity_score == pytest.approx(python_candidate.score)


@pytest.mark.skipif(not rust_backend_available(), reason="TandemX Rust extension is not installed")
def test_two_base_nontrivial_short_repeat_has_backend_parity() -> None:
    monomer = "AACACCAACCCA"
    rng = random.Random(19)
    sequence = "".join(rng.choices("ACGT", k=60)) + monomer * 20 + "".join(rng.choices("ACGT", k=60))
    record = FastaRecord("read1", "read1", sequence)
    python_candidate, _ = find_best_short_periodic_candidate_with_stats(
        record,
        min_period=12,
        max_period=12,
        min_repeat_span=120,
        candidate_index=0,
    )
    rust = scan_read_for_periods(
        sequence,
        k=11,
        min_period=12,
        max_period=12,
        top_periods=3,
        min_seed_occurrences=2,
        min_spacing_support=2,
        max_pairs_per_kmer=100,
        min_repeat_span=120,
    )
    assert python_candidate is not None
    assert rust.status == "accepted"
    assert (rust.best_period, rust.repeat_start, rust.repeat_end) == (
        python_candidate.period_bp,
        python_candidate.read_start,
        python_candidate.read_end,
    )
