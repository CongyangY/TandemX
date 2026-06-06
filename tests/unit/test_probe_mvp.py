from __future__ import annotations

from pathlib import Path

from tandemx.probe.mvp import (
    ArrayRegion,
    candidate_probe_sequences,
    gc_content,
    low_complexity_ratio,
    overlaps_family_array,
    rank_toy_probes,
    simple_tm,
    ProbeConfig,
)


def test_candidate_probe_sequences() -> None:
    probes = candidate_probe_sequences("ACGT" * 100, min_len=80, max_len=120)
    assert probes
    assert all(len(probe) == 120 for probe in probes)


def test_low_complexity_and_tm() -> None:
    assert low_complexity_ratio("A" * 100) == 1.0
    assert gc_content("ACGT") == 0.5
    assert simple_tm("ACGT") == 12.0


def test_overlaps_family_array() -> None:
    arrays = [ArrayRegion("chr1", 10, 50, "TXF000001")]
    assert overlaps_family_array(("chr1", 20, 30), "TXF000001", arrays)
    assert not overlaps_family_array(("chr1", 60, 70), "TXF000001", arrays)
    assert not overlaps_family_array(("chr1", 20, 30), "TXF000002", arrays)


def test_rank_toy_probes_filters_low_complexity(tmp_path: Path) -> None:
    monomers = tmp_path / "monomers.fa"
    assembly = tmp_path / "assembly.fa"
    copy_number = tmp_path / "copy_number.tsv"
    arrays = tmp_path / "arrays.bed"
    outdir = tmp_path / "probe"
    outdir.mkdir()

    monomers.write_text(
        ">family_id=TXF000001;length_bp=120\n" + "A" * 120 + "\n"
        ">family_id=TXF000002;length_bp=120\n" + "ACGT" * 30 + "\n",
        encoding="utf-8",
    )
    assembly.write_text(">chr1\n" + "ACGT" * 80 + "\n", encoding="utf-8")
    copy_number.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
        "haploid_depth\testimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF000001\t120\t10\t10\t1\t10\t1200\thigh\t\n"
        "TXF000002\t120\t10\t10\t1\t10\t1200\thigh\t\n",
        encoding="utf-8",
    )
    arrays.write_text("chr1\t0\t320\tTXF000002\t1000\t.\thigh\t\n", encoding="utf-8")

    candidates, _ = rank_toy_probes(
        ProbeConfig(
            monomers=monomers,
            assembly=assembly,
            copy_number=copy_number,
            arrays=arrays,
            outdir=outdir,
            min_len=80,
            max_len=120,
        )
    )

    assert candidates
    assert all(candidate.family_id != "TXF000001" for candidate in candidates)
    assert candidates == sorted(candidates, key=lambda item: (-item.probe_score, item.family_id, item.probe_id))
