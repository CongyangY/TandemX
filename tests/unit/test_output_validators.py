from __future__ import annotations

from pathlib import Path

import pytest

from tandemx.io.validators import ValidationError, validate_arrays_bed, validate_project, validate_tsv


def write_valid_project(project: Path) -> None:
    project.mkdir()
    (project / "candidate_reads.tsv").write_text(
        "read_id\tcandidate_id\tread_start\tread_end\tstrand\tperiod_bp\trepeat_span_bp\t"
        "unit_count\tscore\tlow_complexity_flag\tconfidence\twarning\n"
        "r1\tTXC000001\t0\t400\t+\t200\t400\t2.0\t0.95\tfalse\thigh\t\n",
        encoding="utf-8",
    )
    (project / "families.tsv").write_text(
        "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\t"
        "support_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\n"
        "TXF000001\tTXM000001\t200\tabc\t0.5\t3\t1200\t0.95\tfalse\thigh\t\n",
        encoding="utf-8",
    )
    (project / "monomers.fa").write_text(
        ">family_id=TXF000001;monomer_id=TXM000001;length_bp=200;confidence=high\nACGT\n",
        encoding="utf-8",
    )
    (project / "copy_number.tsv").write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\thaploid_depth\t"
        "estimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF000001\t200\t20\t10\t2\t5\t1000\thigh\t\n",
        encoding="utf-8",
    )
    (project / "repeat_density.bedgraph").write_text("chr1\t0\t100\t0.5\n", encoding="utf-8")
    (project / "arrays.bed").write_text("chr1\t0\t200\tTXF000001\t900\t.\thigh\t\n", encoding="utf-8")
    (project / "assembly_vs_read_cn.tsv").write_text(
        "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\n"
        "TXF000001\t1000\t200\t0.2\tpossible_collapse\tmedium\t\n",
        encoding="utf-8",
    )
    (project / "probes.fa").write_text(
        ">probe_id=TXP000001;family_id=TXF000001;length_bp=80;probe_score=0.9000;confidence=high\nACGT\n",
        encoding="utf-8",
    )
    (project / "probes.rank.tsv").write_text(
        "probe_id\tfamily_id\tsequence_length\tgc_content\ttm\testimated_copy_number\tarrayiness_score\t"
        "specificity_score\toff_target_hits\tpredicted_regions\tprobe_score\tconfidence\twarning\n"
        "TXP000001\tTXF000001\t80\t0.5\t240\t5\t1\t1\t0\tchr1:0-200\t0.9\thigh\t\n",
        encoding="utf-8",
    )
    (project / "in_silico_fish.tsv").write_text(
        "probe_id\tchrom\tstart\tend\tpredicted_signal\tconfidence\twarning\n"
        "TXP000001\tchr1\t0\t200\t0.9\thigh\t\n",
        encoding="utf-8",
    )


def test_validate_project_accepts_core_outputs(tmp_path: Path) -> None:
    project = tmp_path / "project"
    write_valid_project(project)

    results = validate_project(project)

    names = {result.path.name for result in results}
    assert {
        "candidate_reads.tsv",
        "families.tsv",
        "copy_number.tsv",
        "repeat_density.bedgraph",
        "arrays.bed",
        "assembly_vs_read_cn.tsv",
        "probes.rank.tsv",
        "in_silico_fish.tsv",
        "monomers.fa",
        "probes.fa",
    } <= names


def test_validate_tsv_rejects_missing_required_fields(tmp_path: Path) -> None:
    path = tmp_path / "copy_number.tsv"
    path.write_text("family_id\testimated_copy_number\nTXF000001\t5\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="missing required field"):
        validate_project(tmp_path)


def test_validate_tsv_rejects_non_numeric_values(tmp_path: Path) -> None:
    path = tmp_path / "copy_number.tsv"
    path.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\thaploid_depth\t"
        "estimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF000001\tabc\t1\t1\t1\t1\t100\thigh\t\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="not numeric"):
        validate_project(tmp_path)


def test_validate_bed_rejects_non_half_open_coordinates(tmp_path: Path) -> None:
    path = tmp_path / "arrays.bed"
    path.write_text("chr1\t10\t10\tTXF000001\t900\t.\thigh\t\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="end must be greater than start"):
        validate_arrays_bed(path)


def test_validate_tandemx_fasta_rejects_bad_header(tmp_path: Path) -> None:
    path = tmp_path / "monomers.fa"
    path.write_text(">TXF000001\nACGT\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="invalid TandemX FASTA header"):
        validate_project(tmp_path)


def test_validate_tsv_rejects_header_only_file(tmp_path: Path) -> None:
    path = tmp_path / "copy_number.tsv"
    path.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\thaploid_depth\t"
        "estimated_copy_number\testimated_bp\tconfidence\twarning\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="no records"):
        validate_project(tmp_path)


def test_validate_project_allows_header_only_family_similarity(tmp_path: Path) -> None:
    path = tmp_path / "family_similarity.tsv"
    path.write_text(
        "family_a\tfamily_b\tlength_a_bp\tlength_b_bp\tkmer_jaccard\t"
        "shared_kmer_fraction\tlocal_identity\tlocal_overlap_bp\t"
        "local_overlap_fraction_shorter\tlength_ratio\torientation\t"
        "relationship\tredundant_candidate\tnotes\n",
        encoding="utf-8",
    )

    results = validate_project(tmp_path)

    assert results[0].path.name == "family_similarity.tsv"
    assert results[0].record_count == 0
