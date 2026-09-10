import pytest

from benchmarks.challenge.unified_native_io import (
    load_tandemx_arrays,
    load_tandemx_catalogue,
    mapping_paf_to_unique_arrays,
    mapping_paf_to_all_arrays,
)


def _write(path, text):
    path.write_text(text, encoding="utf-8")


def _discover(tmp_path):
    _write(tmp_path / "monomers.fa", ">family_id=TXF1;monomer_id=TXM1\nACGT\n")
    _write(tmp_path / "candidate_reads.tsv", "read_id\tcandidate_id\tread_start\tread_end\tperiod_bp\nr1\tc1\t10\t30\t4\n")
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\nr1\tc1\tTXF1\tassigned\n")


def test_load_catalogue_and_strict_candidate_membership_join(tmp_path):
    _discover(tmp_path)
    assert load_tandemx_catalogue(tmp_path / "monomers.fa") == {"TXF1": "ACGT"}
    arrays = load_tandemx_arrays(tmp_path)
    assert [(row.read_id, row.start, row.end, row.sequence, row.family_id) for row in arrays] == [("r1", 10, 30, "ACGT", "TXF1")]


@pytest.mark.parametrize("membership", ["r1\tc1\tUNKNOWN\tassigned\n", "r1\tc1\tTXF1\tassigned\nr1\tc1\tTXF1\tassigned\n"])
def test_load_arrays_rejects_unknown_or_multiple_memberships(tmp_path, membership):
    _discover(tmp_path)
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\n" + membership)
    with pytest.raises(ValueError):
        load_tandemx_arrays(tmp_path)


def test_load_arrays_rejects_unmatched_candidate(tmp_path):
    _discover(tmp_path)
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\n")
    with pytest.raises(ValueError, match="Unmatched"):
        load_tandemx_arrays(tmp_path)


def test_load_arrays_preserves_legitimate_unassigned_candidate(tmp_path):
    _discover(tmp_path)
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\nr1\tc1\tNA\tunresolved_sequence\n")
    array = load_tandemx_arrays(tmp_path)[0]
    assert array.sequence == ""
    assert array.family_id == "unassigned_candidate:c1@r1"


def test_load_arrays_rejects_malformed_na_status_and_accepts_empty_catalogue(tmp_path):
    _discover(tmp_path)
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\nr1\tc1\tNA\tassigned\n")
    with pytest.raises(ValueError, match="Unknown"):
        load_tandemx_arrays(tmp_path)
    _write(tmp_path / "monomers.fa", "")
    _write(tmp_path / "candidate_reads.tsv", "read_id\tcandidate_id\tread_start\tread_end\tperiod_bp\n")
    _write(tmp_path / "monomer_membership.tsv", "read_id\tcandidate_id\tfamily_id\tstatus\n")
    assert load_tandemx_arrays(tmp_path) == []


def test_paf_keeps_primary_secondary_then_excludes_cross_family_overlap(tmp_path):
    paf = tmp_path / "maps.paf"
    _write(
        paf,
        "r1\t1000\t0\t200\t+\tN1\t10000\t0\t200\t190\t200\t60\ttp:A:P\n"
        "r1\t1000\t20\t200\t+\tN1\t10000\t0\t180\t171\t180\t20\ttp:A:S\n"
        "r1\t1000\t100\t250\t+\tN2\t10000\t0\t150\t143\t150\t60\ttp:A:S\n",
    )
    catalogue = {"N1": "ACGT", "N2": "TGCA"}
    raw = mapping_paf_to_all_arrays(paf, catalogue)
    arrays, ambiguous_bp = mapping_paf_to_unique_arrays(paf, catalogue)
    assert [(a.family_id, a.start, a.end) for a in raw] == [("N1", 0, 200), ("N1", 20, 200), ("N2", 100, 250)]
    assert ambiguous_bp == 100
    assert [(a.family_id, a.start, a.end) for a in arrays] == [("N1", 0, 100), ("N2", 200, 250)]


def test_paf_rejects_unknown_family_and_bad_query_coordinates(tmp_path):
    paf = tmp_path / "maps.paf"
    _write(paf, "r1\t100\t40\t20\t+\tUNKNOWN\t100\t0\t20\t20\t20\t60\n")
    with pytest.raises(ValueError, match="query coordinates"):
        mapping_paf_to_unique_arrays(paf, {"N1": "ACGT"})
    _write(paf, "r1\t100\t0\t20\t+\tUNKNOWN\t100\t0\t20\t20\t20\t60\n")
    with pytest.raises(ValueError, match="Unknown PAF"):
        mapping_paf_to_unique_arrays(paf, {"N1": "ACGT"})
