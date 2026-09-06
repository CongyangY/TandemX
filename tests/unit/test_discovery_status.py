from pathlib import Path

from tandemx.discover.status import FILES, has_verified_empty_catalog, verified_empty_output, write_discovery_summary


def test_verified_zero_requires_intact_output_and_completion(tmp_path: Path) -> None:
    for name in FILES:
        (tmp_path / name).write_text("" if name.endswith(".fa") else "header\n")
    assert not has_verified_empty_catalog(tmp_path)
    write_discovery_summary(tmp_path, 10, 1000, 0, 0)
    assert has_verified_empty_catalog(tmp_path)
    assert verified_empty_output(tmp_path / "candidate_reads.tsv")
    (tmp_path / "families.tsv").write_text("tampered\n")
    assert not has_verified_empty_catalog(tmp_path)


def test_supported_candidates_below_family_threshold(tmp_path: Path) -> None:
    for name in FILES:
        (tmp_path / name).write_text("" if name.endswith(".fa") else "header\n")
    write_discovery_summary(tmp_path, 10, 1000, 1, 0)
    assert has_verified_empty_catalog(tmp_path)
    assert not verified_empty_output(tmp_path / "candidate_reads.tsv")


def test_broken_receipt_does_not_allow_empty_data(tmp_path: Path) -> None:
    (tmp_path / "monomers.fa").write_text("")
    (tmp_path / "discovery_summary.json").write_text('{"status": "no_families"}')
    assert not has_verified_empty_catalog(tmp_path)
