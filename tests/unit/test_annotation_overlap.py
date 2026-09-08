import json
from pathlib import Path

import pytest

from benchmarks.retrospective.annotation_overlap import audit_annotation_overlap, run_annotation_audit
from benchmarks.retrospective.raster_validate import validate_rasterized_overlap, write_validation_receipt


def _write_fixture(tmp_path: Path) -> tuple[Path, Path]:
    arrays = tmp_path / "arrays.bed"
    arrays.write_text(
        "chr1\t0\t100\tf1\n"
        "chr1\t50\t150\tf1\n"
        "chr1\t200\t250\tf2\n"
        "chr2\t0\t20\tf3\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "repeats.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\tauthor\tcentromere\t51\t120\t.\t+\t.\tID=c1\n"
        "chr1\tauthor\tcentromere\t111\t170\t.\t+\t.\tID=c2\n"
        "chr1\tauthor\t5S_rDNA\t221\t230\t.\t+\t.\tID=r1\n"
        "chr2\tauthor\t45S_rDNA\t1\t5\t.\t+\t.\tID=r2\n"
        "chr2\tauthor\ttelomere\t18\t30\t.\t+\t.\tID=t1\n",
        encoding="utf-8",
    )
    return arrays, annotation


def test_overlap_audit_unions_coordinates_and_retains_all_families(tmp_path: Path) -> None:
    arrays, annotation = _write_fixture(tmp_path)
    classes, pairs, families, summary = audit_annotation_overlap(arrays, annotation)
    by_class = {row["annotation_class"]: row for row in classes}
    by_family = {row["family_id"]: row for row in families}
    assert by_class["centromere"]["annotation_bp"] == 120
    assert by_class["centromere"]["tandemx_overlap_bp"] == 100
    assert by_class["5S_rDNA"]["tandemx_overlap_bp"] == 10
    assert by_family["f1"]["tandemx_array_bp"] == 150
    assert by_family["f1"]["selected_annotation_overlap_bp"] == 100
    assert by_family["f2"]["best_annotation_class"] == "5S_rDNA"
    assert by_family["f3"]["selected_annotation_overlap_bp"] == 8
    assert len(pairs) == 4
    assert summary["tandemx_array_union_bp"] == 220
    assert summary["overlap_union_bp"] == 118


def test_run_writes_outputs_and_refuses_overwrite(tmp_path: Path) -> None:
    arrays, annotation = _write_fixture(tmp_path)
    outdir = tmp_path / "audit"
    summary = run_annotation_audit(arrays, annotation, outdir)
    assert json.loads((outdir / "summary.json").read_text()) == summary
    assert (outdir / "annotation_class_summary.tsv").read_text().startswith("annotation_class\t")
    with pytest.raises(FileExistsError):
        run_annotation_audit(arrays, annotation, outdir)


def test_audit_rejects_missing_selected_class_and_malformed_gff(tmp_path: Path) -> None:
    arrays, annotation = _write_fixture(tmp_path)
    with pytest.raises(ValueError, match="absent"):
        audit_annotation_overlap(arrays, annotation, ("not_present",))
    annotation.write_text("chr1\tauthor\tcentromere\t0\t10\t.\t+\t.\tID=x\n")
    with pytest.raises(ValueError, match="invalid"):
        audit_annotation_overlap(arrays, annotation, ("centromere",))


def test_independent_raster_validation_matches_and_refuses_receipt_overwrite(tmp_path: Path) -> None:
    arrays, annotation = _write_fixture(tmp_path)
    outdir = tmp_path / "audit"
    run_annotation_audit(arrays, annotation, outdir)
    receipt = validate_rasterized_overlap(
        arrays, annotation, outdir, ("centromere", "5S_rDNA", "45S_rDNA", "telomere")
    )
    assert receipt["complete"] is True
    assert receipt["observed"]["overlap_union_bp"] == 118
    receipt_path = outdir / "independent_validation.json"
    write_validation_receipt(receipt_path, receipt)
    with pytest.raises(FileExistsError):
        write_validation_receipt(receipt_path, receipt)
