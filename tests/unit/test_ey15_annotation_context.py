import csv
import json
from pathlib import Path

from benchmarks.scripts.summarize_ey15_annotation_context import (
    fisher_right_tail,
    summarize,
    write_outputs,
)


def write_bed(path: Path, lengths: dict[str, int]) -> None:
    path.write_text(
        "".join(
            f"chr1\t0\t{length}\t{family}\n"
            for family, length in lengths.items()
            if length > 0
        )
    )


def test_annotation_context_retains_all_eligible_families(tmp_path: Path) -> None:
    old = tmp_path / "old.bed"
    new = tmp_path / "new.bed"
    annotation = tmp_path / "annotation.tsv"
    write_bed(old, {"collapse_annotated": 10, "collapse_unannotated": 0, "retained": 100})
    write_bed(new, {"collapse_annotated": 100, "collapse_unannotated": 100, "retained": 100})
    with annotation.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "family_id",
                "selected_annotation_overlap_bp",
                "selected_annotation_fraction_of_array",
                "best_annotation_class",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "family_id": "collapse_annotated",
                    "selected_annotation_overlap_bp": 80,
                    "selected_annotation_fraction_of_array": 0.8,
                    "best_annotation_class": "centromere",
                },
                {
                    "family_id": "collapse_unannotated",
                    "selected_annotation_overlap_bp": 0,
                    "selected_annotation_fraction_of_array": 0,
                    "best_annotation_class": "NA",
                },
                {
                    "family_id": "retained",
                    "selected_annotation_overlap_bp": 0,
                    "selected_annotation_fraction_of_array": 0,
                    "best_annotation_class": "NA",
                },
            ]
        )
    rows, summary = summarize(
        old,
        new,
        annotation,
        collapse_threshold=0.6,
        overexpansion_threshold=1.5,
        min_new_bp=50,
    )
    assert len(rows) == 3
    assert summary["contingency"] == {
        "collapse_annotated": 1,
        "collapse_unannotated": 1,
        "other_annotated": 0,
        "other_unannotated": 1,
    }
    assert summary["analysis_status"].startswith("posthoc_descriptive")


def test_fisher_right_tail_and_write_refuse_overwrite(tmp_path: Path) -> None:
    assert fisher_right_tail(8, 0, 3, 8) == 0.0021830594586012544
    outdir = tmp_path / "result"
    write_outputs([{"family_id": "TXF1"}], {"complete": True}, outdir)
    assert json.loads((outdir / "summary.json").read_text()) == {"complete": True}
    try:
        write_outputs([{"family_id": "TXF1"}], {"complete": True}, outdir)
    except FileExistsError:
        pass
    else:
        raise AssertionError("existing result directory was overwritten")
