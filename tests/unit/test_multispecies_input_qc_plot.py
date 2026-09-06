import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.plot_multispecies_input_qc import aggregate_histogram, plot


def make_library(root: Path, run: str, species: str, material: str, gc: float) -> dict:
    qcdir = root / run
    qcdir.mkdir(parents=True)
    receipt = {
        "complete": True,
        "read_count": 3,
        "total_bases": 400,
        "median_length": 100,
        "read_n50": 200,
        "gc_fraction": gc,
        "input_sha256": run.lower().ljust(64, "0"),
    }
    receipt_path = qcdir / "qc.json"
    receipt_path.write_text(json.dumps(receipt) + "\n")
    (qcdir / "length_histogram.tsv").write_text("length_bp\tread_count\n100\t2\n200\t1\n")
    (qcdir / "joint_distribution.tsv").write_text(
        "length_bin_kb\tgc_bin_percent\tmean_quality_bin_phred\tread_count\n"
        "0\t40\t20\t2\n0\t42\t25\t1\n"
    )
    return {
        "reported_species": species,
        "reported_material": material,
        "run_accession": run,
        "read_count": 3,
        "total_bases": 400,
        "median_read_length": 100,
        "read_n50": 200,
        "gc_fraction": gc,
        "raw_sha256": receipt["input_sha256"],
        "qc_receipt": str(receipt_path.relative_to(root)),
        "qc_receipt_sha256": digest_file(receipt_path),
    }


def test_multispecies_input_qc_plot_has_vector_heatmaps_and_exact_provenance(tmp_path: Path) -> None:
    rows = [
        make_library(tmp_path, "RUN2", "Zea mays", "Mo17", 0.46),
        make_library(tmp_path, "RUN1", "Arabidopsis thaliana", "Col-0", 0.37),
    ]
    cohort = tmp_path / "cohort.tsv"
    write_table(cohort, rows, list(rows[0]))

    outdir = tmp_path / "figure"
    plot(cohort, outdir, tmp_path)

    provenance = json.loads((outdir / "figure_provenance.json").read_text())
    assert provenance["complete"]
    assert provenance["library_order"] == ["RUN1", "RUN2"]
    assert provenance["read_count"] == 6 and provenance["total_bases"] == 800
    svg = (outdir / "multispecies_input_qc.svg").read_text()
    assert "<image" not in svg and "Cross-cohort QC" in svg
    with pytest.raises(ValueError, match="already exists"):
        plot(cohort, outdir, tmp_path)


def test_multispecies_input_qc_rejects_hash_drift_and_bins_overflow(tmp_path: Path) -> None:
    row = make_library(tmp_path, "RUN1", "Oryza sativa", "Nipponbare", 0.44)
    row["qc_receipt_sha256"] = "0" * 64
    cohort = tmp_path / "cohort.tsv"
    write_table(cohort, [row], list(row))
    with pytest.raises(ValueError, match="hash differs"):
        plot(cohort, tmp_path / "figure", tmp_path)

    bins, counts = aggregate_histogram(
        {1: 2, 20: 3, 22: 4, 99: 5}, width=2, lower=20, upper=24
    )
    assert bins == [20, 22, 24] and counts == [5, 4, 5]
