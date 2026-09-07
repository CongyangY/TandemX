import csv
import json
from pathlib import Path
import subprocess
import sys

from tandemx.io.validators import validate_project


COPY_HEADER = (
    "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
    "haploid_depth\testimated_copy_number\testimated_bp\tdepth_mad\t"
    "copy_number_interval_low\tcopy_number_interval_high\tconfidence\twarning\n"
)
COMPARE_HEADER = (
    "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\t"
    "status\tconfidence\twarning\n"
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_cohort_cli_builds_pan_catalogue_and_na_aware_matrices(tmp_path: Path) -> None:
    a = "ACGTGCACTGATCGTACGATCGTAGCTAGCTA"
    b = "TTGACCGTACCGATGCTAGGCTAACGTTCGGA"
    c = "GGCATTCGATACCGTTAGCACGATTCGAGTCA"
    s1 = tmp_path / "s1"
    s2 = tmp_path / "s2"
    s1.mkdir(); s2.mkdir()
    (s1 / "monomers.fa").write_text(f">family_id=A\n{a}\n>family_id=B\n{b}\n")
    (s2 / "monomers.fa").write_text(f">family_id=A2\n{a[7:] + a[:7]}\n>family_id=C\n{c}\n")
    for folder, names, values in ((s1, ("A", "B"), (3200, 1600)), (s2, ("A2", "C"), (3000, 900))):
        lines = [COPY_HEADER]
        for name, value in zip(names, values):
            lines.append(f"{name}\t32\t20\t10\t5\t{value/32}\t{value}\t1\t{value/40}\t{value/25}\tmedium\ttest\n")
        (folder / "copy_number.tsv").write_text("".join(lines))
    (s1 / "comparison.tsv").write_text(
        COMPARE_HEADER
        + "A\t3200\t1600\t0.5\tpossible_collapse\tmedium\t\n"
        + "B\t1600\t1600\t1\tconsistent\tmedium\t\n"
    )
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tmonomers\tcopy_number\tcomparison\n"
        "s1\ts1/monomers.fa\ts1/copy_number.tsv\ts1/comparison.tsv\n"
        "s2\ts2/monomers.fa\ts2/copy_number.tsv\tNA\n"
    )
    out = tmp_path / "out"
    result = subprocess.run([
        sys.executable, "-m", "tandemx.cli", "cohort", "--manifest", str(manifest),
        "--backend", "python", "--outdir", str(out),
    ], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    summary = json.loads((out / "cohort_summary.json").read_text())
    assert summary["complete"] and summary["sample_count"] == 2
    assert summary["local_family_count"] == 4 and summary["pan_family_count"] == 3
    membership = rows(out / "family_membership.tsv")
    pan_a = next(row["pan_family_id"] for row in membership if row["local_family_id"] == "A")
    assert next(row["pan_family_id"] for row in membership if row["local_family_id"] == "A2") == pan_a
    abundance = rows(out / "sample_family_abundance.tsv")
    assert next(row["estimated_bp"] for row in abundance if row["sample_id"] == "s1" and row["pan_family_id"] == pan_a) == "3200.0000"
    assert next(row["estimated_bp"] for row in abundance if row["sample_id"] == "s2" and row["pan_family_id"] == pan_a) == "3000.0000"
    representation = rows(out / "sample_family_representation.tsv")
    assert next(row["assembly_read_ratio"] for row in representation if row["sample_id"] == "s1" and row["pan_family_id"] == pan_a) == "0.500000"
    assert all(row["assembly_read_ratio"] == "NA" for row in representation if row["sample_id"] == "s2")
    assert 'status: "cohort_completed"' in (out / "run_config.yaml").read_text()
    low_matrix = rows(out / "abundance_interval_low_matrix.tsv")
    high_matrix = rows(out / "abundance_interval_high_matrix.tsv")
    assert next(row["s1"] for row in low_matrix if row["pan_family_id"] == pan_a) == "2560.0000"
    assert next(row["s1"] for row in high_matrix if row["pan_family_id"] == pan_a) == "4096.0000"
    input_qc = rows(out / "cohort_input_qc.tsv")
    assert [row["status"] for row in input_qc] == ["complete", "complete"]
    assert input_qc[1]["comparison_sha256"] == "NA"
    assert input_qc[1]["warning"] == "comparison_not_provided"
    assert len(validate_project(out)) == 10


def test_cohort_cli_rejects_copy_number_family_mismatch(tmp_path: Path) -> None:
    sequence = "ACGTGCACTGATCGTACGATCGTAGCTAGCTA"
    for sample_id in ("s1", "s2"):
        folder = tmp_path / sample_id
        folder.mkdir()
        (folder / "monomers.fa").write_text(f">family_id=A\n{sequence}\n")
        family_id = "UNKNOWN" if sample_id == "s1" else "A"
        (folder / "copy_number.tsv").write_text(
            COPY_HEADER
            + f"{family_id}\t32\t20\t10\t5\t100\t3200\t1\t80\t120\tmedium\ttest\n"
        )
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tmonomers\tcopy_number\tcomparison\n"
        "s1\ts1/monomers.fa\ts1/copy_number.tsv\tNA\n"
        "s2\ts2/monomers.fa\ts2/copy_number.tsv\tNA\n"
    )
    result = subprocess.run([
        sys.executable, "-m", "tandemx.cli", "cohort", "--manifest", str(manifest),
        "--backend", "python", "--outdir", str(tmp_path / "out"),
    ], text=True, capture_output=True)
    assert result.returncode == 2
    assert "families absent from its monomer catalogue: UNKNOWN" in result.stderr
