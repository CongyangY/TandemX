from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from tandemx.importers.tidehunter import iter_tidehunter_f2
from tandemx.io.validators import validate_project


MONOMER = "ACGTGCACTGAA"


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        text=True,
        capture_output=True,
        check=False,
    )


def write_reads(path: Path) -> Path:
    rotated = MONOMER[4:] + MONOMER[:4]
    path.write_text(
        f">r1;source=test_a\n{MONOMER * 3}ACGT\n"
        f">r2;source=test_b\n{rotated * 3}TGCA\n",
        encoding="utf-8",
    )
    return path


def write_native(path: Path) -> Path:
    rotated = MONOMER[4:] + MONOMER[:4]
    path.write_text(
        f"r1;source=test_a rep0 3 40 1 36 12 99.0 3 1,13,25 {MONOMER}\n"
        f"r2;source=test_b rep0 3 40 1 36 12 98.0 3 1,13,25 {rotated}\n",
        encoding="utf-8",
    )
    return path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_import_tidehunter_builds_a_valid_downstream_catalogue(tmp_path: Path) -> None:
    reads = write_reads(tmp_path / "reads.fa")
    native = write_native(tmp_path / "tidehunter.tsv")
    out = tmp_path / "imported"
    result = run_cli(
        "import",
        "tidehunter",
        "--input",
        str(native),
        "--reads",
        str(reads),
        "--min-support-reads",
        "2",
        "--backend",
        "python",
        "--outdir",
        str(out),
    )
    assert result.returncode == 0, result.stderr
    assert "imported 2 candidates into 1 families" in result.stdout
    candidates = rows(out / "candidate_reads.tsv")
    assert [(row["read_start"], row["read_end"]) for row in candidates] == [
        ("0", "36"),
        ("0", "36"),
    ]
    assert {row["warning"].split(";")[0] for row in candidates} == {
        "external_detector=tidehunter"
    }
    audit = rows(out / "tidehunter_import.tsv")
    assert [row["read_id"] for row in audit] == ["r1", "r2"]
    assert [row["native_read_id"] for row in audit] == [
        "r1;source=test_a", "r2;source=test_b",
    ]
    assert [row["native_start"] for row in audit] == ["1", "1"]
    assert [row["normalized_start"] for row in audit] == ["0", "0"]
    families = rows(out / "families.tsv")
    assert len(families) == 1
    assert families[0]["support_read_count"] == "2"
    assert "external_detector=tidehunter" in families[0]["warning"]
    assert rows(out / "family_hierarchy.tsv") == []
    membership = rows(out / "monomer_membership.tsv")
    assert {row["status"] for row in membership} == {"assigned"}
    summary = json.loads((out / "import_summary.json").read_text())
    assert summary["complete"] and summary["source_detector"] == "TideHunter"
    assert summary["read_count"] == 2 and summary["total_bases"] == 80
    assert summary["candidate_count"] == 2 and summary["family_count"] == 1
    assert len(summary["reads_semantic_sha256"]) == 64
    validated = {result.path.name for result in validate_project(out)}
    assert {"candidate_reads.tsv", "families.tsv", "family_hierarchy.tsv", "tidehunter_import.tsv"} <= validated
    assert 'status: "tidehunter_import_completed"' in (out / "run_config.yaml").read_text()

    quantify = run_cli(
        "quantify",
        "--reads",
        str(reads),
        "--catalog",
        str(out / "monomers.fa"),
        "--genome-size",
        "80",
        "--k",
        "5",
        "--outdir",
        str(tmp_path / "quantify"),
    )
    assert quantify.returncode == 0, quantify.stderr
    assert rows(tmp_path / "quantify" / "copy_number.tsv")


def test_import_tidehunter_empty_native_output_is_a_verified_negative(
    tmp_path: Path,
) -> None:
    reads = write_reads(tmp_path / "reads.fa")
    native = tmp_path / "tidehunter.tsv"
    native.write_text("", encoding="utf-8")
    out = tmp_path / "imported"
    result = run_cli(
        "import",
        "tidehunter",
        "--input",
        str(native),
        "--reads",
        str(reads),
        "--min-support-reads",
        "1",
        "--backend",
        "python",
        "--outdir",
        str(out),
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads((out / "import_summary.json").read_text())
    assert summary["candidate_count"] == 0 and summary["family_count"] == 0
    assert (out / "monomers.fa").stat().st_size == 0
    assert validate_project(out)


@pytest.mark.parametrize(
    "native_row, message",
    [
        (
            f"missing rep0 3 40 1 36 12 99 3 1,13,25 {MONOMER}\n",
            "read absent from source inputs",
        ),
        (
            f"r1;source=test_a rep0 3 41 1 36 12 99 3 1,13,25 {MONOMER}\n",
            "read length differs from source input",
        ),
        (
            f"r1;source=test_a rep0 3 40 1 36 11 99 3 1,13,25 {MONOMER}\n",
            "consensus length differs from field 7",
        ),
        (
            f"r1;source=test_a rep0 3 40 5 36 12 99 3 1,13,25 {MONOMER}\n",
            "Invalid TideHunter -f 2 values",
        ),
        ("broken row\n", "11 fields"),
    ],
)
def test_import_tidehunter_rejects_invalid_or_unverifiable_calls(
    tmp_path: Path, native_row: str, message: str
) -> None:
    reads = write_reads(tmp_path / "reads.fa")
    native = tmp_path / "tidehunter.tsv"
    native.write_text(native_row, encoding="utf-8")
    out = tmp_path / "imported"
    result = run_cli(
        "import",
        "tidehunter",
        "--input",
        str(native),
        "--reads",
        str(reads),
        "--backend",
        "python",
        "--outdir",
        str(out),
    )
    assert result.returncode == 2
    assert message in result.stderr
    assert not list(out.glob(".tidehunter-read-index-*.sqlite3"))
    assert not (out / "candidate_reads.tsv").exists()


def test_tidehunter_parser_rejects_late_malformed_rows(tmp_path: Path) -> None:
    native = write_native(tmp_path / "tidehunter.tsv")
    with native.open("a", encoding="utf-8") as handle:
        handle.write("late malformed row\n")
    iterator = iter_tidehunter_f2(native)
    assert next(iterator).read_id == "r1;source=test_a"
    assert next(iterator).read_id == "r2;source=test_b"
    with pytest.raises(ValueError, match="11 fields"):
        next(iterator)
