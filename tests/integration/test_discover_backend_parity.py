from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest
import tandemx

from tandemx.cli import main
from tandemx.discover.rust_backend import rust_backend_available
from tandemx.utils.threads import discover_thread_limit


pytestmark = pytest.mark.skipif(
    not rust_backend_available(),
    reason="TandemX Rust extension is not installed",
)


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_python_and_rust_discover_backend_parity(tmp_path: Path) -> None:
    toy = tmp_path / "toy"
    simulated = run_cli(
        "simulate",
        "toy",
        "--outdir",
        str(toy),
        "--seed",
        "909",
        "--num-reads",
        "300",
        "--read-length",
        "1800",
        "--monomer-lengths",
        "421,729",
        "--copies",
        "8,6",
        "--error-rate",
        "0.005",
    )
    assert simulated.returncode == 0, simulated.stderr

    family_lengths: dict[str, list[int]] = {}
    candidate_rows: dict[str, list[dict[str, str]]] = {}
    family_rows: dict[str, list[dict[str, str]]] = {}
    for backend in ("python", "rust"):
        outdir = tmp_path / backend
        result = run_cli(
            "discover",
            "--reads",
            str(toy / "reads.fa"),
            "--outdir",
            str(outdir),
            "--min-period",
            "300",
            "--max-period",
            "760",
            "--min-repeat-span",
            "800",
            "--min-support-reads",
            "3",
            "--kmer-backend",
            backend,
        )
        assert result.returncode == 0, result.stderr
        candidate_rows[backend] = read_rows(outdir / "candidate_reads.tsv")
        family_rows[backend] = read_rows(outdir / "families.tsv")
        family_lengths[backend] = sorted(
            int(row["monomer_length_bp"])
            for row in read_rows(outdir / "families.tsv")
        )

    assert candidate_rows["python"] == candidate_rows["rust"]
    assert family_rows["python"] == family_rows["rust"]
    assert family_lengths["python"] == family_lengths["rust"]


def test_rust_discover_threads_preserve_outputs(tmp_path: Path) -> None:
    if discover_thread_limit() < 2:
        pytest.skip("host thread policy allows only one discover thread")

    reads = tmp_path / "reads.fa"
    monomer = "ACGTTCAGGACTAACCGTGA"
    sequence = monomer * 30
    reads.write_text(
        "".join(f">read_{index}\n{sequence}\n" for index in range(1, 81)),
        encoding="utf-8",
    )

    outputs: dict[str, list[dict[str, str]]] = {}
    for label, threads in (("single", "1"), ("parallel", "2")):
        outdir = tmp_path / label
        result = run_cli(
            "discover",
            "--reads",
            str(reads),
            "--outdir",
            str(outdir),
            "--min-period",
            "20",
            "--max-period",
            "30",
            "--min-support-reads",
            "2",
            "--min-repeat-span",
            "100",
            "--kmer-backend",
            "rust",
            "--threads",
            threads,
            "--chunk-size",
            "7",
        )
        assert result.returncode == 0, result.stderr
        outputs[label] = read_rows(outdir / "candidate_reads.tsv")

    assert outputs["parallel"] == outputs["single"]
    parallel_log = (tmp_path / "parallel" / "run.log").read_text(encoding="utf-8")
    assert "scan_threads requested=2 effective=2 parallel=true parallel_files=false backend=rust" in parallel_log


def test_rust_discover_multiple_files_preserve_outputs(tmp_path: Path) -> None:
    if discover_thread_limit() < 2:
        pytest.skip("host thread policy allows only one discover thread")

    first = tmp_path / "reads_a.fa"
    second = tmp_path / "reads_b.fa"
    merged = tmp_path / "reads_merged.fa"
    monomer = "ACGTTCAGGACTAACCGTGA"
    sequence = monomer * 30
    first.write_text(
        "".join(f">read_a_{index}\n{sequence}\n" for index in range(1, 41)),
        encoding="utf-8",
    )
    second.write_text(
        "".join(f">read_b_{index}\n{sequence}\n" for index in range(1, 41)),
        encoding="utf-8",
    )
    merged.write_text(first.read_text(encoding="utf-8") + second.read_text(encoding="utf-8"), encoding="utf-8")

    outputs: dict[str, list[dict[str, str]]] = {}
    for label, reads in (("merged", (str(merged),)), ("split", (str(first), str(second)))):
        outdir = tmp_path / label
        result = run_cli(
            "discover",
            "--reads",
            *reads,
            "--outdir",
            str(outdir),
            "--min-period",
            "20",
            "--max-period",
            "30",
            "--min-support-reads",
            "2",
            "--min-repeat-span",
            "100",
            "--kmer-backend",
            "rust",
            "--threads",
            "2",
            "--chunk-size",
            "7",
        )
        assert result.returncode == 0, result.stderr
        outputs[label] = read_rows(outdir / "candidate_reads.tsv")

    assert outputs["split"] == outputs["merged"]
    split_log = (tmp_path / "split" / "run.log").read_text(encoding="utf-8")
    assert "parallel_files=false" in split_log


def test_discover_auto_backend_uses_rust_threads_by_default(tmp_path: Path) -> None:
    if discover_thread_limit() < 2:
        pytest.skip("host thread policy allows only one discover thread")

    reads = tmp_path / "reads.fa"
    sequence = "ACGTTCAGGACTAACCGTGA" * 20
    reads.write_text(
        "".join(f">read_{index}\n{sequence}\n" for index in range(1, 8)),
        encoding="utf-8",
    )
    outdir = tmp_path / "discover"

    result = run_cli(
        "discover",
        "--reads",
        str(reads),
        "--outdir",
        str(outdir),
        "--min-period",
        "20",
        "--max-period",
        "30",
        "--min-support-reads",
        "2",
        "--min-repeat-span",
        "100",
    )

    assert result.returncode == 0, result.stderr
    log = (outdir / "run.log").read_text(encoding="utf-8")
    assert "scan_threads requested=" in log
    assert "parallel=true parallel_files=false backend=rust" in log
    config = (outdir / "run_config.yaml").read_text(encoding="utf-8")
    assert 'kmer_backend: "rust"' in config


def test_rust_cli_reports_unavailable_extension(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read_1\nACGTTCAGGACTAACCGTGAACGTTCAGGACTAACCGTGA\n", encoding="utf-8")
    outdir = tmp_path / "discover"
    monkeypatch.setitem(sys.modules, "tandemx._rust_core", None)
    monkeypatch.delattr(tandemx, "_rust_core", raising=False)

    with pytest.raises(SystemExit) as error:
        main(
            [
                "discover",
                "--reads",
                str(reads),
                "--outdir",
                str(outdir),
                "--min-period",
                "20",
                "--max-period",
                "30",
                "--min-support-reads",
                "1",
                "--min-repeat-span",
                "20",
                "--kmer-backend",
                "rust",
            ]
        )

    assert error.value.code == 2
    assert "Rust backend is unavailable" in capsys.readouterr().err
    assert "rust_backend_unavailable=" in (outdir / "run.log").read_text(encoding="utf-8")
