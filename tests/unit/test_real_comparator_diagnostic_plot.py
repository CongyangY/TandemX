import argparse
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.plot_real_comparator_diagnostics import load_run, parse_run, plot


def make_archive(root: Path, run_id: str, total_bases: int, read_count: int) -> Path:
    path = root / run_id
    path.mkdir()
    environment = {
        "accuracy": "not_assessed_without_curated_independent_truth",
        "resources": "direct-child wait4; acquisition may overlap; not publication ranking",
        "threads": 1,
        "repetitions": 1,
        "input": {"total_bases": total_bases, "read_count": read_count},
    }
    (path / "environment.json").write_text(json.dumps(environment) + "\n")
    rows = []
    for index, tool in enumerate(("tandemx", "trf", "tidehunter"), 1):
        union = 100 * index
        rows.append({
            "tool": tool,
            "exit_code": 0,
            "runtime_seconds": 10 * index,
            "peak_rss_mib": 100 * index,
            "timed_out": False,
            "observed_in_scope_calls": 2 * index,
            "observed_positive_reads": index,
            "observed_union_bp": union,
            "observed_union_base_fraction": union / total_bases,
            "normalization": "ok",
        })
    write_table(path / "summary.tsv", rows, list(rows[0]))
    manifest = [
        {"file": name, "sha256": digest_file(path / name), "bytes": (path / name).stat().st_size}
        for name in ("environment.json", "summary.tsv")
    ]
    (path / "archive_manifest.json").write_text(json.dumps(manifest) + "\n")
    return path


def test_real_diagnostic_plot_is_vector_and_source_backed(tmp_path: Path) -> None:
    small = make_archive(tmp_path, "small", 1_000, 10)
    large = make_archive(tmp_path, "large", 10_000, 100)
    outdir = tmp_path / "figure"

    plot([("Material A", small), ("Material A", large)], outdir)

    provenance = json.loads((outdir / "figure_provenance.json").read_text())
    assert provenance["complete"] and provenance["run_count"] == 2
    svg = (outdir / "real_comparator_diagnostics.svg").read_text()
    assert "<image" not in svg and "One-thread real-read diagnostics" in svg
    assert len((outdir / "panel_source.tsv").read_text().splitlines()) == 7
    with pytest.raises(ValueError, match="already exists"):
        plot([("Material A", small), ("Material A", large)], outdir)


def test_real_diagnostic_plot_rejects_hash_drift_and_invalid_spec(tmp_path: Path) -> None:
    archive = make_archive(tmp_path, "run", 1_000, 10)
    (archive / "summary.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="hash/size differs"):
        load_run("Material", archive)
    with pytest.raises(argparse.ArgumentTypeError):
        parse_run("missing-separator")
