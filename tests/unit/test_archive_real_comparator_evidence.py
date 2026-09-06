import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import write_table
from benchmarks.scripts.archive_real_comparator_evidence import ROOT_FILES, TOOL_FILES, archive


def make_run(root: Path) -> None:
    root.mkdir()
    (root / "environment.json").write_text(
        json.dumps({"accuracy": "not_assessed_without_curated_independent_truth"}) + "\n"
    )
    rows = [
        {
            "tool": tool,
            "exit_code": 0,
            "runtime_seconds": index + 1,
            "peak_rss_mib": 100 + index,
            "timed_out": False,
            "normalization": "ok",
        }
        for index, tool in enumerate(TOOL_FILES)
    ]
    write_table(root / "summary.tsv", rows, list(rows[0]))
    for tool, names in TOOL_FILES.items():
        for name in names:
            path = root / tool / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"{tool}:{name}\n")


def test_archive_real_comparator_copies_only_declared_compact_evidence(tmp_path: Path) -> None:
    source = tmp_path / "source"
    make_run(source)
    (source / "reads.fa").write_text(">large\nACGT\n")
    outdir = tmp_path / "archive"

    manifest = archive(source, outdir)

    expected = len(ROOT_FILES) + sum(len(names) for names in TOOL_FILES.values())
    assert len(manifest) == expected
    assert not (outdir / "reads.fa").exists()
    assert json.loads((outdir / "archive_manifest.json").read_text()) == manifest
    with pytest.raises(ValueError, match="already exists"):
        archive(source, outdir)


def test_archive_real_comparator_rejects_incomplete_run(tmp_path: Path) -> None:
    source = tmp_path / "source"
    make_run(source)
    rows = [
        {
            "tool": "tandemx",
            "exit_code": 1,
            "runtime_seconds": 1,
            "peak_rss_mib": 100,
            "timed_out": False,
            "normalization": "not_run",
        }
    ]
    write_table(source / "summary.tsv", rows, list(rows[0]))
    with pytest.raises(ValueError, match="three successful"):
        archive(source, tmp_path / "archive")
