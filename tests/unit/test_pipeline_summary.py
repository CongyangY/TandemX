from __future__ import annotations

import csv
import json
from pathlib import Path

from tandemx.pipeline import PipelineConfig, StepRecord, finalize_run_outputs, write_summaries


def pipeline_config(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        reads=(tmp_path / "reads.fa",),
        assembly=None,
        genome_size=1_000_000,
        haploid_depth=None,
        outdir=tmp_path / "run",
        max_reads=10,
        max_read_bases=None,
        kmer_backend="python",
        steps=("discover", "quantify"),
        min_period=20,
        max_period=2000,
        top_periods=5,
        threads=1,
        resume=False,
        force=False,
        profile=False,
    )


def test_pipeline_summary_writes_tsv_and_json(tmp_path: Path) -> None:
    config = pipeline_config(tmp_path)
    records = [
        StepRecord(
            run_id="run-1",
            input_reads=str(config.reads),
            input_assembly="",
            max_reads=10,
            max_read_bases=None,
            kmer_backend="python",
            step="discover",
            command="python -m tandemx.cli discover",
            start_time="2026-01-01T00:00:00+00:00",
            end_time="2026-01-01T00:00:01+00:00",
            runtime_seconds=1.25,
            exit_status=0,
            output_dir=str(config.outdir / "discover"),
            output_validated=True,
            notes="",
        ),
        StepRecord(
            run_id="run-1",
            input_reads=str(config.reads),
            input_assembly="",
            max_reads=10,
            max_read_bases=None,
            kmer_backend="python",
            step="quantify",
            command="python -m tandemx.cli quantify",
            start_time="2026-01-01T00:00:01+00:00",
            end_time="2026-01-01T00:00:02+00:00",
            runtime_seconds=0.75,
            exit_status=2,
            output_dir=str(config.outdir / "quantify"),
            output_validated=False,
            notes="failed_step",
        ),
    ]

    write_summaries(config, records)

    with (config.outdir / "pipeline_summary.tsv").open(encoding="utf-8") as handle:
        tsv_rows = list(csv.DictReader(handle, delimiter="\t"))
    json_rows = json.loads((config.outdir / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert [row["step"] for row in tsv_rows] == ["discover", "quantify"]
    assert float(tsv_rows[0]["runtime_seconds"]) == 1.25
    assert tsv_rows[1]["exit_status"] == "2"
    assert json_rows[1]["notes"] == "failed_step"
    assert json_rows[1]["output_validated"] is False


def test_run_report_summarizes_existing_repeat_annotation(tmp_path: Path) -> None:
    config = pipeline_config(tmp_path)
    (config.outdir / "discover").mkdir(parents=True)
    (config.outdir / "repeat_annotation.tsv").write_text(
        "family_id\tmonomer_length\tbest_known_id\tbest_known_length\tbest_orientation\t"
        "shared_kmer_fraction\tjaccard\tdice\tcontainment_discovered_in_known\t"
        "containment_known_in_discovered\tlocal_identity\tlocal_overlap_bp\tannotation_status\tnotes\n"
        "TXF000001\t120\tknown\t120\tforward\t1\t1\t1\t1\t1\t1\t120\tstrong_known_match\tpost hoc\n",
        encoding="utf-8",
    )
    records = [
        StepRecord(
            run_id="run-1",
            input_reads=str(config.reads),
            input_assembly="",
            max_reads=10,
            max_read_bases=None,
            kmer_backend="python",
            step="discover",
            command="python -m tandemx.cli discover",
            start_time="2026-01-01T00:00:00+00:00",
            end_time="2026-01-01T00:00:01+00:00",
            runtime_seconds=1.0,
            exit_status=0,
            output_dir=str(config.outdir / "discover"),
            output_validated=True,
            notes="",
        )
    ]

    finalize_run_outputs(config, records)

    report = (config.outdir / "run_report.md").read_text(encoding="utf-8")
    assert "Repeat annotation summary: strong_known_match=1" in report
    manifest = (config.outdir / "output_manifest.tsv").read_text(encoding="utf-8")
    assert "repeat_annotation.tsv\ttrue" in manifest
