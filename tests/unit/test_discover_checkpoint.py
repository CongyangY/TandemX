from __future__ import annotations

import json
import hashlib
import gzip
from dataclasses import replace
from pathlib import Path

import pytest

from tandemx.discover import checkpoint
from tandemx.discover.mvp import DiscoverConfig, discover_toy_repeats


def config_for(tmp_path: Path, *, fastq_gzip: bool = False) -> DiscoverConfig:
    reads = tmp_path / ("reads.fq.gz" if fastq_gzip else "reads.fa")
    if fastq_gzip:
        with gzip.open(reads, "wt", encoding="utf-8") as handle:
            handle.write("".join(f"@r{i}\n{'ACGT' * 12}\n+\n{'I' * 48}\n" for i in range(7)))
    else:
        reads.write_text("".join(f">r{i}\n{'ACGT' * 12}\n" for i in range(7)), encoding="utf-8")
    return DiscoverConfig(
        reads=reads, outdir=tmp_path / "resumed", min_monomer_len=4,
        max_monomer_len=8, min_support_reads=2, min_repeat_span=16,
        chunk_size=2, checkpoint_every=2, threads=1,
    )


def interrupt_once(config: DiscoverConfig, monkeypatch: pytest.MonkeyPatch) -> None:
    real_save = checkpoint.save

    def after_save(path: Path, body: dict[str, object]) -> None:
        real_save(path, body)
        raise KeyboardInterrupt()

    with monkeypatch.context() as patch:
        patch.setattr(checkpoint, "save", after_save)
        with pytest.raises(KeyboardInterrupt):
            discover_toy_repeats(config)


@pytest.mark.parametrize("fastq_gzip", [False, True])
def test_interrupted_scan_resumes_to_byte_identical_public_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fastq_gzip: bool) -> None:
    config = config_for(tmp_path, fastq_gzip=fastq_gzip)
    interrupt_once(config, monkeypatch)
    saved = config.outdir / "discover_scan.checkpoint.json"
    assert json.loads(saved.read_text())["body"]["processed_reads"] == 2

    discover_toy_repeats(config)
    fresh = replace(config, outdir=tmp_path / "fresh", checkpoint_every=None)
    discover_toy_repeats(fresh)
    public = (
        "candidate_reads.tsv", "candidate_monomers.fa", "monomers.fa", "families.tsv",
        "family_similarity.tsv", "family_hierarchy.tsv", "family_audit_summary.json",
        "discovery_summary.json",
    )
    for filename in public:
        assert (config.outdir / filename).read_bytes() == (fresh.outdir / filename).read_bytes(), filename
    assert not saved.exists()


def test_changed_input_tail_and_corrupt_checkpoint_fail_closed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = config_for(tmp_path)
    interrupt_once(config, monkeypatch)
    reads = config.reads
    assert isinstance(reads, Path)
    original = reads.read_text(encoding="utf-8")
    reads.write_text(original[:-2] + "AA", encoding="utf-8")
    with pytest.raises(ValueError, match="input or configuration changed"):
        discover_toy_repeats(config)
    reads.write_text(original, encoding="utf-8")
    saved = config.outdir / "discover_scan.checkpoint.json"
    original_checkpoint = saved.read_text(encoding="utf-8")
    saved.write_text('{"body":', encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid discover scan checkpoint"):
        discover_toy_repeats(config)
    saved.write_text(original_checkpoint, encoding="utf-8")
    altered = json.loads(original_checkpoint)
    altered["body"]["prefix_sha256"] = "0" * 64
    altered["body_sha256"] = hashlib.sha256(checkpoint._encoded(altered["body"])).hexdigest()
    saved.write_text(json.dumps(altered), encoding="utf-8")
    with pytest.raises(ValueError, match="read prefix changed"):
        discover_toy_repeats(config)
    saved.write_text(original_checkpoint, encoding="utf-8")
    with (config.outdir / "candidate_reads.tsv").open("a", encoding="utf-8") as handle:
        handle.write("uncommitted tail\n")
    discover_toy_repeats(config)
    assert "uncommitted tail" not in (config.outdir / "candidate_reads.tsv").read_text()


def test_candidate_corruption_and_completed_output_fail_closed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = config_for(tmp_path)
    interrupt_once(config, monkeypatch)
    table = config.outdir / "candidate_reads.tsv"
    table.write_text(table.read_text().replace("TXC000001", "TXC999999"), encoding="utf-8")
    with pytest.raises(ValueError, match="candidate table changed"):
        discover_toy_repeats(config)
    config2 = replace(config, outdir=tmp_path / "another")
    interrupt_once(config2, monkeypatch)
    (config2.outdir / "discovery_summary.json").write_text("partial", encoding="utf-8")
    with pytest.raises(ValueError, match="coexists"):
        discover_toy_repeats(config2)


def test_checkpoint_rejects_unsupported_modes(tmp_path: Path) -> None:
    config = config_for(tmp_path)
    with pytest.raises(ValueError, match="one read file"):
        discover_toy_repeats(replace(config, reads=[config.reads, config.reads]))
    with pytest.raises(ValueError, match="sample-rate 1"):
        discover_toy_repeats(replace(config, sample_rate=0.5))


def test_failed_later_snapshot_preserves_previous_commit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = config_for(tmp_path)
    real_save = checkpoint.save
    calls = 0

    def fail_second(path: Path, body: dict[str, object]) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated snapshot I/O failure")
        real_save(path, body)

    with monkeypatch.context() as patch:
        patch.setattr(checkpoint, "save", fail_second)
        with pytest.raises(OSError, match="simulated snapshot"):
            discover_toy_repeats(config)
    saved = config.outdir / "discover_scan.checkpoint.json"
    assert json.loads(saved.read_text())["body"]["processed_reads"] == 2
    discover_toy_repeats(config)
    fresh = replace(config, outdir=tmp_path / "fresh", checkpoint_every=None)
    discover_toy_repeats(fresh)
    assert (config.outdir / "families.tsv").read_bytes() == (fresh.outdir / "families.tsv").read_bytes()
