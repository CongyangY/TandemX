from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from tandemx.quantify import mvp
from tandemx.quantify.mvp import QuantifyConfig, quantify_toy_copy_number


def fixture_config(tmp_path: Path, *, fastq: bool = False) -> QuantifyConfig:
    reads = tmp_path / ("reads.fq" if fastq else "reads.fa")
    sequence = "ACGTTCAGGACACGTTCAGGAC"
    if fastq:
        reads.write_text(
            "".join(f"@r{i}\n{sequence}\n+\n{'I' * len(sequence)}\n" for i in range(8)),
            encoding="utf-8",
        )
    else:
        reads.write_text(
            "".join(f">r{i}\n{sequence}\n" for i in range(8)), encoding="utf-8"
        )
    monomers = tmp_path / "monomers.fa"
    monomers.write_text(
        ">family_id=TXF000001;length_bp=11\nACGTTCAGGAC\n", encoding="utf-8"
    )
    return QuantifyConfig(
        reads=reads, monomers=monomers, genome_size=100,
        outdir=tmp_path / "resumed", k=5, haploid_depth=1.0,
        checkpoint_every=2,
    )


@pytest.mark.parametrize("fastq", [False, True])
@pytest.mark.parametrize("backend", ["python", "rust"])
def test_interrupted_scan_resumes_to_identical_public_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fastq: bool, backend: str,
) -> None:
    config = replace(fixture_config(tmp_path, fastq=fastq), kmer_backend=backend)
    original = mvp._save_scan_checkpoint
    saved = 0

    def interrupt_after_second_checkpoint(path: Path, body: dict[str, object]) -> None:
        nonlocal saved
        original(path, body)
        saved += 1
        if saved == 2:
            raise KeyboardInterrupt()

    with monkeypatch.context() as patch:
        patch.setattr(mvp, "_save_scan_checkpoint", interrupt_after_second_checkpoint)
        with pytest.raises(KeyboardInterrupt):
            quantify_toy_copy_number(config)
    checkpoint = config.outdir / "quantify_scan.checkpoint.json"
    assert checkpoint.is_file()
    assert not (config.outdir / "copy_number.tsv").exists()
    assert json.loads(checkpoint.read_text())["body"]["read_count"] == 4

    quantify_toy_copy_number(config)
    fresh = replace(config, outdir=tmp_path / "fresh", checkpoint_every=None)
    quantify_toy_copy_number(fresh)
    assert (config.outdir / "copy_number.tsv").read_bytes() == (
        fresh.outdir / "copy_number.tsv"
    ).read_bytes()
    assert not checkpoint.exists()


def test_corrupt_checkpoint_fails_closed(tmp_path: Path) -> None:
    config = fixture_config(tmp_path)
    config.outdir.mkdir()
    checkpoint = config.outdir / "quantify_scan.checkpoint.json"
    checkpoint.write_text('{"body":', encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid quantify scan checkpoint"):
        quantify_toy_copy_number(config)
    assert not (config.outdir / "copy_number.tsv").exists()


def test_changed_input_tail_and_stale_output_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = fixture_config(tmp_path)
    original = mvp._save_scan_checkpoint

    def interrupt(path: Path, body: dict[str, object]) -> None:
        original(path, body)
        raise KeyboardInterrupt()

    with monkeypatch.context() as patch:
        patch.setattr(mvp, "_save_scan_checkpoint", interrupt)
        with pytest.raises(KeyboardInterrupt):
            quantify_toy_copy_number(config)
    reads = config.reads
    assert isinstance(reads, Path)
    text = reads.read_text(encoding="utf-8")
    reads.write_text(text[:-2] + "AA", encoding="utf-8")
    with pytest.raises(ValueError, match="input or scientific configuration changed"):
        quantify_toy_copy_number(config)
    reads.write_text(text, encoding="utf-8")
    checkpoint = config.outdir / "quantify_scan.checkpoint.json"
    original_checkpoint = checkpoint.read_text(encoding="utf-8")
    tampered = json.loads(original_checkpoint)
    tampered["body"]["counts"]["ACGTT"] = 999999
    checkpoint.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        quantify_toy_copy_number(config)
    checkpoint.write_text(original_checkpoint, encoding="utf-8")
    (config.outdir / "copy_number.tsv").write_text("partial\n", encoding="utf-8")
    with pytest.raises(ValueError, match="coexists with copy_number.tsv"):
        quantify_toy_copy_number(config)
    assert (config.outdir / "copy_number.tsv").read_text() == "partial\n"
