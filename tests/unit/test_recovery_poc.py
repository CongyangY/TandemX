from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from tandemx.recovery.poc import FROZEN_RULES, evaluate, recruit, sha256, validate_rules


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _paf(query: str, query_start: int, query_end: int, target: str, target_length: int,
         target_start: int, target_end: int, *, strand: str = "+", matches: int | None = None) -> str:
    block = query_end - query_start
    return "\t".join(map(str, (
        query, 10_000, query_start, query_end, strand, target, target_length,
        target_start, target_end, block if matches is None else matches, block, 60,
    ))) + "\n"


def _state(tmp_path: Path, *, anchors: list[dict]) -> tuple[dict, Path]:
    reads = _write(tmp_path / "reads.fa", "".join(f">read{i}\n{'A' * 10000}\n" for i in range(1, 4)))
    copy_number = _write(tmp_path / "copy_number.tsv", "family_id\testimated_bp\nF1\t12000\n")
    state = {
        "config": {"families": ["F1"], "inputs": {"reads": str(reads), "copy_number": str(copy_number)}},
        "loci": [{"locus_id": "L0001", "family_id": "F1", "chromosome": "chr1", "start": 100, "end": 900}],
        "anchors": anchors,
        "catalog_lengths": {"F1": 10},
    }
    return state, reads


def _anchors() -> list[dict]:
    return [
        {"anchor_id": "L0001_left_0", "locus_id": "L0001", "family_id": "F1", "side": "left", "offset": 0, "start": 0, "end": 2000, "eligible": True},
        {"anchor_id": "L0001_right_0", "locus_id": "L0001", "family_id": "F1", "side": "right", "offset": 0, "start": 3000, "end": 5000, "eligible": True},
    ]


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_evaluate_retains_negative_no_anchor_outputs(tmp_path: Path) -> None:
    state, _reads = _state(tmp_path, anchors=[])
    paf = _write(tmp_path / "reads_to_targets.paf", "".join(
        _paf(f"read{i}", 1000, 1600, "repeat_F1", 30_000, 0, 600) for i in range(1, 4)
    ))

    evaluate(tmp_path, state, paf)

    row = _rows(tmp_path / "recovery_candidates.tsv")[0]
    assert row["recovery_status"] == "unresolved_no_unique_anchor"
    assert row["recovered_bp"] == "NA"
    assert (tmp_path / "recovered_sequences.fasta").read_text(encoding="utf-8") == ""
    assert json.loads((tmp_path / "candidate_lock.json").read_text(encoding="utf-8"))["expansion_decision"] == "stop_recovery_expansion_negative_poc"
    assert len(_rows(tmp_path / "recruited_reads.tsv")) == 3


def test_evaluate_writes_only_partially_resolved_observed_span(tmp_path: Path) -> None:
    state, _reads = _state(tmp_path, anchors=_anchors())
    rows: list[str] = []
    for index in range(1, 4):
        read = f"read{index}"
        rows.extend((
            _paf(read, 1_000, 3_000, "L0001_left_0", 2_000, 0, 2_000),
            _paf(read, 7_000, 9_000, "L0001_right_0", 2_000, 0, 2_000),
            _paf(read, 4_000, 4_600, "repeat_F1", 30_000, 0, 600),
        ))
    paf = _write(tmp_path / "reads_to_targets.paf", "".join(rows))

    evaluate(tmp_path, state, paf)

    row = _rows(tmp_path / "recovery_candidates.tsv")[0]
    assert row["recovery_status"] == "partially_resolved"
    assert row["dual_flank_read_count"] == "3"
    assert row["recovered_bp"] == "NA"
    assert row["candidate_span_bp"] == "4000"
    assert row["maximum_read_span_bp"] == "4000"
    assert row["maximum_recruited_read_length_bp"] == "10000"
    fasta = (tmp_path / "recovered_sequences.fasta").read_text(encoding="utf-8")
    assert "status=partially_resolved;unpolished=true" in fasta
    assert fasta.count("\n" + "A" * 4000 + "\n") == 1


def test_evaluate_rejects_unexpected_target_without_silent_output(tmp_path: Path) -> None:
    state, _reads = _state(tmp_path, anchors=[])
    paf = _write(tmp_path / "reads_to_targets.paf", _paf("read1", 1000, 1600, "newer_proxy", 30_000, 0, 600))

    with pytest.raises(ValueError, match="Unexpected recruitment target: newer_proxy"):
        evaluate(tmp_path, state, paf)
    assert not (tmp_path / "recovery_candidates.tsv").exists()


def test_frozen_rules_reject_any_enrollment_change() -> None:
    validate_rules({"rules": FROZEN_RULES})
    changed = dict(FROZEN_RULES)
    changed["flank_bp"] = 2001
    with pytest.raises(ValueError, match="frozen implemented contract"):
        validate_rules({"rules": changed})


def test_evaluate_rejects_selected_read_absent_from_source_without_lock(tmp_path: Path) -> None:
    state, _reads = _state(tmp_path, anchors=_anchors())
    paf_rows: list[str] = []
    for read in ("read2", "read3", "absent_read"):
        paf_rows.extend((
            _paf(read, 1_000, 3_000, "L0001_left_0", 2_000, 0, 2_000),
            _paf(read, 7_000, 9_000, "L0001_right_0", 2_000, 0, 2_000),
            _paf(read, 4_000, 4_600, "repeat_F1", 30_000, 0, 600),
        ))
    paf = _write(tmp_path / "reads_to_targets.paf", "".join(paf_rows))

    with pytest.raises(ValueError, match="source read absent"):
        evaluate(tmp_path, state, paf)
    assert not (tmp_path / "candidate_lock.json").exists()


def test_evaluate_rejects_selected_family_absent_from_abundance_table(tmp_path: Path) -> None:
    state, _reads = _state(tmp_path, anchors=[])
    _write(Path(state["config"]["inputs"]["copy_number"]), "family_id\testimated_bp\n")
    paf = _write(tmp_path / "reads_to_targets.paf", _paf("read1", 1000, 1600, "repeat_F1", 30_000, 0, 600))

    with pytest.raises(ValueError, match="absent from the frozen abundance table"):
        evaluate(tmp_path, state, paf)
    assert not (tmp_path / "recovery_candidates.tsv").exists()


def _locked_recruitment_dir(path: Path) -> tuple[Path, dict]:
    """Create the minimal hash-locked state needed to reach mapping import checks."""
    path.mkdir()
    reads = _write(path.parent / "enrolled_reads.fa", ">r1\nAAAA\n")
    abundance = _write(path.parent / "enrolled_copy_number.tsv", "family_id\testimated_bp\nF1\t4\n")
    minimap = _write(path.parent / "pinned_minimap2", "pinned executable bytes\n")
    targets = _write(path / "recruitment_targets.fa", ">repeat_F1\nAAAA\n")
    audit = _write(path / "flank_audit.tsv", "anchor_id\n")
    inputs = {}
    for key, source in (("reads", reads), ("copy_number", abundance)):
        inputs[key] = {
            "path": str(source.resolve()),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
        }
    config = {
        "families": ["F1"],
        "inputs": {key: entry["path"] for key, entry in inputs.items()},
        "minimap2": str(minimap.resolve()),
        "rules": dict(FROZEN_RULES),
    }
    state = {
        "config": config, "loci": [], "anchors": [], "catalog_lengths": {"F1": 4},
        "targets_sha256": sha256(targets), "anchor_audit_sha256": sha256(audit),
    }
    manifest = {"inputs": inputs, "minimap2_sha256": sha256(minimap)}
    (path / "prepared.json").write_text(json.dumps(state), encoding="utf-8")
    (path / "input_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (path / "prepare_lock.json").write_text(json.dumps({
        "prepared_sha256": sha256(path / "prepared.json"),
        "input_manifest_sha256": sha256(path / "input_manifest.json"),
    }), encoding="utf-8")
    return path, manifest


@pytest.mark.parametrize("kind", ("receipt_command", "receipt_hash"))
def test_recruit_rejects_non_equivalent_mapping_receipt_before_import(tmp_path: Path, kind: str) -> None:
    current, manifest = _locked_recruitment_dir(tmp_path / "current")
    previous, _ = _locked_recruitment_dir(tmp_path / "previous")
    paf = _write(previous / "reads_to_targets.paf", "")
    command = [
        str((tmp_path / "pinned_minimap2").resolve()), "-x", "map-hifi", "-c", "-N", "50",
        "-p", "0.5", "--secondary=yes", "-K", "50M", "-t", "4",
        str(previous / "recruitment_targets.fa"), str((tmp_path / "enrolled_reads.fa").resolve()),
    ]
    receipt = {"exit_code": 0, "command": command, "sha256": sha256(paf)}
    if kind == "receipt_command":
        receipt["command"] = ["tampered"]
    else:
        receipt["sha256"] = "0" * 64
    (previous / "reads_to_targets.receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
    _write(previous / "reads_to_targets.log", "prior mapper output\n")
    assert json.loads((previous / "input_manifest.json").read_text(encoding="utf-8")) == manifest

    with pytest.raises(ValueError, match="not exact-input/command equivalent"):
        recruit(current, reuse_mapping_from=previous)

    assert not (current / "reads_to_targets.paf").exists()
    assert not (current / "mapping_import.json").exists()


def test_only_bounded_resource_rules_may_change() -> None:
    allowed = dict(FROZEN_RULES)
    allowed.update(maximum_recruited_ids=250_000, maximum_retained_alignment_rows=5_000_000)
    validate_rules({"rules": allowed})

    changed_scientific = dict(allowed)
    changed_scientific["repeat_read_identity_min"] = 0.86
    with pytest.raises(ValueError, match="frozen implemented contract"):
        validate_rules({"rules": changed_scientific})

    changed_mapping = dict(allowed)
    changed_mapping["mapping_timeout_seconds"] = 3600
    with pytest.raises(ValueError, match="frozen implemented contract"):
        validate_rules({"rules": changed_mapping})
