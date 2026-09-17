"""Guard the archived technical N/A against accidental accuracy promotion."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "benchmarks/evidence/cendetecthor_171_precondition_20260917"
INPUT = ROOT / "benchmarks/assembly_hor_comparator_v1/control/control.chr1.fasta"


def test_fixed_input_receipts_and_native_failure_evidence() -> None:
    protocol = json.loads((ARCHIVE / "protocol.json").read_text())
    assessment = json.loads((ARCHIVE / "assessment.json").read_text())
    assert hashlib.sha256(INPUT.read_bytes()).hexdigest() == protocol["input_sha256"]
    assert assessment["status"] == "technical_incomplete"
    assert assessment["accuracy"] is None
    assert [item["name"] for item in protocol["fixed_conditions_in_order"]] == [
        "default_replay", "documented_k8_w5000"
    ]

    for name, expected_error in (
        ("default_replay_corrected", "IndexError: list index out of range"),
        ("documented_k8_w5000_corrected", "NameError: name 'cons' is not defined"),
    ):
        manifest = json.loads((ARCHIVE / f"{name}.manifest.json").read_text())
        command = manifest["stages"][0]["command"]
        assert "consFile=no" in command
        assert not any(token.startswith("expMonSize=") for token in command)
        assert assessment["conditions"][name]["input_sha256"] == protocol["input_sha256"]
        assert assessment["conditions"][name]["exit_code"] != 0
        assert not assessment["conditions"][name]["native_final_decomposition_exists"]
        assert not assessment["conditions"][name]["native_hor_tree_exists"]
        with (ARCHIVE / name / "stages.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1 and int(rows[0]["exit_code"]) != 0
        assert float(rows[0]["peak_process_tree_rss_mib"]) > 0
        stderr = (ARCHIVE / name / rows[0]["stderr_log"]).read_text()
        assert expected_error in stderr

    default_bed = ARCHIVE / "default_replay_corrected/native/results/wind2analize/control.chr1.windows.bed"
    diagnostic_bed = ARCHIVE / "documented_k8_w5000_corrected/native/results/wind2analize/control.FULLchr.windows.filtered.bed"
    consensus = ARCHIVE / "documented_k8_w5000_corrected/native/monomers/control_chr0cons.fasta"
    assert default_bed.stat().st_size == 0
    assert diagnostic_bed.read_text().strip() == "control\t0\t36371\t1710"
    assert consensus.stat().st_size == 0

    for line in (ARCHIVE / "SHA256SUMS.txt").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        assert hashlib.sha256((ARCHIVE / relative).read_bytes()).hexdigest() == digest
