"""Guard the compact equal-input competitor replay against mismatched artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "benchmarks/evidence/common_assembly_5700_v1"
INPUT = ROOT / "benchmarks/evidence/competitor_local_smoke_20260916/inputs/cendetecthor_pipeline.fa"


def test_common_assembly_native_outputs_and_resource_states() -> None:
    assessment = json.loads((ARCHIVE / "assessment.json").read_text())
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    assert digest == assessment["input_sha256"]
    assert assessment["engineered_truth_array_0_based_half_open"] == [0, 5700]
    for name in ("cendetecthor", "tidehunter", "tidecluster"):
        manifest = json.loads((ARCHIVE / "manifests" / f"{name}_rss.json").read_text())
        command = manifest["stages"][0]["command"]
        assert any("toy.chr1.fasta" in argument for argument in command)
    with (ARCHIVE / "profiles/combined_rss/stages.tsv").open() as handle:
        stages = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["stage"] for row in stages] == [
        "cendetecthor_full_pipeline", "tidehunter", "tidecluster_clustering"
    ]
    assert all(row["exit_code"] == "0" for row in stages)
    assert float(stages[0]["peak_process_tree_rss_mib"]) > 0
    assert float(stages[2]["peak_process_tree_rss_mib"]) > 0
    assert (ARCHIVE / "native/cendetecthor_selected.bed").read_text().strip() == "toy\t0\t5700\t60"
    gff = (ARCHIVE / "native/tidecluster.gff3").read_text().splitlines()[1].split("\t")
    assert gff[:5] == ["toy:0-5700", "TideCluster", "tandem_repeat", "1", "5700"]
    assert assessment["cendetecthor_selected_interval_0_based_half_open"] == [0, 5700]
    assert assessment["tidecluster_selected_interval_0_based_half_open"] == [0, 5700]


def test_common_assembly_archive_hash_manifest() -> None:
    for line in (ARCHIVE / "SHA256SUMS.txt").read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = ARCHIVE / relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
