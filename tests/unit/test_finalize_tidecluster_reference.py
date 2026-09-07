from pathlib import Path
import json

import pytest

from benchmarks.scripts.finalize_tidecluster_reference import validate_completed_stages


GNU_TIME = """\
User time (seconds): 1.0
System time (seconds): 0.1
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.2
Maximum resident set size (kbytes): 1024
Exit status: 0
"""


def make_stage_evidence(path: Path, complete: bool = True) -> None:
    profile = path / "profile"
    profile.mkdir(parents=True)
    (profile / "receipt.json").write_text(
        json.dumps(
            {
                "complete": complete,
                "requested_stage_count": 2,
                "completed_stage_count": 2,
            }
        )
    )
    (profile / "stages.tsv").write_text(
        "stage\texit_code\n" "tidehunter\t0\n" "clustering\t0\n"
    )
    for stage in ("tidehunter", "clustering"):
        (path / f"{stage}.gnu_time.txt").write_text(GNU_TIME)
    for name in (
        "tc_tidehunter.gff3",
        "tc_clustering.gff3_1.gff3",
        "tc_clustering.gff3",
        "tc_cmd_args.json",
        "tc_consensus/consensus_sequences_all.fasta",
    ):
        target = path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("fixture\n")


def test_validate_completed_external_stages(tmp_path: Path) -> None:
    make_stage_evidence(tmp_path)
    result = validate_completed_stages(tmp_path)
    assert result["profile"]["complete"] is True
    assert result["internal_gnu_time"]["clustering"]["wall_seconds"] == 1.2


def test_validate_rejects_incomplete_stage_receipt(tmp_path: Path) -> None:
    make_stage_evidence(tmp_path, complete=False)
    with pytest.raises(ValueError, match="not complete"):
        validate_completed_stages(tmp_path)
