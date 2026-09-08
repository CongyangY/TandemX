import csv
import json
import math
from pathlib import Path

from benchmarks.scripts.independent_verify_tidecluster_factorial import (
    base_union_metrics,
    canonical_monomer,
    maximum_matching,
    reaches_cyclic_threshold,
    same_value,
    sha256,
    verify,
)
from benchmarks.tidecluster.factorial_continue import (
    CELL_FATE_FIELDS,
    STAGE_FATE_FIELDS,
    blank_summary,
)
from benchmarks.tidecluster.factorial_run import SUMMARY_FIELDS


def test_independent_interval_union_and_matching() -> None:
    predicted = [
        {"chrom": "chr1", "start": 10, "end": 30},
        {"chrom": "chr1", "start": 20, "end": 40},
    ]
    truth = [{"chrom": "chr1", "start": 15, "end": 35}]
    recall, precision = base_union_metrics(predicted, truth)
    assert recall == 1.0
    assert precision == 20 / 30
    assert maximum_matching([[0], [0]]) == {0: 0}


def test_independent_cyclic_threshold_handles_rotation_and_reverse_complement() -> None:
    truth = "AACCGT"
    assert canonical_monomer(truth) == canonical_monomer("CGTAAC")
    assert reaches_cyclic_threshold(truth, "CGTAAC")
    assert reaches_cyclic_threshold(truth, "ACGGTT")
    assert not reaches_cyclic_threshold(truth, "TTTTTT")


def test_same_value_handles_json_null_as_nan() -> None:
    assert same_value(None, math.nan)
    assert same_value("0.5", 0.5)
    assert not same_value("0.6", 0.5)


def _write_tsv(path: Path, fields: tuple[str, ...] | list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_verifier_accepts_explicit_unavailable_accuracy_without_zero(
    tmp_path: Path,
) -> None:
    genome = tmp_path / "genome"
    genome.mkdir()
    (genome / "genome.fa").write_text(">chr\nACGT\n")
    (genome / "catalogue.fa").write_text(">f1\nACGT\n")
    (genome / "truth_copy_number.tsv").write_text("chrom\tstart\tend\tperiod\n")
    config = {
        "datasets": [
            {
                "seed": 1,
                "genome_dir": str(genome),
                "genome_sha256": sha256(genome / "genome.fa"),
                "catalogue_sha256": sha256(genome / "catalogue.fa"),
                "truth_sha256": sha256(genome / "truth_copy_number.tsv"),
            }
        ],
        "settings": {"default_primary": {}},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    row = blank_summary(
        1,
        "default_primary",
        "external_resource_failure",
        "accuracy_unavailable_not_zero",
        {"wall_seconds": 2.0, "maximum_rss_kb": 100},
    )
    summary_path = tmp_path / "summary.tsv"
    _write_tsv(summary_path, SUMMARY_FIELDS, [row])
    cell_fates = tmp_path / "cell_fates.tsv"
    _write_tsv(
        cell_fates,
        CELL_FATE_FIELDS,
        [
            {
                "seed": 1,
                "setting": "default_primary",
                "status": "external_resource_failure",
                "tidehunter_status": "failed_exit_137",
                "clustering_status": "not_started_dependency_failure",
                "evaluation_status": "not_started_dependency_failure",
                "source": "test",
                "reason": "resource_failure",
            }
        ],
    )
    stage_fates = tmp_path / "stage_fates.tsv"
    _write_tsv(
        stage_fates,
        STAGE_FATE_FIELDS,
        [
            {
                "stage": "s1_default_primary_tidehunter",
                "seed": 1,
                "setting": "default_primary",
                "component": "tidehunter",
                "status": "failed",
                "execution_source": "test",
                "exit_code": 1,
                "reason": "resource_failure",
            },
            {
                "stage": "s1_default_primary_clustering",
                "seed": 1,
                "setting": "default_primary",
                "component": "clustering",
                "status": "not_started_dependency_failure",
                "execution_source": "test",
                "exit_code": "",
                "reason": "tidehunter_failed",
            },
        ],
    )
    (tmp_path / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": True,
                "run_count": 1,
                "successful_run_count": 0,
                "accuracy_complete": False,
                "summary_sha256": sha256(summary_path),
                "cell_fates_sha256": sha256(cell_fates),
                "stage_fates_sha256": sha256(stage_fates),
            }
        )
    )
    payload = verify(config_path, tmp_path, tmp_path / "verification.json")
    assert payload["verification_passed"] is True
    assert payload["unavailable_accuracy_run_count"] == 1
