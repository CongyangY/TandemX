from benchmarks.scripts.analyze_external_tool_comparison import (
    build_macro_rows,
    build_pairwise_rows,
    build_chart_rows,
    build_dataset_table_rows,
    validate_inputs,
)


def summary_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for tool, runtime, memory, precision in (
        ("tandemx", 2.0, 40.0, 1.0),
        ("trf", 4.0, 20.0, 1.0),
        ("tidehunter", 6.0, 100.0, 0.99),
    ):
        rows.append(
            {
                "dataset_id": "long_noisy_421",
                "tool": tool,
                "total_bases": "1000000",
                "monomer_length": "171",
                "error_rate": "0.05",
                "median_runtime_seconds": str(runtime),
                "median_peak_memory_mb": str(memory),
                "runtime_cv": "0.02",
                "recall": "1.0",
                "precision": str(precision),
                "false_positive_rate": "0.0",
                "monomer_length_mae_bp": "0.0",
                "successful_runs": "3",
                "deterministic_predictions": "true",
            }
        )
    return rows


def manifest_rows() -> list[dict[str, str]]:
    return [
        {
            "dataset_id": "long_noisy_421",
            "source_mode": "reused_existing_files",
            "truth_id_coverage": "1.0",
            "observed_read_count": "100",
            "truth_row_count": "100",
            "expected_read_length": "10000",
            "positive_reads": "75",
        }
    ]


def test_validate_inputs_accepts_complete_deterministic_matrix() -> None:
    validation = validate_inputs(summary_rows(), manifest_rows())
    assert validation["status"] == "passed"
    assert validation["successful_runs"] == 9


def test_validate_inputs_rejects_empty_inputs() -> None:
    validation = validate_inputs([], [])
    assert validation["status"] == "failed"
    assert "summary_is_empty" in validation["issues"]
    assert "dataset_manifest_is_empty" in validation["issues"]


def test_build_pairwise_rows_calculates_speed_and_memory_tradeoff() -> None:
    rows = build_pairwise_rows(summary_rows())
    trf = next(row for row in rows if row["competitor"] == "trf")
    tidehunter = next(row for row in rows if row["competitor"] == "tidehunter")
    assert trf["tandemx_speedup"] == "2.000000"
    assert trf["memory_winner"] == "trf"
    assert tidehunter["tandemx_speedup"] == "3.000000"
    assert tidehunter["tandemx_memory_reduction_fraction"] == "0.600000"


def test_build_macro_rows_preserves_accuracy_difference() -> None:
    rows = {row["tool"]: row for row in build_macro_rows(summary_rows())}
    assert rows["tandemx"]["macro_precision"] == "1.000000"
    assert rows["tidehunter"]["macro_precision"] == "0.990000"


def test_build_chart_rows_executes_sql_selection() -> None:
    rows = build_chart_rows(summary_rows())
    assert [row["tool"] for row in rows] == ["TandemX", "TRF", "TideHunter"]
    assert rows[0]["runtime_seconds"] == 2.0
    assert rows[-1]["peak_memory_mb"] == 100.0


def test_build_dataset_table_rows_executes_sql_selection() -> None:
    rows = build_dataset_table_rows(manifest_rows())
    assert rows == [
        {
            "dataset_id": "long_noisy_421",
            "observed_read_count": 100,
            "expected_read_length": 10000,
            "positive_reads": 75,
            "truth_id_coverage": 1.0,
            "source_mode": "reused_existing_files",
        }
    ]
