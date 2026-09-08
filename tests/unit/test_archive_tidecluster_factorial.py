from benchmarks.scripts.archive_tidecluster_factorial import aggregate_settings


def test_aggregate_settings_preserves_seed_range_and_resource_values() -> None:
    rows = [
        {
            "seed": "1",
            "setting": "default",
            "array_recall": "0.8",
            "array_precision": "0.5",
            "base_union_recall": "0.7",
            "base_union_precision": "0.6",
            "matched_boundary_mae_bp": "2",
            "matched_period_mae_bp": "1",
            "cyclic_monomer_recall": "0.4",
            "homologous_consensus_fraction": "0.3",
            "tidehunter_wall_seconds": "10",
            "tidehunter_maximum_rss_kb": "100",
            "clustering_wall_seconds": "20",
            "clustering_maximum_rss_kb": "200",
        },
        {
            "seed": "2",
            "setting": "default",
            "array_recall": "1.0",
            "array_precision": "0.7",
            "base_union_recall": "0.9",
            "base_union_precision": "0.8",
            "matched_boundary_mae_bp": "4",
            "matched_period_mae_bp": "3",
            "cyclic_monomer_recall": "0.6",
            "homologous_consensus_fraction": "0.5",
            "tidehunter_wall_seconds": "12",
            "tidehunter_maximum_rss_kb": "120",
            "clustering_wall_seconds": "22",
            "clustering_maximum_rss_kb": "220",
        },
    ]
    result = aggregate_settings(rows)["default"]
    assert result["run_count"] == 2
    assert result["seeds"] == [1, 2]
    assert result["metrics"]["array_recall"] == {
        "mean": 0.9,
        "minimum": 0.8,
        "maximum": 1.0,
    }
    assert result["metrics"]["clustering_maximum_rss_kb"]["mean"] == 210
