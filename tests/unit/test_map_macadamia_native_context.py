"""The Macadamia source gate must keep record and molecule support distinct."""

import gzip

from benchmarks.scripts.map_macadamia_native_context import score


def test_three_sra_records_do_not_prove_three_zmws(tmp_path):
    config = {
        "context_name": "ctg.000105F:115439-124594",
        "array_start_0_in_context": 3000,
        "array_end_0_in_context": 6155,
        "minimum_natural_flank_bp_each_side": 1000,
        "minimum_nmatch_over_alignment_columns": 0.99,
        "minimum_distinct_original_records": 3,
    }
    target = config["context_name"]

    def paf_row(name: str, matches: int, tag: str = "tp:A:P") -> str:
        return f"{name}\t10000\t500\t8000\t+\t{target}\t9155\t1900\t7200\t{matches}\t5300\t60\t{tag}\n"

    path = tmp_path / "mappings.paf.gz"
    with gzip.open(path, "wt") as handle:
        for identifier in ("SRR13557763.1", "SRR13557763.2", "SRR13557762.1"):
            handle.write(paf_row(identifier, 5290))
        handle.write(paf_row("SRR13557763.low_identity", 5000))
        handle.write(paf_row("SRR13557763.secondary", 5290, "tp:A:S"))
    result = score(path, config, tmp_path / "score")
    assert result["qualifying_distinct_record_count"] == 3
    assert result["record_gate_pass"] is True
    assert result["distinct_zmw_count"] is None
    assert result["molecule_gate_pass"] is False
    assert result["primary_geometric_spanner_alignment_count"] == 4
