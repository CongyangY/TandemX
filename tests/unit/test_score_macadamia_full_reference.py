"""Full-reference source audit reports unmapped and wrong-locus reads."""

import gzip
import json

from benchmarks.scripts.map_macadamia_native_context import sha256
from benchmarks.scripts.score_macadamia_full_reference import score


def test_full_reference_denominator_includes_unmapped_and_wrong_locus(tmp_path):
    reads = tmp_path / "selected.fastq.gz"
    with gzip.open(reads, "wb") as handle:
        for name in ("SRR.1", "SRR.2", "SRR.3"):
            handle.write(f"@{name}\nACGT\n+\nIIII\n".encode())
    extraction = tmp_path / "extraction.json"
    extraction.write_text(json.dumps({"status": "complete", "runs": [{
        "selected_fastq": str(reads), "selected_fastq_sha256": sha256(reads),
        "reads": [{"id": name} for name in ("SRR.1", "SRR.2", "SRR.3")],
    }]}))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "status": "frozen_before_full_reference_mapping",
        "extraction_receipt": str(extraction), "expected_selected_record_count": 3,
        "contig": "ctg", "array_start_0": 3000, "array_end_0": 6000,
        "minimum_natural_flank_bp_each_side": 1000,
        "minimum_identity": 0.99, "minimum_mapq": 20,
    }))
    paf = tmp_path / "hits.paf"
    paf.write_text(
        "SRR.1\t10000\t0\t9000\t+\tctg\t20000\t1000\t8000\t6990\t7000\t60\ttp:A:P\n"
        "SRR.2\t10000\t0\t9000\t+\tother\t20000\t1000\t8000\t6990\t7000\t60\ttp:A:P\n"
    )
    result = score(config, paf, tmp_path / "out")
    assert result["selected_record_count"] == 3
    assert result["full_reference_primary_spanner_count"] == 1
    assert result["no_primary_count"] == 1
    table = (tmp_path / "out/per_read.tsv").read_text()
    assert "SRR.2\tnot_eligible" in table
    assert "SRR.3\tno_primary" in table
