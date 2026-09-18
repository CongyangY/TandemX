"""Score selected Macadamia source records against the complete assembly."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.scripts.map_macadamia_native_context import sha256
from benchmarks.scripts.score_native_genomewide_mapping import parse_paf_line


def score(protocol: Path, paf: Path, outdir: Path) -> dict:
    config = json.loads(protocol.read_text())
    if config.get("status") != "frozen_before_full_reference_mapping":
        raise ValueError("full-reference mapping protocol is not frozen")
    archive = json.loads(Path(config["extraction_receipt"]).read_text())
    if archive["status"] != "complete":
        raise ValueError("selected source reads are not complete")
    expected = {}
    for run in archive["runs"]:
        fastq = Path(run["selected_fastq"])
        run_id = fastq.name.split("_", 1)[0]
        if ("selected_fastq_sha256" in config and
                config["selected_fastq_sha256"].get(run_id) != run["selected_fastq_sha256"]):
            raise ValueError("selected FASTQ changed after full-reference freeze")
        if sha256(fastq) != run["selected_fastq_sha256"]:
            raise ValueError("selected source FASTQ changed")
        for item in run["reads"]:
            if item["id"] in expected:
                raise ValueError("duplicate selected read ID")
            expected[item["id"]] = []
    if len(expected) != config["expected_selected_record_count"]:
        raise ValueError("selected read count changed")
    with paf.open() as handle:
        for line in handle:
            row = parse_paf_line(line)
            if row["read_id"] not in expected:
                raise ValueError("full-reference PAF has an undeclared read")
            expected[row["read_id"]].append(row)
    rows = []
    for read_id, hits in sorted(expected.items()):
        primaries = [hit for hit in hits if hit["primary"]]
        if len(primaries) > 1:
            raise ValueError(f"multiple reported primaries: {read_id}")
        if not primaries:
            rows.append({"read_id": read_id, "status": "no_primary", "target": "",
                         "target_start_0": "", "target_end_0": "", "identity": "",
                         "mapq": "", "left_flank_bp": "", "right_flank_bp": "",
                         "reported_secondary_count": len(hits)})
            continue
        hit = primaries[0]
        left = config["array_start_0"] - hit["target_start0"]
        right = hit["target_end0"] - config["array_end_0"]
        passed = (hit["target"] == config["contig"] and
                  left >= config["minimum_natural_flank_bp_each_side"] and
                  right >= config["minimum_natural_flank_bp_each_side"] and
                  hit["identity"] >= config["minimum_identity"] and
                  hit["mapq"] >= config["minimum_mapq"])
        rows.append({"read_id": read_id,
                     "status": "full_reference_primary_spanner" if passed else "not_eligible",
                     "target": hit["target"], "target_start_0": hit["target_start0"],
                     "target_end_0": hit["target_end0"], "identity": f"{hit['identity']:.9f}",
                     "mapq": hit["mapq"], "left_flank_bp": left,
                     "right_flank_bp": right,
                     "reported_secondary_count": len(hits) - 1})
    outdir.mkdir(parents=True, exist_ok=False)
    with (outdir / "per_read.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "selected_record_count": len(rows),
        "full_reference_primary_spanner_count": sum(
            row["status"] == "full_reference_primary_spanner" for row in rows),
        "no_primary_count": sum(row["status"] == "no_primary" for row in rows),
        "total_reported_secondary_count": sum(row["reported_secondary_count"] for row in rows),
        "distinct_zmw_verified": False,
        "source_pairing": "paper_same_sample_exact_extraction_unverified",
        "physical_copy_truth": "unavailable",
        "protocol_sha256": sha256(protocol),
        "paf_sha256": sha256(paf),
        "per_read_sha256": sha256(outdir / "per_read.tsv"),
    }
    (outdir / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(score(args.protocol, args.paf, args.outdir), sort_keys=True))


if __name__ == "__main__":
    main()
