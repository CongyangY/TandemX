"""Preserve the exact source FASTQ records that passed the frozen context gate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.scripts.map_macadamia_native_context import sha256
from benchmarks.scripts.prepare_native_pair_audit import extract_selected


def extract(protocol: Path, score_dir: Path, output_dir: Path) -> dict:
    config = json.loads(protocol.read_text())
    mapping = json.loads((score_dir / "run_receipt.json").read_text())
    if mapping["status"] != "completed" or not mapping["record_gate_pass"]:
        raise ValueError("frozen context record gate did not pass")
    if mapping["protocol_sha256"] != sha256(protocol):
        raise ValueError("mapping protocol changed")
    table = score_dir / "qualifying_record_alignments.tsv"
    if mapping["qualifying_tsv_sha256"] != sha256(table):
        raise ValueError("qualifying record table changed")
    with table.open(newline="") as handle:
        selected = {row["read_id"] for row in csv.DictReader(handle, delimiter="\t")}
    if len(selected) != mapping["qualifying_distinct_record_count"]:
        raise ValueError("qualified read count differs from frozen score")
    runs = {item["run"]: item for item in config["read_inputs"]}
    if any(read_id.split(".", 1)[0] not in runs for read_id in selected):
        raise ValueError("qualified ID has an unexpected run prefix")
    output_dir.mkdir(parents=True, exist_ok=True)
    receipts = []
    for run, item in runs.items():
        ids = {read_id for read_id in selected if read_id.startswith(run + ".")}
        if not ids:
            continue
        receipts.append(extract_selected(Path(item["path"]), item["sha256"], ids,
                                         output_dir / f"{run}_qualified.fastq.gz"))
    result = {
        "status": "complete",
        "mapping_protocol_sha256": sha256(protocol),
        "mapping_receipt_sha256": sha256(score_dir / "run_receipt.json"),
        "qualifying_record_count": len(selected),
        "distinct_zmw_verified": False,
        "runs": receipts,
    }
    (output_dir / "extraction_receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(extract(args.protocol, args.score_dir, args.output_dir), sort_keys=True))


if __name__ == "__main__":
    main()
