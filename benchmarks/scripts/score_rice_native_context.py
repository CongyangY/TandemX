"""Apply the frozen natural-flank spanner rule to rice v2 context PAF."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def parse_paf(line: str) -> dict:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("PAF line has fewer than 12 fields")
    tags = {item.split(":", 2)[0]: item for item in fields[12:] if ":" in item}
    row = {
        "read_id": fields[0], "read_length": int(fields[1]),
        "query_start_0": int(fields[2]), "query_end_0": int(fields[3]),
        "strand": fields[4], "target": fields[5], "target_length": int(fields[6]),
        "target_start_0": int(fields[7]), "target_end_0": int(fields[8]),
        "matching_bases": int(fields[9]), "alignment_columns": int(fields[10]),
        "mapq": int(fields[11]), "primary": tags.get("tp") == "tp:A:P",
    }
    if not (0 <= row["query_start_0"] < row["query_end_0"] <= row["read_length"]):
        raise ValueError("invalid query coordinates")
    if not (0 <= row["target_start_0"] < row["target_end_0"] <= row["target_length"]):
        raise ValueError("invalid target coordinates")
    if not (0 < row["matching_bases"] <= row["alignment_columns"]):
        raise ValueError("invalid alignment length")
    row["identity"] = row["matching_bases"] / row["alignment_columns"]
    return row


def qualifies(row: dict, context_name: str, array_start: int, array_end: int,
              flank: int, identity_min: float) -> bool:
    return (row["primary"] and row["target"] == context_name and
            row["target_start_0"] <= array_start - flank and
            row["target_end_0"] >= array_end + flank and
            row["identity"] >= identity_min)


def score(protocol: Path, selected: Path, context: Path, paf: Path, outdir: Path) -> dict:
    config = json.loads(protocol.read_text())
    chosen = json.loads(selected.read_text())
    context_hash = hashlib.sha256(context.read_bytes()).hexdigest()
    if context_hash != chosen["selected_context.fa_sha256"]:
        raise ValueError("selected context SHA-256 mismatch")
    item = chosen["selected"]
    if item is None:
        raise ValueError("no frozen selected interval")
    name = context.open().readline()[1:].split()[0]
    flank = config["candidate_generation_before_any_read_alignment"]["natural_flank_bp_each_side"]
    array_start = flank
    array_end = flank + item["array_end_0"] - item["array_start_0"]
    read_rule = config["read_alignment_only_after_selected_context_committed"]
    minflank = read_rule["spanner_min_natural_flank_bp_each_side"]
    identity_min = read_rule["spanner_min_identity_nmatch_over_alignment_columns"]
    rows = [parse_paf(line) for line in paf.open() if line.strip()]
    primary = [row for row in rows if row["primary"]]
    good = [row for row in rows if qualifies(row, name, array_start, array_end,
                                             minflank, identity_min)]
    best = {}
    for row in good:
        old = best.get(row["read_id"])
        if old is None or (row["identity"], row["mapq"], row["target_end_0"] - row["target_start_0"]) > (
            old["identity"], old["mapq"], old["target_end_0"] - old["target_start_0"]):
            best[row["read_id"]] = row
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "qualifying_spanners.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["read_id", "read_length", "query_start_0", "query_end_0", "strand", "target_start_0", "target_end_0", "matching_bases", "alignment_columns", "identity", "mapq"])
        for row in sorted(best.values(), key=lambda r: r["read_id"]):
            writer.writerow([row[k] for k in ("read_id", "read_length", "query_start_0", "query_end_0", "strand", "target_start_0", "target_end_0", "matching_bases", "alignment_columns")] + [f"{row['identity']:.9f}", row["mapq"]])
    result = {
        "selected_context_sha256": context_hash,
        "context_paf_sha256": hashlib.sha256(paf.read_bytes()).hexdigest(),
        "reported_alignment_count": len(rows),
        "reported_primary_alignment_count": len(primary),
        "reported_distinct_primary_reads": len({r["read_id"] for r in primary}),
        "qualifying_distinct_spanners": len(best),
        "required_distinct_spanners": read_rule["spanner_min_distinct_original_records"],
        "identity_min": identity_min,
        "natural_flank_bp_each_side_min": minflank,
        "pass_for_full_reference_step": len(best) >= read_rule["spanner_min_distinct_original_records"],
    }
    (outdir / "context_score.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(score(args.protocol, args.selected, args.context, args.paf, args.outdir)))


if __name__ == "__main__":
    main()
