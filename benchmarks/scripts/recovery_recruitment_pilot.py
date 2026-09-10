"""Non-inferential first-1000-read cost gate for an enrolled recovery target bank."""
from __future__ import annotations

import argparse
import json
from itertools import islice
from pathlib import Path

from tandemx.io.sequences import read_sequence_records
from tandemx.recovery.poc import run_mapping, sha256, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    args = parser.parse_args()
    root = args.prepared
    state = json.loads((root / "prepared.json").read_text())
    target = root / "recruitment_targets.fa"
    if sha256(target) != state["targets_sha256"]:
        raise ValueError("Frozen target bank changed")
    out = root / "recruitment_cost_pilot"
    out.mkdir(exist_ok=False)
    reads = out / "prefix1000.fa"
    config = state["config"]
    count = bases = 0
    with reads.open("w") as handle:
        for record in islice(read_sequence_records(Path(config["inputs"]["reads"])), 1000):
            count += 1
            bases += len(record.sequence)
            handle.write(f">{record.id}\n{record.sequence}\n")
    write_json(out / "input.json", dict(reads=count, bases=bases, sha256=sha256(reads),
               selection="first_1000_records_noninferential_resource_pilot", target_sha256=sha256(target)))
    run_mapping([config["minimap2"], "-x", "map-hifi", "-c", "-N", "50", "-p", "0.5", "--secondary=yes",
                 "-K", "50M", "-t", str(config.get("threads", 4)), str(target), str(reads)], out / "pilot.paf", timeout=120)


if __name__ == "__main__":
    main()
