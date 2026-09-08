#!/usr/bin/env python3
"""Run the frozen Ey15 old-to-new minimap2 audit with resource sampling."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.profile_stage_resources import (
    SAMPLE_FIELDS,
    STAGE_FIELDS,
    run_stage,
    write_table,
)


def run(config_path: Path, outdir: Path, interval: float) -> dict[str, object]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("status") != "frozen_before_alignment_output_inspection":
        raise ValueError("alignment audit is not frozen before output inspection")
    tool = config["tool"]
    binary = Path(tool["binary"])
    inputs = config["inputs"]
    old = Path(inputs["old_assembly"]["path"])
    new = Path(inputs["new_assembly"]["path"])
    for label, path in (("old_assembly", old), ("new_assembly", new)):
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = digest_file(path)
        if observed != inputs[label]["sha256"]:
            raise ValueError(f"{label} SHA-256 differs from preregistration")
    version = subprocess.run(
        [str(binary), "--version"], check=True, text=True, capture_output=True
    ).stdout.strip()
    source_commit = subprocess.run(
        ["git", "-C", str(binary.parent), "rev-parse", "HEAD"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    if version != tool["version"] or source_commit != tool["source_commit"]:
        raise ValueError("minimap2 version or source commit differs from preregistration")
    command = [str(binary), *config["arguments"], str(old), str(new)]
    outdir.mkdir(parents=True)
    (outdir / "environment.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "config_path": str(config_path.resolve()),
                "config_sha256": digest_file(config_path),
                "command": command,
                "tool_version": version,
                "tool_source_commit": source_commit,
                "input_sha256": {
                    "old_assembly": digest_file(old),
                    "new_assembly": digest_file(new),
                },
                "resource_scope": "aggregate_live_process_tree_sampled_on_host",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    resource_dir = outdir / "resources"
    resource_dir.mkdir()
    stage = {
        "name": "old_new_minimap2",
        "command": command,
        "cwd": None,
        "scratch_dir": outdir,
    }
    result, samples = run_stage(stage, resource_dir, interval)
    raw_stdout = resource_dir / "old_new_minimap2.stdout.log"
    paf = outdir / ("old_new.paf" if result["exit_code"] == 0 else "old_new.partial.paf")
    raw_stdout.replace(paf)
    result["stdout_log"] = paf.relative_to(outdir).as_posix()
    write_table(resource_dir / "samples.tsv", SAMPLE_FIELDS, samples)
    write_table(resource_dir / "stages.tsv", STAGE_FIELDS, [result])
    receipt = {
        "schema_version": 1,
        "complete": result["exit_code"] == 0,
        "exit_code": result["exit_code"],
        "paf": paf.name,
        "paf_bytes": paf.stat().st_size,
        "paf_sha256": digest_file(paf),
        "resource": result,
        "boundary": config["boundary"],
    }
    (outdir / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    receipt = run(args.config, args.outdir, args.interval)
    return 0 if receipt["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

