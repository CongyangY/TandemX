#!/usr/bin/env python3
"""Independently verify and inventory the enrolled ZH13-T2T assembly."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


def inspect(path: Path) -> tuple[dict[str, Any], list[tuple[str, int]]]:
    md5 = hashlib.md5()  # noqa: S324 - publisher checksum verification
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            size += len(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    names: set[str] = set()
    lengths: list[tuple[str, int]] = []
    name: str | None = None
    length = acgt = other = 0
    total_acgt = total_other = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if length == 0:
                        raise ValueError(f"empty FASTA record before line {line_number}")
                    lengths.append((name, length))
                    total_acgt += acgt
                    total_other += other
                name = line[1:].split()[0]
                if not name or name in names:
                    raise ValueError(f"empty or duplicate FASTA ID at line {line_number}")
                names.add(name)
                length = acgt = other = 0
            else:
                if name is None:
                    raise ValueError(f"sequence before FASTA header at line {line_number}")
                upper = line.upper()
                line_acgt = sum(upper.count(base) for base in "ACGT")
                length += len(line)
                acgt += line_acgt
                other += len(line) - line_acgt
    if name is None or length == 0:
        raise ValueError("empty or incomplete FASTA")
    lengths.append((name, length))
    total_acgt += acgt
    total_other += other
    total_bases = total_acgt + total_other
    summary = {
        "bytes": size,
        "md5": md5.hexdigest(),
        "sha256": sha256.hexdigest(),
        "record_count": len(lengths),
        "base_count": total_bases,
        "acgt_bases": total_acgt,
        "other_bases": total_other,
    }
    return summary, lengths


def verify(config_path: Path, result_dir: Path) -> dict[str, Any]:
    output = result_dir / "independent_verification.json"
    lengths_path = result_dir / "sequence_lengths.tsv"
    if output.exists() or lengths_path.exists():
        raise FileExistsError("refusing to overwrite ZH13 independent verification")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt = json.loads((result_dir / "run_receipt.json").read_text(encoding="utf-8"))
    file_record = receipt.get("compressed_file", {})
    payload = result_dir / file_record.get("file", "")
    summary, lengths = inspect(payload)
    if (
        receipt.get("complete") is not True
        or receipt.get("fate") != "source_enrollment_passed"
        or summary["bytes"] != config["assembly"]["expected_bytes"]
        or summary["md5"] != config["assembly"]["expected_md5"]
        or summary != {
            "bytes": file_record["bytes"],
            "md5": file_record["md5"],
            "sha256": file_record["sha256"],
            "record_count": receipt["fasta_qc"]["record_count"],
            "base_count": receipt["fasta_qc"]["base_count"],
            "acgt_bases": receipt["fasta_qc"]["acgt_bases"],
            "other_bases": receipt["fasta_qc"]["other_bases"],
        }
    ):
        raise ValueError("ZH13 assembly differs between independent and enrollment checks")
    with lengths_path.open("x", encoding="utf-8") as handle:
        handle.write("sequence_id\tlength_bp\n")
        for name, length in lengths:
            handle.write(f"{name}\t{length}\n")
    result = {
        "schema_version": 1,
        "complete": True,
        "independent_of_enrollment_implementation": True,
        "assembly": summary,
        "sequence_lengths_file": {
            "file": lengths_path.name,
            "bytes": lengths_path.stat().st_size,
            "sha256": hashlib.sha256(lengths_path.read_bytes()).hexdigest(),
        },
        "boundary": (
            "This verifies source bytes, publisher MD5, gzip/FASTA readability and "
            "sequence lengths only. It is not unitFinder output or accuracy evidence."
        ),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--result-dir", required=True, type=Path)
    args = parser.parse_args()
    verify(args.config, args.result_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
