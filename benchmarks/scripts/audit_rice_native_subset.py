"""Read back the frozen bounded rice FASTQ and record structural integrity."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


class HashingFile:
    def __init__(self, path: Path) -> None:
        self.file = path.open("rb")
        self.digest = hashlib.sha256()
        self.bytes_read = 0

    def read(self, n: int = -1) -> bytes:
        data = self.file.read(n)
        self.digest.update(data)
        self.bytes_read += len(data)
        return data

    def close(self) -> None:
        self.file.close()


def audit(path: Path, expected_sha256: str) -> dict:
    source = HashingFile(path)
    names: set[bytes] = set()
    count = total_bases = 0
    try:
        with gzip.GzipFile(fileobj=source, mode="rb") as stream:
            while header := stream.readline():
                sequence = stream.readline().rstrip(b"\r\n")
                plus = stream.readline()
                quality = stream.readline().rstrip(b"\r\n")
                if not header.startswith(b"@") or not plus.startswith(b"+"):
                    raise ValueError(f"invalid FASTQ framing at record {count + 1}")
                if len(sequence) != len(quality) or not sequence:
                    raise ValueError(f"invalid FASTQ sequence/quality at record {count + 1}")
                name = header[1:].split()[0]
                if name in names:
                    raise ValueError(f"duplicate read ID at record {count + 1}")
                names.add(name)
                count += 1
                total_bases += len(sequence)
        # GzipFile may stop at the final member before requesting an EOF read.
        for _ in iter(lambda: source.read(1024 * 1024), b""):
            pass
        digest = source.digest.hexdigest()
        if digest != expected_sha256:
            raise ValueError(f"SHA-256 mismatch: {digest}")
        return {"path": str(path), "size_bytes": path.stat().st_size,
                "sha256": digest, "gzip_trailer_valid": True,
                "fastq_structure_valid": True, "read_count": count,
                "total_bases": total_bases, "unique_read_ids": len(names)}
    finally:
        source.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.protocol.read_text())
    receipt = audit(Path(config["original_read_subset_path"]), config["subset_sha256_expected"])
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
