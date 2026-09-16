"""Opt-in exact scan snapshots for direct discovery."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, fields
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tandemx.discover.mvp import CandidateRepeat, DiscoverConfig


def _encoded(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(config: DiscoverConfig, reads: Path) -> dict[str, object]:
    # Include operational parameters too: a changed command must never silently reuse state.
    parameters = {}
    for field in fields(config):
        if field.name == "reads":
            continue
        value = getattr(config, field.name)
        parameters[field.name] = str(value) if isinstance(value, Path) else value
    return {
        "parameters": parameters,
        "reads_path": str(reads.resolve()),
        "reads_size": reads.stat().st_size,
        "reads_sha256": _sha_file(reads),
    }


def update_prefix(digest: object, read_id: str, sequence: str) -> None:
    # Length framing removes ambiguities between adjacent IDs and sequences.
    for value in (read_id, sequence):
        data = value.encode("utf-8")
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)


def save(path: Path, body: dict[str, object]) -> None:
    envelope = {"body": body, "body_sha256": hashlib.sha256(_encoded(body)).hexdigest()}
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(_encoded(envelope) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load(path: Path, expected_identity: dict[str, object], candidate_path: Path) -> tuple[dict[str, object], list[CandidateRepeat]]:
    from tandemx.discover.mvp import CandidateRepeat, candidate_reads_header, format_candidate_read

    try:
        envelope = json.loads(path.read_text(encoding="utf-8"))
        body = envelope["body"]
        if envelope["body_sha256"] != hashlib.sha256(_encoded(body)).hexdigest():
            raise ValueError("checksum mismatch")
        if body["version"] != 1 or body["identity"] != expected_identity:
            raise ValueError("input or configuration changed")
        for name in ("processed_reads", "processed_bases", "skipped_short_reads", "skipped_short_kmer", "skipped_low_complexity", "seed_overflow_count"):
            if type(body[name]) is not int or body[name] < 0:
                raise ValueError(f"invalid {name}")
        if body["processed_reads"] == 0 or not isinstance(body["prefix_sha256"], str) or len(body["prefix_sha256"]) != 64:
            raise ValueError("invalid read prefix")
        raw_candidates = body["candidates"]
        if not isinstance(raw_candidates, list):
            raise ValueError("invalid candidates")
        names = {field.name for field in fields(CandidateRepeat)}
        if any(not isinstance(item, dict) or set(item) != names for item in raw_candidates):
            raise ValueError("invalid candidate fields")
        candidates = [CandidateRepeat(**item) for item in raw_candidates]
        if any(candidate.candidate_id != f"TXC{index:06d}" for index, candidate in enumerate(candidates, 1)):
            raise ValueError("invalid candidate numbering")
        table_digest = hashlib.sha256()
        with candidate_path.open("rb") as handle:
            for expected in chain((candidate_reads_header() + "\n",), (format_candidate_read(item) + "\n" for item in candidates)):
                actual = handle.readline()
                if actual != expected.encode("utf-8"):
                    raise ValueError("candidate table changed or truncated")
                table_digest.update(actual)
        if body["candidate_prefix_sha256"] != table_digest.hexdigest():
            raise ValueError("candidate table checksum mismatch")
        return body, candidates
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid discover scan checkpoint {path}: {exc}") from exc


def snapshot_body(identity_value: dict[str, object], prefix_sha256: str, candidate_path: Path,
                  candidates: list[CandidateRepeat], **counts: int) -> dict[str, object]:
    return {
        "version": 1,
        "identity": identity_value,
        "prefix_sha256": prefix_sha256,
        "candidate_prefix_sha256": _sha_file(candidate_path),
        "candidates": [asdict(candidate) for candidate in candidates],
        **counts,
    }
