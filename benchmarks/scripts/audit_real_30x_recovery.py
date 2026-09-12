"""Read-only integrity audit for interrupted real-read 30x material blocks.

This is deliberately a recovery *gate*, not a resume controller.  It does not
modify the material root or attempt to reconstruct a receipt.  In particular,
a directory with useful nonzero chunks but a zero-byte ``block_receipt.json``
is ``recoverable_invalid``: its data may help a later repair, but it is not a
valid benchmark input or result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


RECEIPT_NAMES = ("block_receipt.json", "sampling_receipt.json", "run_receipt.json",
                 "execution.json", "workflow_receipt.json")
SHA256_LENGTH = 64
MAX_RECEIPT_BYTES = 10 << 20


def _sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == SHA256_LENGTH and all(
        char in "0123456789abcdef" for char in value.lower()
    )


def _contains_sha256(value: Any) -> bool:
    if isinstance(value, dict):
        return any("sha256" in str(key).lower() and _sha256(item)
                   for key, item in value.items()) or any(_contains_sha256(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_sha256(item) for item in value)
    return False


def _positive_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0


def _find_first_positive(mapping: dict[str, Any], names: tuple[str, ...]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _stream_sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    """Hash a file without retaining a potentially multi-GB artifact in RAM."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_file_reasons(receipt: dict[str, Any], directory: Path) -> list[str]:
    """Check explicitly declared required artifacts without guessing a schema."""
    declared = receipt.get("required_files", receipt.get("required_artifacts", []))
    if declared is None:
        return []
    if not isinstance(declared, list):
        return ["declared_required_files_not_list"]
    reasons: list[str] = []
    for item in declared:
        if isinstance(item, str):
            relative, expected_sha = item, None
        elif isinstance(item, dict):
            relative = item.get("path", item.get("file"))
            expected_sha = item.get("sha256")
        else:
            reasons.append("declared_required_file_invalid_entry")
            continue
        if not isinstance(relative, str) or not relative:
            reasons.append("declared_required_file_missing_path")
            continue
        candidate = (directory / relative).resolve()
        if directory.resolve() not in candidate.parents:
            reasons.append(f"declared_required_file_escapes_block:{relative}")
        elif not candidate.is_file() or candidate.stat().st_size == 0:
            reasons.append(f"declared_required_file_missing_or_zero:{relative}")
        elif expected_sha is not None:
            observed = _stream_sha256(candidate)
            if not _sha256(expected_sha) or observed != expected_sha:
                reasons.append(f"declared_required_file_sha256_mismatch:{relative}")
    return reasons


def _receipt_gate(receipt: dict[str, Any], directory: Path) -> tuple[bool, list[str], dict[str, Any]]:
    """Validate only portable completion facts; unknown schemas stay invalid."""
    reasons: list[str] = []
    complete = receipt.get("complete") is True or receipt.get("status") in {"complete", "completed", "ok"}
    if not complete:
        reasons.append("receipt_does_not_declare_complete")
    reads = _find_first_positive(receipt, ("read_count", "record_count", "records", "selected_reads"))
    bases = _find_first_positive(receipt, ("total_bases", "base_count", "bases", "selected_bases"))
    if not _positive_number(reads):
        reasons.append("missing_positive_record_total")
    if not _positive_number(bases):
        reasons.append("missing_positive_base_total")
    if not _contains_sha256(receipt):
        reasons.append("missing_valid_sha256_field")
    expected_reads = _find_first_positive(receipt, ("expected_read_count", "planned_read_count"))
    expected_bases = _find_first_positive(receipt, ("expected_total_bases", "planned_total_bases", "target_bases"))
    if _positive_number(expected_reads) and _positive_number(reads) and expected_reads != reads:
        reasons.append("declared_record_total_mismatch")
    if _positive_number(expected_bases) and _positive_number(bases) and expected_bases != bases:
        reasons.append("declared_base_total_mismatch")
    reasons.extend(_required_file_reasons(receipt, directory))
    return not reasons, reasons, {"read_count": reads, "total_bases": bases,
                                   "expected_read_count": expected_reads,
                                   "expected_total_bases": expected_bases}


def _read_receipt(path: Path) -> tuple[dict[str, Any] | None, list[str]]:
    if not path.exists():
        return None, ["receipt_missing"]
    try:
        if path.stat().st_size == 0:
            return None, ["receipt_zero_bytes"]
        if path.stat().st_size > MAX_RECEIPT_BYTES:
            return None, ["receipt_exceeds_bounded_max_bytes"]
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return None, [f"receipt_unreadable_or_invalid_json:{type(exc).__name__}"]
    if not isinstance(data, dict):
        return None, ["receipt_not_json_object"]
    return data, []


def _payload_files(directory: Path, receipt: Path | None, *, max_depth: int,
                   max_entries: int) -> tuple[list[Path], str | None]:
    """Bounded local payload inspection; never walk an entire mounted results tree."""
    payload: list[Path] = []
    pending: list[tuple[Path, int]] = [(directory, 0)]
    visited = 0
    try:
        while pending:
            parent, depth = pending.pop()
            for path in parent.iterdir():
                visited += 1
                if visited > max_entries:
                    return payload, "payload_entry_limit_exceeded"
                if path.is_symlink():
                    return payload, "payload_symlink_not_followed"
                if path.is_file() and path != receipt and path.stat().st_size > 0:
                    payload.append(path)
                elif path.is_dir() and depth < max_depth:
                    pending.append((path, depth + 1))
    except OSError as exc:
        return payload, f"payload_scan_oserror:{type(exc).__name__}"
    return payload, None


def audit_block(directory: Path, root: Path, *, payload_max_depth: int = 2,
                payload_max_entries: int = 10_000) -> dict[str, Any]:
    receipts = [directory / name for name in RECEIPT_NAMES if (directory / name).exists()]
    receipt = next((path for path in receipts if path.name == "block_receipt.json"), receipts[0] if receipts else None)
    payload, payload_scan_issue = _payload_files(directory, receipt, max_depth=payload_max_depth,
                                                  max_entries=payload_max_entries)
    result: dict[str, Any] = {
        "relative_path": str(directory.relative_to(root)),
        "receipt": str(receipt.relative_to(root)) if receipt else None,
        "receipt_size_bytes": receipt.stat().st_size if receipt and receipt.exists() else None,
        "nonzero_payload_file_count": len(payload),
        "nonzero_payload_bytes": sum(path.stat().st_size for path in payload),
        "status": "missing",
        "reasons": [],
        "payload_scan_status": payload_scan_issue or "complete",
    }
    if payload_scan_issue and receipt is None:
        result["status"] = "blocked"
        result["reasons"] = [payload_scan_issue, "receipt_missing_payload_state_not_fully_observed"]
        return result
    if receipt is None:
        result["status"] = "recoverable_invalid" if payload else "missing"
        result["reasons"] = ["receipt_missing"] + (["nonzero_payload_requires_regenerated_receipt"] if payload else [])
        return result
    parsed, parse_reasons = _read_receipt(receipt)
    if parsed is None:
        result["status"] = "blocked" if payload_scan_issue else "recoverable_invalid" if payload else "corrupt"
        result["reasons"] = parse_reasons + ([payload_scan_issue] if payload_scan_issue else []) + (["nonzero_payload_requires_regenerated_receipt"] if payload else [])
        return result
    valid, reasons, totals = _receipt_gate(parsed, directory)
    result["declared_totals"] = totals
    if valid:
        result["status"] = "valid"
        result["receipt_sha256"] = _stream_sha256(receipt)
        return result
    result["status"] = "blocked" if payload_scan_issue else "recoverable_invalid" if payload else "corrupt"
    result["reasons"] = reasons + ([payload_scan_issue] if payload_scan_issue else []) + (["nonzero_payload_requires_regenerated_receipt"] if payload else [])
    return result


def audit_root(root: Path, blocks: list[Path] | tuple[Path, ...] = (), manifest: Path | None = None,
               *, payload_max_depth: int = 2, payload_max_entries: int = 10_000) -> dict[str, Any]:
    """Return serializable material/block states without modifying ``root``."""
    result: dict[str, Any] = {"schema_version": 1, "root": str(root), "root_status": "ok", "blocks": []}
    if not root.exists():
        result["root_status"] = "blocked"
        result["root_reason"] = "input_root_missing_or_not_mounted"
        return result
    if not root.is_dir():
        result["root_status"] = "blocked"
        result["root_reason"] = "input_root_is_not_directory"
        return result
    selected: list[Path] = []
    for raw_path in blocks:
        directory = raw_path.resolve() if raw_path.is_absolute() else (root / raw_path).resolve()
        if root not in directory.parents and directory != root:
            result["blocks"].append({"relative_path": str(raw_path), "status": "blocked",
                                     "reasons": ["explicit_block_escapes_root"]})
        elif not directory.is_dir():
            result["blocks"].append({"relative_path": str(raw_path), "status": "missing",
                                     "reasons": ["explicit_block_directory_missing"]})
        else:
            selected.append(directory)
    result["blocks"].extend(audit_block(directory, root, payload_max_depth=payload_max_depth,
                                         payload_max_entries=payload_max_entries) for directory in selected)
    if not blocks:
        result["root_status"] = "missing"
        result["root_reason"] = "no_explicit_block_paths_provided"
    if manifest is not None:
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            result["material_inventory"] = [
                {"id": item.get("id"), "target_30x_bp": item.get("target_30x_bp"), "state": item.get("state")}
                for item in data.get("materials", []) if isinstance(item, dict)
            ]
        except (OSError, json.JSONDecodeError) as exc:
            result["manifest_status"] = f"blocked:{type(exc).__name__}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("tmp/real_30x_recovery_input"),
                        help="Read-only root containing material blocks; default is an internal placeholder, never T7.")
    parser.add_argument("--manifest", type=Path, default=Path("docs/evidence/real_30x_recovery/recovery_manifest_20260913.json"))
    parser.add_argument("--block", type=Path, action="append", default=[],
                        help="Explicit block directory, relative to --root; repeat as needed. No root-wide discovery occurs.")
    parser.add_argument("--payload-max-depth", type=int, default=2,
                        help="Maximum depth below each explicit block for payload evidence (default: 2).")
    parser.add_argument("--payload-max-entries", type=int, default=10_000,
                        help="Maximum metadata entries inspected below each explicit block (default: 10000).")
    parser.add_argument("--out", type=Path, default=Path("tmp/real_30x_recovery_audit.json"),
                        help="Output JSON; must not be inside --root.")
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.out.resolve()
    if out == root or root == out.parent or root in out.parents:
        raise ValueError("--out must be outside the read-only --root")
    out.parent.mkdir(parents=True, exist_ok=True)
    if args.payload_max_depth < 0 or args.payload_max_entries < 1:
        raise ValueError("payload bounds must be positive (depth may be zero)")
    result = audit_root(root, args.block, args.manifest.resolve() if args.manifest.exists() else None,
                        payload_max_depth=args.payload_max_depth, payload_max_entries=args.payload_max_entries)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"root_status": result["root_status"], "block_count": len(result["blocks"]), "out": str(out)}))


if __name__ == "__main__":
    main()
