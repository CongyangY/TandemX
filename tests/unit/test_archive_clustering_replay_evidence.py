import hashlib
import json

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_clustering_replay_evidence import LABELS, SOURCE_FILES, archive


def build_replay(tmp_path):
    source = tmp_path / "replay"
    snapshot = source / "snapshot"
    snapshot.mkdir(parents=True)
    hashes = {}
    helpers = {}
    for name in SOURCE_FILES:
        path = snapshot / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"source={name}\n")
        if name.startswith("benchmarks/"):
            helpers[name] = digest_file(path)
        else:
            hashes[name] = digest_file(path)
    source_digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    environment = {
        "comparison": "native_sequence_interface_ablation",
        "source_snapshot": str(snapshot),
        "file_hashes": hashes,
        "helper_hashes": helpers,
        "source_digest": source_digest,
    }
    (source / "environment.json").write_text(json.dumps(environment))
    payload = b'{"families":[],"membership":[]}'
    payload_hash = hashlib.sha256(payload).hexdigest()
    measurements = []
    for label in LABELS:
        (source / f"{label}.json").write_bytes(payload)
        receipt = {"output_sha256": payload_hash, "index_backend": label}
        (source / f"{label}.receipt.json").write_text(json.dumps(receipt))
        (source / f"{label}.stdout.log").write_text("")
        (source / f"{label}.stderr.log").write_text("")
        measurements.append({"label": label, "index_backend": label, "exit_code": 0,
                             "timed_out": False, "output_sha256": payload_hash})
    (source / "validation.json").write_text(json.dumps(
        {"complete": True, "exact_output_parity": True, "measurements": measurements}
    ))
    return source


def test_archive_validates_payloads_sources_and_excludes_large_json(tmp_path):
    source = build_replay(tmp_path)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert len(manifest) == 2 + 3 * len(LABELS) + len(SOURCE_FILES)
    assert not (outdir / "word_bridge.json").exists()
    assert (outdir / "source_snapshot/rust-core/src/representative_index.rs").is_file()
    for row in manifest:
        assert digest_file(outdir / row["file"]) == row["sha256"]
    with pytest.raises(ValueError, match="already exists"):
        archive(source, outdir)


@pytest.mark.parametrize("fault", ["incomplete", "payload", "source"])
def test_archive_rejects_incomplete_or_changed_evidence(tmp_path, fault):
    source = build_replay(tmp_path)
    if fault == "incomplete":
        validation = json.loads((source / "validation.json").read_text())
        validation["complete"] = False
        (source / "validation.json").write_text(json.dumps(validation))
    elif fault == "payload":
        (source / "sequence_native.json").write_text("changed")
    else:
        (source / "snapshot/tandemx/discover/clustering.py").write_text("changed")
    with pytest.raises(ValueError):
        archive(source, tmp_path / "archive")
