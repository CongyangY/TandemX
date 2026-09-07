from __future__ import annotations

import hashlib
from pathlib import Path

import yaml


def test_cascade_heldout_config_freezes_unseen_seeds_and_development_evidence() -> None:
    root = Path(__file__).resolve().parents[2]
    config = yaml.safe_load((root / "benchmarks/configs/cascade_native_screen_heldout_v1.yaml").read_text())
    seeds = config["seeds"]
    assert set(seeds["heldout"]).isdisjoint(seeds["development"] + seeds["validation"])
    assert seeds["heldout"] == [3101, 3102, 3103]
    assert config["discovery_method"] == "cascade"
    assert set(config["tools"]) == {"tandemx", "trf", "tidehunter"}
    provenance = config["selection_provenance"]
    for name, expected in (
        ("archive_manifest.json", provenance["development_archive_manifest_sha256"]),
        ("aggregate.json", provenance["development_aggregate_sha256"]),
    ):
        observed = hashlib.sha256(
            (root / "paper/evidence/cascade_native_screen_development" / name).read_bytes()
        ).hexdigest()
        assert observed == expected
    assert config["acceptance_gates"]["warning"].startswith("promotion_requires_all_gates")
