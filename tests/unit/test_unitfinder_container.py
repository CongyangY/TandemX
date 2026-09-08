from pathlib import Path


def test_unitfinder_container_pins_source_and_preserves_upstream_code() -> None:
    root = Path(__file__).resolve().parents[2]
    dockerfile = (root / "benchmarks/containers/unitfinder/Dockerfile").read_text()
    readme = (root / "benchmarks/containers/unitfinder/README.md").read_text()
    assert "UNITFINDER_COMMIT=e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28" in dockerfile
    assert "trf=4.09.1" in dockerfile
    assert "mummer4=4.0.1" in dockerfile
    assert "clustalo=1.2.4" in dockerfile
    assert "git -C /opt/unitFinder status --short" in dockerfile
    assert "(cd /opt/unitFinder" in dockerfile
    assert "COPY" not in dockerfile
    assert "sed -i" not in dockerfile
    assert "-t tandemx/unitfinder:e80bff38-v2 benchmarks/containers/unitfinder" in readme
    assert "original `e80bff38` image-tag attempt is retained" in readme
