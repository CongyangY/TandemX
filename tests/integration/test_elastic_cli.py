import csv
import json
import random
import subprocess
import sys

import pytest

from tandemx.discover.rust_backend import rust_backend_available
from tandemx.io.validators import validate_project


@pytest.mark.parametrize("backend", ["python", "rust"])
def test_multi_array_cli_outputs_are_deterministic_and_valid(tmp_path, backend):
    if backend == "rust" and not rust_backend_available():
        pytest.skip("compiled extension unavailable")
    rng = random.Random(45)
    a = "".join(rng.choices("ACGT", k=61))
    b = "".join(rng.choices("ACGT", k=93))
    read = a * 8 + "N" * 120 + b * 5
    reads = tmp_path / "reads.fa"
    reads.write_text("".join(f">r{i}\n{read}\n" for i in range(6)))
    artifacts = []
    for threads in [1, 2]:
        out = tmp_path / str(threads)
        result = subprocess.run([sys.executable, "-m", "tandemx.cli", "run", "--reads", str(reads),
            "--outdir", str(out), "--steps", "discover,validate", "--min-period", "30", "--max-period", "120",
            "--discovery-method", "elastic", "--kmer-backend", backend, "--threads", str(threads)],
            text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
        folder = out / "discover"
        rows = list(csv.DictReader((folder / "candidate_reads.tsv").open(), delimiter="\t"))
        assert len(rows) == 12
        assert len({r["candidate_id"] for r in rows}) == 12
        assert {int(r["period_bp"]) for r in rows} == {61, 93}
        summary = json.loads((folder / "discovery_summary.json").read_text())
        assert summary["processed_reads"] == 6
        assert summary["candidate_count"] == 12
        assert summary["family_count"] == 2
        families = list(csv.DictReader((folder / "families.tsv").open(), delimiter="\t"))
        assert all("uncalibrated_confidence" in row["warning"] for row in families)
        assert validate_project(folder)
        assert 'discovery_method: "elastic"' in (folder / "run_config.yaml").read_text()
        artifacts.append([(folder / name).read_bytes() for name in
                          ["candidate_reads.tsv", "monomers.fa", "families.tsv", "family_similarity.tsv"]])
    assert artifacts[0] == artifacts[1]
