import csv
import json
import random
import subprocess
import sys

import pytest

from tandemx.discover.rust_backend import rust_backend_available
from tandemx.io.validators import validate_project
from tandemx.utils.threads import discover_thread_limit


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
    for replicate, threads in enumerate([1, min(2, discover_thread_limit())]):
        out = tmp_path / str(replicate)
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
        members = list(csv.DictReader((folder / "monomer_membership.tsv").open(), delimiter="\t"))
        assert len(members) == 12
        assert all(row["status"] == "assigned" and float(row["similarity_lower_bound"]) >= 0.95 for row in members)
        assert (folder / "candidate_monomers.fa").read_text().count(">candidate_id=") == 12
        assert validate_project(folder)
        assert 'discovery_method: "elastic"' in (folder / "run_config.yaml").read_text()
        artifacts.append([(folder / name).read_bytes() for name in
                          ["candidate_reads.tsv", "candidate_monomers.fa", "monomer_membership.tsv",
                           "monomers.fa", "families.tsv", "family_similarity.tsv"]])
    assert artifacts[0] == artifacts[1]


def test_related_family_audit_pipeline_preserves_catalog_and_exposes_omissions(tmp_path):
    rng = random.Random(6205)
    a = ''.join(rng.choices('ACGT', k=61))
    b = ''.join(rng.choices('ACGT', k=93))
    reads = tmp_path/'reads.fa'
    reads.write_text(''.join(f'>r{i}\n{a*8}'+ 'N'*100+f'{b*5}\n' for i in range(3)))
    outputs = []
    for mode in ('full', 'related'):
        out = tmp_path/mode
        result = subprocess.run([sys.executable, '-m', 'tandemx.cli', 'run', '--reads', str(reads),
            '--outdir', str(out), '--steps', 'discover,validate', '--min-period', '30', '--max-period', '120',
            '--discovery-method', 'elastic', '--kmer-backend', 'rust', '--threads', '1', '--family-audit', mode],
            text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
        folder = out/'discover'
        outputs.append(folder)
        summary = json.loads((folder/'family_audit_summary.json').read_text())
        assert summary['mode'] == mode and summary['complete']
        assert f'family_audit: "{mode}"' in (folder/'run_config.yaml').read_text()
        assert 'family_audit_summary.json' in json.loads((folder/'discovery_summary.json').read_text())['output_sha256']
    for filename in ['candidate_reads.tsv','candidate_monomers.fa','monomers.fa','families.tsv','monomer_membership.tsv']:
        assert (outputs[0]/filename).read_bytes() == (outputs[1]/filename).read_bytes()
    full = list(csv.DictReader((outputs[0]/'family_similarity.tsv').open(), delimiter='\t'))
    related = list(csv.DictReader((outputs[1]/'family_similarity.tsv').open(), delimiter='\t'))
    assert related == [r for r in full if r['relationship'] != 'distinct']
