import json
import random
import subprocess
import sys

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.replay_discovery import replay
from tandemx.discover.rust_backend import rust_backend_available


@pytest.mark.skipif(not rust_backend_available(), reason='Requires native discovery')
def test_live_discovery_profile_replay_and_changed_input_rejection(tmp_path):
    previous = tmp_path/'previous'
    previous.mkdir()
    rng = random.Random(71022)
    unit = ''.join(rng.choices('ACGT', k=61))
    fasta = previous/'reads.fa'
    fasta.write_text(''.join(f'>r{i}\n{unit*12}\n' for i in range(3)))
    output = previous/'tandemx/discover'
    command = [sys.executable, '-m', 'tandemx.cli', 'discover', '--reads', str(fasta), '--outdir', str(output),
               '--min-period', '30', '--max-period', '1000', '--min-repeat-span', '100', '--min-support-reads', '1',
               '--min-read-length', '1', '--kmer-backend', 'rust', '--threads', '1', '--no-progress',
               '--discovery-method', 'elastic', '--clustering-method', 'sequence', '--family-audit', 'related']
    subprocess.run(command, check=True, capture_output=True)
    (previous/'environment.json').write_text(json.dumps(dict(input=dict(fasta_sha256=digest_file(fasta)))))
    (previous/'tandemx/execution.json').write_text(json.dumps(dict(command=command, exit_code=0, timed_out=False)))
    result = replay(previous, tmp_path/'replay', 60, profile=True)
    assert result['complete'] and len(result['products']) == 7
    assert all(r['byte_identical'] for r in result['products'].values())
    assert (tmp_path/'replay/discovery.prof').stat().st_size > 0
    fasta.write_text('>changed\nACGT\n')
    with pytest.raises(ValueError, match='no longer matches'):
        replay(previous, tmp_path/'invalid')
    assert not (tmp_path/'invalid').exists()
