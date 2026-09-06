import json
from pathlib import Path
import random

import pytest

from benchmarks.challenge.run import source_manifest
from benchmarks.scripts.replay_clustering_isolated import replay
from tandemx.discover.mvp import CandidateRepeat, write_candidate_reads
from tandemx.discover.rust_backend import rust_backend_available


@pytest.mark.skipif(not rust_backend_available(), reason='Native kernel required')
def test_isolated_replay_runs_frozen_children_and_rejects_changed_baseline(tmp_path):
    root = Path(__file__).resolve().parents[2]
    previous = tmp_path/'previous'
    previous.mkdir()
    environment = source_manifest(root, previous/'source_snapshot')
    (previous/'environment.json').write_text(json.dumps(environment))
    folder = previous/'tandemx/discover'
    folder.mkdir(parents=True)
    rng = random.Random(6101)
    monomer = ''.join(rng.choices('ACGT', k=61))
    variants = [monomer, monomer[7:]+monomer[:7], ''.join(rng.choices('ACGT', k=61))]
    candidates = [CandidateRepeat(f'r{i}', f'TXC{i:06d}', seq, 0, 610, '+', 61, 610, 10, 1.0,
                                 False, 'medium', '') for i, seq in enumerate(variants)]
    write_candidate_reads(folder/'candidate_reads.tsv', candidates)
    (folder/'candidate_monomers.fa').write_text(''.join(f'>candidate_id={c.candidate_id}\n{c.sequence}\n' for c in candidates))
    result = replay(previous, previous, tmp_path/'result', 60)
    assert result['complete'] and result['exact_output_parity']
    assert all(r['candidate_count'] == 3 and r['family_count'] == 2 and r['peak_rss_mib'] > 0 for r in result['measurements'])
    ablation = replay(previous, None, tmp_path/'ablation', 60)
    assert ablation['complete'] and ablation['exact_output_parity']
    assert [r['index_backend'] for r in ablation['measurements']] == ['python', 'native']
    assert ablation['measurements'][0]['clustering_source_sha256'] == ablation['measurements'][1]['clustering_source_sha256']
    interface = replay(previous, None, tmp_path/'interface', 60, interface_ablation=True)
    assert interface['complete'] and interface['exact_output_parity']
    assert [r['index_backend'] for r in interface['measurements']] == ['word_bridge', 'sequence_native']
    assert interface['measurements'][0]['output_sha256'] == interface['measurements'][1]['output_sha256']
    (previous/'source_snapshot/tandemx/discover/distance.py').write_text('changed')
    with pytest.raises(ValueError, match='snapshot changed'):
        replay(previous, previous, tmp_path/'invalid', 60)
