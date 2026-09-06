import json

import pytest

from benchmarks.scripts.generate_factorial_scale import run as generate
from benchmarks.abundance.run_stream_quantify import run
from benchmarks.scripts.evaluate_multik_factorial import worker as multik_worker
from benchmarks.challenge.schema import read_table


def test_factorial_quantification_runs_real_cli_with_truth_only_in_evaluator(tmp_path):
    histogram=tmp_path/'lengths.tsv';histogram.write_text('length_bp\tread_count\n200\t1\n')
    config=dict(seeds={'development':[6361],'heldout':[7361]},genome_bp=3000,background_gc=.45,
        periods=[31],copies=[3,7],gc_fractions=[.5],unit_substitution_rates=[0,.05],coverages=[3],
        read_error_models=[dict(label='clean',substitution_rate=0,insertion_rate=0,deletion_rate=0)])
    config_path=tmp_path/'config.json';config_path.write_text(json.dumps(config))
    dataset=tmp_path/'dataset'
    generate(config_path,histogram,dataset,6361,100000)
    out=tmp_path/'evaluated';run([dataset],out,30)
    receipt=json.loads((out/'validation.json').read_text())
    assert receipt['complete'] and receipt['executions']==receipt['successful']==1
    assert receipt['copy_number_family_conditions']==4
    metrics=read_table(out/'copy_number_metrics.tsv')
    assert {r['truth_copies'] for r in metrics}=={'3','7'}
    assert {r['unit_substitution_rate'] for r in metrics}=={'0.0','0.05'}
    command=json.loads((out/'runs/s6361/condition_001/execution.json').read_text())['command']
    assert '--genome-size' in command and '--catalog' in command
    assert not any('truth_copy' in arg or 'sampling.tsv' in arg or 'substitution' in arg for arg in command)
    replay=tmp_path/'multik';replay.mkdir()
    multik_worker(out,replay)
    replay_receipt=json.loads((replay/'validation.json').read_text())
    assert replay_receipt['complete'] and replay_receipt['family_conditions']==4
    assert replay_receipt['paired_method_rows']==12 and not replay_receipt['heldout_used']
    assert len(read_table(replay/'per_k.tsv'))==16
    assert {r['method'] for r in read_table(replay/'metrics.tsv')}=={
        'native_median_k21','mean_k21_exposure','multik_loglinear'}
    with pytest.raises(ValueError,match='unique'):
        run([dataset,dataset],tmp_path/'duplicate')
    (dataset/'genome/catalogue.fa').write_text('changed catalogue')
    with pytest.raises(ValueError,match='hash mismatch'):
        run([dataset],tmp_path/'changed')


def test_stream_quantifier_refuses_incomplete_or_reserved_generation(tmp_path):
    source=tmp_path/'dataset';source.mkdir()
    receipt=source/'generation_receipt.json'
    for record in [dict(seed=1,complete=False,split='development',heldout_used=False,conditions_completed=[{}]),
                   dict(seed=1,complete=True,split='heldout',heldout_used=True,conditions_completed=[{}]),
                   dict(seed=1,complete=True,split='development',heldout_used=False,conditions_completed=[])]:
        receipt.write_text(json.dumps(record))
        with pytest.raises(ValueError,match='development'):
            run([source],tmp_path/'invalid')
        assert not (tmp_path/'invalid').exists()
    with pytest.raises(ValueError,match='positive timeout'):
        run([source],tmp_path/'invalid',float('nan'))
