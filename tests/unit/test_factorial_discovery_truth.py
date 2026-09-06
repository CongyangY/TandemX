import random

import pytest

from benchmarks.abundance.discovery_truth import DiscoveryTruth, fasta_lengths, score_predictions
from benchmarks.challenge.schema import ArrayRecord
from benchmarks.challenge.sequence_metrics import cyclic_edit_similarity, cyclic_reaches_threshold, score_threshold_recovery


def test_threshold_decision_matches_independent_exact_edit_endpoint():
    rng=random.Random(87324)
    pairs=[('N'*20,'N'*20),('ACGT'*10,'ACGT'*8),('GATTACA','TGTAATC')]
    for _ in range(35):
        sequence=''.join(rng.choices('ACGTN',k=rng.randrange(8,30)))
        rotated=sequence[3:]+sequence[:3]
        mutated=rotated[:5]+'A'+rotated[6:]
        pairs.extend([(sequence,mutated),(sequence,rotated+'C'),(sequence,'A'*len(sequence))])
    for truth,prediction in pairs:
        exact=cyclic_edit_similarity(truth,prediction)
        for threshold in (.5,.8,.9,1):
            assert cyclic_reaches_threshold(truth,prediction,threshold)==(exact+1e-12>=threshold)
    metrics,rows=score_threshold_recovery(['GATTACA','TGTAATC'],{'F':'GATTACA','G':'GATTACA'},1)
    assert metrics['distinct_consensus_count']==1 and metrics['recovered_family_count']==1
    assert sum(r['recovered'] for r in rows)==1
    with pytest.raises(ValueError):
        score_threshold_recovery([],{},0)


def test_truth_only_fragment_exclusion_preserves_full_base_accounting():
    full=ArrayRecord('positive',10,210,61,family_id='F')
    fragment=ArrayRecord('partial',20,70,61,family_id='F')
    truth=DiscoveryTruth({'positive':250,'partial':100,'negative':100},[full,fragment],[full],{'partial'},{'F':'A'*61},{})
    predictions=[full,ArrayRecord('partial',0,100,61),ArrayRecord('negative',0,100,61)]
    metrics,details=score_predictions(predictions,truth)
    assert metrics['eligible_read_array_recall']==1
    assert metrics['eligible_read_array_precision']==.5
    assert metrics['eligible_read_negative_read_count']==1
    assert metrics['eligible_read_false_positive_read_count']==1
    assert metrics['all_input_base_union_recall']==1
    assert metrics['all_input_base_union_precision']==250/400
    assert metrics['predictions_on_excluded_reads']==1 and len(details)==2
    with pytest.raises(ValueError):
        score_predictions([ArrayRecord('unknown',0,200,61)],truth)


def test_streamed_fasta_length_validation(tmp_path):
    path=tmp_path/'reads.fa'
    path.write_text('>r1 comment\nACGT\nNN\n>r2\nacgt\n')
    assert fasta_lengths(path)=={'r1':6,'r2':4}
    for value in ('','ACGT\n','>\nACGT\n','>r\n>r\nACGT\n','>r\n','>r\nABCD\n'):
        path.write_text(value)
        with pytest.raises(ValueError):
            fasta_lengths(path)
