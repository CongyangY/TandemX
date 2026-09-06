from dataclasses import asdict
import random

import pytest

from tandemx.discover import mvp
from tandemx.discover.mvp import RepeatFamily, best_local_identity, compare_families, compare_family_pair
from tandemx.utils.kmers import reverse_complement


def test_native_all_offset_identity_matches_python_with_ties_and_ambiguous_bases():
    generator = random.Random(6201)
    pairs = [('AAAA', 'CCCC'), ('NNNN', 'NNNN'), ('ACGACG', 'CGTCGT'), ('acgN', 'ACGn')]
    for _ in range(160):
        a = ''.join(generator.choices('ACGTN', k=generator.choice([1, 2, 11, 49, 50, 51, 70, 120, 171, 300])))
        b = ''.join(generator.choices('ACGTN', k=generator.choice([1, 2, 11, 49, 50, 51, 70, 120, 171, 300])))
        pairs.append((a, b))
    a = ''.join(generator.choices('ACGT', k=171))
    pairs += [(a, a), (a, reverse_complement(a)), (a, a+a), (a, a[25:]+a[:25]),
              (a, a[:50]+'AAAA'+a[54:]), ('A'*120, 'A'*200)]
    for a, b in pairs:
        assert best_local_identity(a, b, 'rust') == best_local_identity(a, b, 'python')
    for backend in ('python', 'rust'):
        with pytest.raises(ValueError, match='nonempty'):
            best_local_identity('', 'ACG', backend)
    with pytest.raises(ValueError, match='ASCII'):
        best_local_identity('ACé', 'ACG', 'rust')


def test_cached_pair_audit_preserves_all_fields_and_builds_each_sketch_once(monkeypatch):
    generator = random.Random(6202)
    a = ''.join(generator.choices('ACGT', k=80))
    sequences = [a, a[:35]+'T'+a[36:], a+a, reverse_complement(a), 'N'*80,
                 ''.join(generator.choices('ACGT', k=61))]
    families = [RepeatFamily(f'f{i}', f'm{i}', s, len(s), 1, len(s), .9, False, 'medium', '')
                for i, s in enumerate(sequences)]
    expected = [asdict(compare_family_pair(a, b, 11)) for i, a in enumerate(families) for b in families[i+1:]]
    original = mvp.canonical_kmer_set
    calls = []
    def counted(sequence, k):
        calls.append(sequence)
        return original(sequence, k)
    monkeypatch.setattr(mvp, 'canonical_kmer_set', counted)
    assert [asdict(r) for r in compare_families(families, 11, 'rust')] == expected
    assert calls == sequences
    assert [asdict(r) for r in compare_families(families, 11, 'python')] == expected


def test_streaming_audit_preserves_table_warning_order_and_optional_collapse(tmp_path):
    from tandemx.discover.family_audit import write_family_audit
    generator = random.Random(6203)
    a = ''.join(generator.choices('ACGT', k=100))
    sequences = [a, a[:44]+'A'+a[45:], a+a, ''.join(generator.choices('ACGT', k=100))]
    families = [RepeatFamily(f'f{i}', f'm{i}', s, len(s), 1, len(s), .9, False, 'medium', 'initial_warning')
                for i, s in enumerate(sequences)]
    expected = compare_families(families, 11, 'python')
    expected_families = mvp.annotate_family_redundancy(families, expected)
    mvp.write_family_similarity(tmp_path/'expected.tsv', expected)
    for keep in (True, False):
        observed, pairs = write_family_audit(tmp_path/'observed.tsv', families, k=11, backend='rust', keep_redundant=keep)
        assert observed == expected_families
        assert (tmp_path/'observed.tsv').read_bytes() == (tmp_path/'expected.tsv').read_bytes()
        assert not (tmp_path/'observed.tsv.partial').exists()
        assert len(pairs) == (sum(p.relationship == 'likely_redundant' for p in expected) if keep else 0)
        if keep:
            assert mvp.collapse_redundant_families(observed, pairs) == mvp.collapse_redundant_families(expected_families, expected)
    assert write_family_audit(tmp_path/'empty.tsv', [], k=11, backend='rust') == ([], [])
    assert (tmp_path/'empty.tsv').read_text() == mvp.family_similarity_header()+'\n'
