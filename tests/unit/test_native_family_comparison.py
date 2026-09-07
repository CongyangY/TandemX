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
    from tandemx.discover.hierarchy import write_family_hierarchy
    generator = random.Random(6203)
    a = ''.join(generator.choices('ACGT', k=100))
    sequences = [a, a[:44]+'A'+a[45:], a+a, ''.join(generator.choices('ACGT', k=100))]
    families = [RepeatFamily(f'f{i}', f'm{i}', s, len(s), 1, len(s), .9, False, 'medium', 'initial_warning')
                for i, s in enumerate(sequences)]
    expected = compare_families(families, 11, 'python')
    expected_families = mvp.annotate_family_redundancy(families, expected)
    mvp.write_family_similarity(tmp_path/'expected.tsv', expected)
    write_family_hierarchy(tmp_path/'expected_hierarchy.tsv', expected)
    for keep in (True, False):
        observed, pairs = write_family_audit(tmp_path/'observed.tsv', families, k=11, backend='rust', keep_redundant=keep)
        assert observed == expected_families
        assert (tmp_path/'observed.tsv').read_bytes() == (tmp_path/'expected.tsv').read_bytes()
        assert (tmp_path/'family_hierarchy.tsv').read_bytes() == (tmp_path/'expected_hierarchy.tsv').read_bytes()
        assert not (tmp_path/'observed.tsv.partial').exists()
        assert not (tmp_path/'family_hierarchy.tsv.partial').exists()
        assert len(pairs) == (sum(p.relationship == 'likely_redundant' for p in expected) if keep else 0)
        if keep:
            assert mvp.collapse_redundant_families(observed, pairs) == mvp.collapse_redundant_families(expected_families, expected)
    assert write_family_audit(tmp_path/'empty.tsv', [], k=11, backend='rust') == ([], [])
    assert (tmp_path/'empty.tsv').read_text() == mvp.family_similarity_header()+'\n'
    assert (tmp_path/'family_hierarchy.tsv').read_text().count('\n') == 1


@pytest.mark.parametrize('k', [1, 11, 31, 200])
def test_related_index_preserves_every_non_distinct_pair_and_collapse(tmp_path, k):
    import json
    from tandemx.discover.family_audit import write_family_audit
    rng = random.Random(6204)
    sequences = []
    for _ in range(8):
        a = ''.join(rng.choices('ACGT', k=rng.choice([31, 50, 71, 100])))
        sequences += [a, a+a, reverse_complement(a), a[:20]+'A'+a[21:], a[:len(a)//2]]
    sequences += ['N'*80, 'N'*100, 'A'*120, 'A'*160]
    families = [RepeatFamily(f'f{i}', f'm{i}', s, len(s), 1, len(s), .9, False, 'medium', 'initial')
                for i, s in enumerate(sequences)]
    full = compare_families(families, k, 'rust')
    expected = [p for p in full if p.relationship != 'distinct']
    mvp.write_family_similarity(tmp_path/'expected.tsv', expected)
    annotated, collapse = write_family_audit(tmp_path/'related.tsv', families, k=k, backend='rust', mode='related', keep_redundant=True)
    assert (tmp_path/'related.tsv').read_bytes() == (tmp_path/'expected.tsv').read_bytes()
    assert annotated == mvp.annotate_family_redundancy(families, full)
    assert mvp.collapse_redundant_families(annotated, collapse) == mvp.collapse_redundant_families(annotated, full)
    summary = json.loads((tmp_path/'family_audit_summary.json').read_text())
    assert summary['related_pairs'] == summary['emitted_pairs'] == len(expected)
    assert summary['hierarchy_edges'] == sum(
        pair.relationship == 'possible_higher_order_or_partial' for pair in expected
    )
    assert summary['putative_period_multiple_edges'] <= summary['hierarchy_edges']
    assert summary['omitted_distinct_pairs'] == len(full)-len(expected)
    assert summary['possible_pairs'] == summary['pairs_scored']+summary['pairs_pruned_by_kmer_gate']
    if k == 200:
        assert summary['pairs_scored'] == 0
    with pytest.raises(ValueError, match='full/related'):
        write_family_audit(tmp_path/'invalid.tsv', [], k=k, backend='rust', mode='skip')


@pytest.mark.parametrize('k', [1, 2, 11, 31, 200])
def test_compact_audit_tokens_match_string_canonical_kmers(k):
    from tandemx.discover.family_audit import canonical_kmer_tokens
    rng = random.Random(953107 + k)
    sequences = [
        ''.join(rng.choices('ACGTN', k=size))
        for size in (1, 2, 10, 31, 211, 400)
    ]
    for sequence in sequences:
        strings = mvp.canonical_kmer_set(sequence, k)
        encoded = canonical_kmer_tokens(sequence, k)
        assert len(encoded) == len(strings)
        for word in strings:
            forward = reverse = 0
            for base in word:
                forward = (forward << 3) | {"A": 0, "C": 1, "G": 2, "N": 3, "T": 4}[base]
            for base in word.translate(str.maketrans("ACGT", "TGCA"))[::-1]:
                reverse = (reverse << 3) | {"A": 0, "C": 1, "G": 2, "N": 3, "T": 4}[base]
            assert min(forward, reverse) in encoded
    with pytest.raises(ValueError, match='ACGTN'):
        canonical_kmer_tokens('ACG?', min(k, 4))
