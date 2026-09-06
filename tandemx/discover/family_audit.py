"""Stream the exhaustive representative audit without retaining distinct pairs."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    from tandemx.discover.mvp import RepeatFamily, FamilySimilarity


def related_pair_candidates(families: Sequence[RepeatFamily], k: int, backend: str):
    """Exact upper-bound filter for the current redundancy rules, not a sketch.

    Every non-distinct rule requires shared canonical k-mers. An inverted index
    counts exact intersections; identity/overlap set to one is an upper bound.
    Only one representative's candidate-counter is retained at a time.
    """
    from tandemx.discover.mvp import canonical_kmer_set, classify_family_relationship, compare_family_pair

    kmers = [canonical_kmer_set(f.monomer_sequence, k) for f in families]
    postings = defaultdict(list)
    for index, words in enumerate(kmers):
        for word in words:
            postings[word].append(index)
    for index, words in enumerate(kmers):
        shared_counts = Counter(other for word in words for other in postings[word] if other > index)
        for other, shared in sorted(shared_counts.items()):
            a, b = families[index], families[other]
            relation, _, _ = classify_family_relationship(
                kmer_jaccard=shared/(len(words)+len(kmers[other])-shared),
                shared_kmer_fraction=shared/min(len(words), len(kmers[other])),
                local_identity=1.0, local_overlap_fraction_shorter=1.0,
                length_ratio=max(a.monomer_length_bp, b.monomer_length_bp)/min(a.monomer_length_bp, b.monomer_length_bp))
            if relation != 'distinct':
                yield compare_family_pair(a, b, k, backend, _kmers=(words, kmers[other]))


def write_family_audit(path: Path, families: Sequence[RepeatFamily], *, k: int,
                       backend: str, keep_redundant: bool = False,
                       logger: logging.Logger | None = None, mode: str = 'full') -> tuple[list[RepeatFamily], list[FamilySimilarity]]:
    from tandemx.discover.mvp import iter_family_similarities, family_similarity_header, format_family_similarity

    logger = logger or logging.getLogger('tandemx.discover')
    if mode not in {'full', 'related'} or backend not in {'python', 'rust'} or k < 1:
        raise ValueError('Family audit requires full/related mode, python/rust backend and positive k')
    if len({f.family_id for f in families}) != len(families):
        raise ValueError('Family audit requires unique family identifiers')
    pair_count = len(families)*(len(families)-1)//2
    logger.info('family_audit mode=%s backend=%s families=%s possible_pairs=%s', mode, backend, len(families), pair_count)
    warnings = {family.family_id: [] for family in families}
    redundant = []
    scored = emitted = related = 0
    similarities = (iter_family_similarities(families, k, backend) if mode == 'full'
                    else related_pair_candidates(families, k, backend))
    temporary = path.with_suffix(path.suffix+'.partial')
    with temporary.open('w', encoding='utf-8') as handle:
        handle.write(family_similarity_header()+'\n')
        for index, similarity in enumerate(similarities, 1):
            scored += 1
            if mode == 'full' or similarity.relationship != 'distinct':
                handle.write(format_family_similarity(similarity)+'\n')
                emitted += 1
            if similarity.relationship != 'distinct':
                related += 1
                warning = f'{similarity.relationship}:{similarity.family_a}-{similarity.family_b}'
                warnings[similarity.family_a].append(warning)
                warnings[similarity.family_b].append(warning)
            if keep_redundant and similarity.relationship == 'likely_redundant':
                redundant.append(similarity)
            if index % 50_000 == 0:
                logger.info('family_audit compared_pairs=%s total_pairs=%s', index, pair_count)
    temporary.replace(path)
    receipt = dict(schema_version=1, complete=True, mode=mode, backend=backend,
                   family_count=len(families), possible_pairs=pair_count, pairs_scored=scored,
                   pairs_pruned_by_kmer_gate=pair_count-scored, emitted_pairs=emitted,
                   related_pairs=related, omitted_distinct_pairs=pair_count-emitted,
                   warning='heuristic_relationships_not_biological_truth;dense_catalogues_can_remain_quadratic')
    summary = path.with_name('family_audit_summary.json')
    partial_summary = summary.with_suffix('.json.partial')
    partial_summary.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    partial_summary.replace(summary)
    annotated = [replace(family, warning=';'.join(([family.warning] if family.warning else [])+warnings[family.family_id]))
                 for family in families]
    logger.info('family_audit scored_pairs=%s emitted_pairs=%s related_pairs=%s retained_collapse_pairs=%s', scored, emitted, related, len(redundant))
    return annotated, redundant
