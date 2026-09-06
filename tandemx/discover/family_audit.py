"""Stream the exhaustive representative audit without retaining distinct pairs."""
from __future__ import annotations

from dataclasses import replace
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    from tandemx.discover.mvp import RepeatFamily, FamilySimilarity


def write_family_audit(path: Path, families: Sequence[RepeatFamily], *, k: int,
                       backend: str, keep_redundant: bool = False,
                       logger: logging.Logger | None = None) -> tuple[list[RepeatFamily], list[FamilySimilarity]]:
    from tandemx.discover.mvp import iter_family_similarities, family_similarity_header, format_family_similarity

    logger = logger or logging.getLogger('tandemx.discover')
    pair_count = len(families)*(len(families)-1)//2
    logger.info('family_audit backend=%s families=%s exhaustive_pairs=%s', backend, len(families), pair_count)
    warnings = {family.family_id: [] for family in families}
    redundant = []
    temporary = path.with_suffix(path.suffix+'.partial')
    with temporary.open('w', encoding='utf-8') as handle:
        handle.write(family_similarity_header()+'\n')
        for index, similarity in enumerate(iter_family_similarities(families, k, backend), 1):
            handle.write(format_family_similarity(similarity)+'\n')
            if similarity.relationship != 'distinct':
                warning = f'{similarity.relationship}:{similarity.family_a}-{similarity.family_b}'
                warnings[similarity.family_a].append(warning)
                warnings[similarity.family_b].append(warning)
            if keep_redundant and similarity.relationship == 'likely_redundant':
                redundant.append(similarity)
            if index % 50_000 == 0:
                logger.info('family_audit compared_pairs=%s total_pairs=%s', index, pair_count)
    temporary.replace(path)
    annotated = [replace(family, warning=';'.join(([family.warning] if family.warning else [])+warnings[family.family_id]))
                 for family in families]
    logger.info('family_audit completed_pairs=%s retained_collapse_pairs=%s', pair_count, len(redundant))
    return annotated, redundant
