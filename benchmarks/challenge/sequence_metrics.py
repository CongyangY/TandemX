"""Independent, exhaustive cyclic edit-distance endpoints using edlib.

No TandemX alignment or clustering implementation is imported here. Every
rotation of the predicted monomer and its reverse complement is compared by
global unit-cost Levenshtein distance. Ns are never treated as matching evidence.
"""
from __future__ import annotations

from functools import lru_cache
import math
from statistics import fmean

from .evaluate import canonical_monomer, maximum_matching


@lru_cache(maxsize=8192)
def cyclic_edit_similarity(truth: str, prediction: str) -> float:
    if not truth or not prediction:
        return 0.0
    if set((truth + prediction).upper()) - set("ACGTN"):
        raise ValueError("Monomers must contain only ACGTN")
    try:
        import edlib
    except ImportError as exc:
        raise RuntimeError("Install the independent evaluator: pip install -e '.[benchmark]'") from exc
    truth, prediction = truth.upper(), prediction.upper()
    # Different placeholders prevent N==N from counting as known sequence.
    reference = truth.replace("N", "X")
    reverse = prediction.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    lower_bound = abs(len(truth) - len(prediction))
    best = max(len(truth), len(prediction))
    for strand in (prediction, reverse):
        strand = strand.replace("N", "Y")
        for offset in range(len(strand)):
            rotated = strand[offset:] + strand[:offset]
            distance = edlib.align(reference, rotated, mode="NW", task="distance", k=best)["editDistance"]
            if distance >= 0:
                best = min(best, distance)
            if best == lower_bound:
                return 1 - best / max(len(truth), len(prediction))
    return 1 - best / max(len(truth), len(prediction))


def score_cyclic_recovery(predictions: list[str], truth: dict[str, str], threshold: float = 0.9
                          ) -> tuple[dict, list[dict]]:
    """One-to-one planted-monomer recovery plus distinct-consensus purity.

    Homologous-consensus fraction allows several allelic/error variants to match
    one truth. It is deliberately not named family precision. One-to-one recall
    still prevents a single merged consensus from recovering two planted units.
    """
    if not 0 < threshold <= 1:
        raise ValueError("Edit-similarity threshold must be in (0,1]")
    sequences = sorted(set(canonical_monomer(s) for s in predictions if s))
    identifiers = sorted(truth)
    scores = [[cyclic_edit_similarity(truth[name], seq) for seq in sequences] for name in identifiers]
    edges = [[j for j, score in enumerate(row) if score + 1e-12 >= threshold] for row in scores]
    matching = maximum_matching(edges)
    supported = {j for row in edges for j in row}
    rows = [{"truth_id": name, "best_cyclic_edit_similarity": max(scores[i], default=0.0),
             "recovered": int(i in matching), "assigned_sequence_index": matching.get(i, "NA"),
             "criterion": f"cyclic_global_levenshtein_similarity_ge_{threshold:g}"}
            for i, name in enumerate(identifiers)]
    metrics = {"cyclic_monomer_recall": len(matching) / len(truth) if truth else math.nan,
               "distinct_consensus_count": len(sequences),
               "homologous_consensus_fraction": len(supported) / len(sequences) if sequences else math.nan,
               "unmatched_distinct_consensus_count": len(sequences) - len(supported),
               "mean_best_cyclic_edit_similarity": fmean(max(row, default=0.0) for row in scores) if truth else math.nan}
    return metrics, rows


@lru_cache(maxsize=8192)
def cyclic_reaches_threshold(truth: str, prediction: str, threshold: float = .9) -> bool:
    """Exact threshold decision with bounded edlib, not an approximate score.

    Avoid calculating below-threshold exact identities when only recovery is
    requested. Length bounds and early success preserve the original decision.
    """
    if not 0<threshold<=1 or not truth or not prediction or set((truth+prediction).upper())-set('ACGTN'):
        raise ValueError('Require nonempty ACGTN sequences and threshold in (0,1]')
    maximum=max(len(truth),len(prediction))
    limit=math.floor((1-threshold+1e-12)*maximum)
    if abs(len(truth)-len(prediction))>limit:
        return False
    import edlib
    reference=truth.upper().replace('N','X')
    prediction=prediction.upper()
    reverse=prediction.translate(str.maketrans('ACGT','TGCA'))[::-1]
    for strand in (prediction,reverse):
        strand=strand.replace('N','Y')
        for offset in range(len(strand)):
            if edlib.align(reference,strand[offset:]+strand[:offset],mode='NW',task='distance',k=limit)['editDistance']>=0:
                return True
    return False


def score_threshold_recovery(predictions: list[str], truth: dict[str,str], threshold: float=.9
                             ) -> tuple[dict,list[dict]]:
    if not 0<threshold<=1:
        raise ValueError('Threshold must be in (0,1]')
    sequences=sorted(set(canonical_monomer(s) for s in predictions if s))
    identifiers=sorted(truth)
    edges=[[j for j,sequence in enumerate(sequences) if cyclic_reaches_threshold(truth[name],sequence,threshold)]
           for name in identifiers]
    matching=maximum_matching(edges)
    supported={j for row in edges for j in row}
    rows=[dict(truth_id=name,recovered=int(i in matching),assigned_sequence_index=matching.get(i,'NA'),
               threshold=threshold,criterion='exact_cyclic_levenshtein_threshold_no_below_threshold_score')
          for i,name in enumerate(identifiers)]
    return dict(truth_family_count=len(truth),recovered_family_count=len(matching),
                cyclic_monomer_recall=len(matching)/len(truth) if truth else math.nan,
                distinct_consensus_count=len(sequences),
                homologous_consensus_fraction=len(supported)/len(sequences) if sequences else math.nan,
                unmatched_distinct_consensus_count=len(sequences)-len(supported)),rows
