"""Bounded read collector for experimental joint multi-k sampling uncertainty."""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from tandemx.discover.rust_backend import RustBackendUnavailable
from tandemx.quantify.joint_moments import JointEstimate, JointReadMoments
from tandemx.quantify.multik import DEFAULT_K_VALUES, validate_k_values
from tandemx.quantify.mvp import MonomerRecord, family_kmer_membership, monomer_kmer_counts
from tandemx.utils.kmers import (canonical_kmer_code, is_low_complexity_kmer,
                                  iter_canonical_kmer_codes)


@dataclass(frozen=True)
class JointMultiKResult:
    estimates: list[JointEstimate]
    per_k: list[dict]
    read_count: int
    total_bases: int
    ambiguous_bases: int
    warning: str


def estimate_joint_multik(
    sequences: Iterable[str], catalogue: Sequence[MonomerRecord], genome_size: int,
    k_values: Sequence[int] = DEFAULT_K_VALUES, *, backend: str = 'rust',
    batch_bases: int = 8_000_000, minimum_effective_reads: int = 20,
) -> JointMultiKResult:
    """Use each read as one sampling cluster, retaining all its k axes together.

    The point estimator is mathematically identical to ``estimate_multik``.
    Per-word weights are 1/(diagnostic count * circular multiplicity). Memory is
    O(target words + families*K² + bounded batch); one long read can exceed the
    batch-base target. Reads themselves are never retained after a batch.
    """
    ks = validate_k_values(k_values)
    if (type(genome_size) is not int or genome_size <= 0 or type(batch_bases) is not int
            or batch_bases <= 0 or backend not in {'python', 'rust'}
            or type(minimum_effective_reads) is not int or minimum_effective_reads < 2
            or not catalogue or len({m.family_id for m in catalogue}) != len(catalogue)):
        raise ValueError('Require positive sizes, supported backend, unique catalogue and at least two effective reads')
    if any(not isinstance(m.family_id, str) or not m.family_id or not m.sequence
           or set(m.sequence.upper())-set('ACGTN') for m in catalogue):
        raise ValueError('Catalogue requires nonempty family IDs and ACGTN sequences')
    names = [m.family_id for m in catalogue]
    diagnostic, banks, counters = {}, [], []
    for k in ks:
        membership = family_kmer_membership(catalogue, k)
        diagnostic[k] = {m.family_id: {word: n for word, n in monomer_kmer_counts(m.sequence, k).items()
                         if len(membership[word]) == 1 and not is_low_complexity_kmer(word)} for m in catalogue}
        targets = [(word, i, 1/(len(words)*n)) for i, name in enumerate(names)
                   for words in [diagnostic[k][name]] for word, n in words.items()]
        if backend == 'rust':
            try:
                from tandemx._rust_core import WeightedKmerCounter
            except ImportError as exc:
                raise RustBackendUnavailable('Rebuild the Rust extension for joint read weighted counts') from exc
            counters.append(WeightedKmerCounter(k, len(names), targets))
        else:
            banks.append({canonical_kmer_code(w): (i, weight) for w, i, weight in targets})
    moments = JointReadMoments(ks, {name: tuple(len(diagnostic[k][name]) for k in ks) for name in names})
    valid = [0]*len(ks)
    ambiguous = buffered = 0
    batch: list[str] = []

    def flush() -> None:
        if not batch:
            return
        by_k = []
        for axis, k in enumerate(ks):
            if backend == 'rust':
                by_k.append(counters[axis].count_sequences(batch))
            else:
                rows = []
                for sequence in batch:
                    values = defaultdict(float)
                    for _, code in iter_canonical_kmer_codes(sequence, k):
                        target = banks[axis].get(code)
                        if target is not None:
                            i, weight = target
                            values[i] += weight
                    rows.append(sorted(values.items()))
                by_k.append(rows)
        for read_index, sequence in enumerate(batch):
            joint: dict[str, list[float]] = {}
            for axis, rows in enumerate(by_k):
                for family_index, value in rows[read_index]:
                    name = names[family_index]
                    joint.setdefault(name, [0.0]*len(ks))[axis] = value
            moments.add(len(sequence), joint)
        batch.clear()

    for sequence in sequences:
        sequence = sequence.upper()
        if not sequence or set(sequence)-set('ACGTN'):
            raise ValueError('Reads require nonempty ACGTN sequences')
        ambiguous += sequence.count('N')
        for axis, k in enumerate(ks):
            valid[axis] += sum(max(0, len(part)-k+1) for part in sequence.split('N'))
        batch.append(sequence)
        buffered += len(sequence)
        if buffered >= batch_bases or len(batch) >= 512:
            flush()
            buffered = 0
    flush()
    estimates = moments.estimates(genome_size, minimum_effective_reads)
    points = [dict(family_id=name, k=k, diagnostic_kmer_count=len(diagnostic[k][name]),
                   total_exposure=moments.sum_x[axis], valid_windows=valid[axis],
                   mean_corrected_count=moments.sum_y[name][axis] if diagnostic[k][name] else None,
                   uncorrected_copy_number=(genome_size*moments.sum_y[name][axis]/moments.sum_x[axis]
                       if diagnostic[k][name] and moments.sum_x[axis] else None))
              for name in names for axis, k in enumerate(ks)]
    warning = 'experimental_sampling_interval_not_total_uncertainty;biological_divergence_and_sequencing_error_not_separable'
    if ambiguous:
        warning += ';ambiguous_base_windows_remain_in_exposure'
    return JointMultiKResult(estimates, points, moments.read_count, moments.total_bases, ambiguous, warning)
