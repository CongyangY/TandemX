"""Experimental multi-k word-survival extrapolation, separate from default CLI.

No error rates or truth are inputs. A log-linear fit is an explicit model, not
proof of sequencing accuracy or a calibrated copy-number confidence interval.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
import math

from tandemx.discover.rust_backend import RustDiagnosticKmerCounter
from tandemx.quantify.mvp import MonomerRecord, family_kmer_membership, monomer_kmer_counts
from tandemx.utils.kmers import is_low_complexity_kmer, iter_linear_canonical_kmers


DEFAULT_K_VALUES = (15, 21, 27, 31)


@dataclass(frozen=True)
class AttenuationFit:
    extrapolated_copy_number: float | None
    log_slope_per_base: float | None
    effective_word_loss_probability: float | None
    max_absolute_log_residual: float | None
    leave_one_k_out_log_range: float | None
    status: str


def validate_k_values(k_values: Sequence[int]) -> tuple[int, ...]:
    values = tuple(k_values)
    if (len(values) < 3 or len(values) != len(set(values))
            or any(type(k) is not int or not 1 <= k <= 31 for k in values)):
        raise ValueError('Require at least three distinct integer k values in 1..31')
    return tuple(sorted(values))


def _line(x: Sequence[int], y: Sequence[float]) -> tuple[float, float]:
    mx, my = math.fsum(x)/len(x), math.fsum(y)/len(y)
    slope = math.fsum((a-mx)*(b-my) for a,b in zip(x,y))/math.fsum((a-mx)**2 for a in x)
    return my-slope*mx, slope


def fit_attenuation(k_values: Sequence[int], copies: Sequence[float | None]) -> AttenuationFit:
    """Fit log(theta_k)=a+b*k; no pseudo-counts, slope clamping or CI from k."""
    ordered = validate_k_values(k_values)
    if len(copies) != len(ordered):
        raise ValueError('Every k must have a copy-number estimate')
    pairs = sorted(zip(k_values, copies))
    values = [value for _, value in pairs]
    if any(value is None for value in values):
        return AttenuationFit(None, None, None, None, None, 'missing_diagnostic_or_exposure')
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Copy-number estimates must be finite and nonnegative')
    if any(v == 0 for v in values):
        return AttenuationFit(None, None, None, None, None, 'zero_support_at_one_or_more_k_not_absence')
    logs = [math.log(v) for v in values]
    intercept, slope = _line(ordered, logs)
    residual = max(abs(y-intercept-slope*k) for k,y in zip(ordered,logs))
    omitted = [_line(ordered[:i]+ordered[i+1:], logs[:i]+logs[i+1:])[0] for i in range(len(ordered))]
    try:
        estimate = math.exp(intercept)
        loss = -math.expm1(slope) if slope <= 0 else None
    except OverflowError:
        return AttenuationFit(None, slope, None, residual, max(omitted)-min(omitted), 'numeric_extrapolation_failure')
    if not math.isfinite(estimate) or estimate == 0:
        return AttenuationFit(None, slope, loss, residual, max(omitted)-min(omitted), 'numeric_extrapolation_failure')
    return AttenuationFit(estimate, slope, loss, residual, max(omitted)-min(omitted),
                          'positive_slope_inconsistent_with_loss' if slope > 1e-12 else 'experimental_log_linear_extrapolation')


@dataclass(frozen=True)
class MultiKResult:
    estimates: list[dict]
    per_k: list[dict]
    read_count: int
    total_bases: int
    ambiguous_bases: int


def estimate_multik(sequences: Iterable[str], catalogue: Sequence[MonomerRecord], genome_size: int,
                    k_values: Sequence[int] = DEFAULT_K_VALUES, *, backend: str = 'rust',
                    batch_bases: int = 8_000_000) -> MultiKResult:
    """Stream reads once; count the declared k values on bounded sequence batches.

    theta_k = G/X_k * mean_j(count_j / multiplicity_j), X_k=sum_i max(0,L_i-k+1).
    Targets are family-exclusive and low-complexity-filtered within the supplied
    catalogue only. One unusually long record can exceed the batch-base target.
    """
    ks = validate_k_values(k_values)
    if (type(genome_size) is not int or genome_size <= 0 or type(batch_bases) is not int
            or batch_bases <= 0 or backend not in {'python','rust'}
            or not catalogue or len({m.family_id for m in catalogue}) != len(catalogue)):
        raise ValueError('Require positive genome/batch sizes, supported backend and unique nonempty catalogue')
    if any(not m.family_id or not m.sequence or set(m.sequence.upper())-set('ACGTN') for m in catalogue):
        raise ValueError('Catalogue requires nonempty family IDs and ACGTN sequences')
    diagnostic, counts, native = {}, {}, {}
    for k in ks:
        membership = family_kmer_membership(catalogue, k)
        diagnostic[k] = {m.family_id: {word:n for word,n in monomer_kmer_counts(m.sequence, k).items()
                         if len(membership[word]) == 1 and not is_low_complexity_kmer(word)} for m in catalogue}
        targets = {word for words in diagnostic[k].values() for word in words}
        counts[k] = Counter()
        native[k] = RustDiagnosticKmerCounter(k, targets) if backend == 'rust' else None
    exposure = dict.fromkeys(ks,0)
    valid = dict.fromkeys(ks,0)
    targets_by_k = {k: {w for d in diagnostic[k].values() for w in d} for k in ks}
    reads = bases = ambiguous = buffered = 0
    batch: list[str] = []

    def flush() -> None:
        for k in ks:
            if native[k] is not None:
                native[k].count_sequences(batch)
            else:
                for sequence in batch:
                    counts[k].update(w for w in iter_linear_canonical_kmers(sequence,k) if w in targets_by_k[k])
        batch.clear()

    for sequence in sequences:
        sequence = sequence.upper()
        if not sequence or set(sequence)-set('ACGTN'):
            raise ValueError('Reads must contain nonempty ACGTN sequences')
        reads += 1; bases += len(sequence); ambiguous += sequence.count('N')
        for k in ks:
            exposure[k] += max(0,len(sequence)-k+1)
            valid[k] += (sum(max(0,len(part)-k+1) for part in sequence.split('N'))
                         if 'N' in sequence else max(0,len(sequence)-k+1))
        batch.append(sequence); buffered += len(sequence)
        if buffered >= batch_bases or len(batch) >= 512:
            flush(); buffered = 0
    if not reads or not exposure[max(ks)]:
        raise ValueError('No reads with exposure at the largest k')
    flush()
    for k in ks:
        if native[k] is not None:
            counts[k].update(native[k].counts())
    estimates, per_k = [], []
    for monomer in catalogue:
        values = []
        warnings = ['experimental_model_not_calibrated', 'no_sampling_confidence_interval',
                    'catalogue_specificity_not_background_verified', 'genome_size_treated_as_fixed',
                    'biological_divergence_and_sequencing_error_not_separable']
        if ambiguous:
            warnings.append('ambiguous_base_windows_remain_in_exposure')
        for k in ks:
            words = diagnostic[k][monomer.family_id]
            mean_depth = math.fsum(counts[k][w]/n for w,n in words.items())/len(words) if words else None
            value = genome_size*mean_depth/exposure[k] if mean_depth is not None and exposure[k] else None
            values.append(value)
            per_k.append(dict(family_id=monomer.family_id, k=k, diagnostic_kmer_count=len(words),
                              total_exposure=exposure[k], valid_windows=valid[k], mean_corrected_count=mean_depth,
                              uncorrected_copy_number=value))
        fit = fit_attenuation(ks,values)
        if fit.extrapolated_copy_number is not None and fit.extrapolated_copy_number*len(monomer.sequence) > genome_size:
            warnings.append('estimated_repeat_bases_exceed_genome_size')
        estimates.append(dict(family_id=monomer.family_id, monomer_length=len(monomer.sequence),
                              extrapolated_copy_number=fit.extrapolated_copy_number,
                              uncorrected_mean_k21=values[ks.index(21)] if 21 in ks else None,
                              log_slope_per_base=fit.log_slope_per_base,
                              effective_word_loss_probability=fit.effective_word_loss_probability,
                              max_absolute_log_residual=fit.max_absolute_log_residual,
                              leave_one_k_out_log_range=fit.leave_one_k_out_log_range,
                              status=fit.status, warning=';'.join(warnings)))
    return MultiKResult(estimates, per_k, reads, bases, ambiguous)
