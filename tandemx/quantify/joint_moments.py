"""Experimental joint-read uncertainty for multi-k log-linear extrapolation.

Only sufficient moments are retained. The first-order sandwich interval treats
reads as independent sampling clusters and genome size as fixed. It does not
calibrate extrapolation bias, background specificity or biological sampling.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

from tandemx.quantify.multik import fit_attenuation, validate_k_values


@dataclass(frozen=True)
class JointEstimate:
    family_id: str
    estimated_copy_number: float | None
    log_sampling_variance: float | None
    sampling_interval_low: float | None
    sampling_interval_high: float | None
    minimum_effective_reads: float
    positive_reads_by_k: tuple[int, ...]
    fit_status: str
    interval_status: str
    warning: str


class JointReadMoments:
    """K-vector first/cross moments, O(families*K²) state independent of reads.

    ``contributions[family][k]`` is the read's mean multiplicity-corrected count
    across that family's diagnostic words at k. Missing families have zero
    contribution. All axes use the supplied increasing k order. Zero-diagnostic
    families remain represented and cannot acquire a sampling interval.
    """
    def __init__(self, k_values: Sequence[int], diagnostic_counts: Mapping[str, Sequence[int]]):
        self.ks = validate_k_values(k_values)
        if tuple(k_values) != self.ks:
            raise ValueError('Joint moment axes require increasing k values')
        self.size = len(self.ks)
        if not diagnostic_counts or any(not isinstance(name, str) or not name for name in diagnostic_counts):
            raise ValueError('Require nonempty diagnostic family IDs')
        self.diagnostic_counts = {f: tuple(v) for f, v in diagnostic_counts.items()}
        if any(len(v) != self.size or any(type(x) is not int or x < 0 for x in v)
               for v in self.diagnostic_counts.values()):
            raise ValueError('Diagnostic counts must be nonnegative integers on every k axis')
        self.read_count = self.total_bases = 0
        self.sum_x = [0] * self.size
        self.sum_xx = [0] * (self.size*self.size)
        self.sum_y = {f: [0.0]*self.size for f in diagnostic_counts}
        self.sum_yy = {f: [0.0]*(self.size*self.size) for f in diagnostic_counts}
        self.sum_yx = {f: [0.0]*(self.size*self.size) for f in diagnostic_counts}
        self.positive = {f: [0]*self.size for f in diagnostic_counts}

    def add(self, length: int, contributions: Mapping[str, Sequence[float]]) -> None:
        if type(length) is not int or length <= 0 or contributions.keys()-self.sum_y.keys():
            raise ValueError('Require positive read length and known families')
        # Validate before mutating any moment, including implicit zero axes.
        for family, values in contributions.items():
            if (len(values) != self.size
                    or any(not isinstance(v, (int, float)) or isinstance(v, bool)
                           or not math.isfinite(v) or v < 0 for v in values)
                    or any(v and not n for v, n in zip(values, self.diagnostic_counts[family]))):
                raise ValueError('Invalid joint read contribution or no diagnostic target')
            if any(value > max(0, length-k+1) for k, value in zip(self.ks, values)):
                raise ValueError('Read contribution exceeds its available k-mer windows')
        x = [max(0, length-k+1) for k in self.ks]
        self.read_count += 1
        self.total_bases += length
        for i, xi in enumerate(x):
            self.sum_x[i] += xi
            for j, xj in enumerate(x):
                self.sum_xx[i*self.size+j] += xi*xj
        for family, values in contributions.items():
            y, yy, yx, positive = (self.sum_y[family], self.sum_yy[family],
                                   self.sum_yx[family], self.positive[family])
            for i, value in enumerate(values):
                y[i] += value
                positive[i] += value > 0
                for j in range(self.size):
                    yy[i*self.size+j] += value*values[j]
                    yx[i*self.size+j] += value*x[j]

    def estimates(self, genome_size: int, minimum_effective_reads: int = 20) -> list[JointEstimate]:
        if type(genome_size) is not int or genome_size <= 0 or type(minimum_effective_reads) is not int or minimum_effective_reads < 2:
            raise ValueError('Require positive integer genome size and at least two effective reads')
        if not self.read_count or not self.sum_x[-1]:
            raise ValueError('No reads with exposure at the largest k')
        mean_k = math.fsum(self.ks)/self.size
        ss_k = math.fsum((k-mean_k)**2 for k in self.ks)
        # OLS intercept is this fixed linear combination of log(theta_k).
        weights = [1/self.size-mean_k*(k-mean_k)/ss_k for k in self.ks]
        results = []
        for family, y in self.sum_y.items():
            yy, yx = self.sum_yy[family], self.sum_yx[family]
            copies = [genome_size*value/exposure if diagnostic and exposure else None
                      for value, exposure, diagnostic in zip(y, self.sum_x, self.diagnostic_counts[family])]
            fit = fit_attenuation(self.ks, copies)
            effective = min((value**2/yy[i*self.size+i] if yy[i*self.size+i] else 0.0)
                            for i, value in enumerate(y))
            variance = low = high = None
            status = 'no_finite_extrapolation'
            warning = ['experimental_joint_read_delta_method', 'k_values_are_correlated_not_replicates',
                       'requires_independent_reads', 'model_bias_not_calibrated',
                       'catalogue_specificity_not_background_verified', 'genome_size_treated_as_fixed']
            if fit.extrapolated_copy_number is not None:
                if self.read_count < 2 or effective < minimum_effective_reads:
                    status = 'insufficient_effective_reads'
                else:
                    a = [w/value for w, value in zip(weights, y)]
                    b = [w/exposure for w, exposure in zip(weights, self.sum_x)]
                    terms = []
                    for i in range(self.size):
                        for j in range(self.size):
                            index = i*self.size+j
                            terms.extend((a[i]*a[j]*yy[index], b[i]*b[j]*self.sum_xx[index],
                                          -2*a[i]*b[j]*yx[index]))
                    residual = math.fsum(terms)
                    resolution = 1e-12*math.fsum(abs(t) for t in terms)
                    if residual < -resolution:
                        raise ArithmeticError('Invalid joint read influence sum of squares')
                    if residual <= resolution:
                        status = 'variance_below_moment_resolution_uncalibrated'
                    else:
                        variance = residual*self.read_count/(self.read_count-1)
                        half = 1.959963984540054*math.sqrt(variance)
                        log_estimate = math.log(fit.extrapolated_copy_number)
                        try:
                            low, high = math.exp(log_estimate-half), math.exp(log_estimate+half)
                            status = 'approximate_95pct_joint_read_log_normal'
                        except OverflowError:
                            low = high = None
                            status = 'numeric_interval_failure'
                        if low == 0 or high is not None and not math.isfinite(high):
                            low = high = None
                            status = 'numeric_interval_failure'
            results.append(JointEstimate(family, fit.extrapolated_copy_number, variance, low, high,
                                          effective, tuple(self.positive[family]), fit.status,
                                          status, ';'.join(warning)))
        return results
