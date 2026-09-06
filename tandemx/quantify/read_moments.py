"""Experimental read-cluster ratio estimator; not the default quantify model.

Keeps family-level moments, never all reads or a read-by-family matrix. The
normal intervals describe sampling under independent reads, not sequence-error,
library-selection, background-specificity or genome-size uncertainty.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import math

from tandemx.utils.kmers import canonical_kmer, canonical_kmer_code, iter_canonical_kmer_codes


@dataclass
class FamilyMoments:
    sum_y: float = 0.0
    sum_y2: float = 0.0
    sum_xy: float = 0.0
    positive_reads: int = 0

    def add(self, exposure: int, contribution: float) -> None:
        if exposure < 0 or contribution < 0 or not math.isfinite(contribution):
            raise ValueError('Read exposure and contribution must be finite and nonnegative')
        self.sum_y += contribution
        self.sum_y2 += contribution * contribution
        self.sum_xy += exposure * contribution
        self.positive_reads += contribution > 0


@dataclass(frozen=True)
class ReadClusterEstimate:
    family_id: str
    diagnostic_kmer_count: int
    estimated_copy_number: float | None
    sampling_standard_error: float | None
    sampling_interval_low: float | None
    sampling_interval_high: float | None
    positive_reads: int
    effective_positive_reads: float
    interval_status: str
    warning: str


@dataclass
class ReadMoments:
    """Raw moments; zero family contributions are implicit for each read."""
    families: dict[str, FamilyMoments] = field(default_factory=dict)
    read_count: int = 0
    total_bases: int = 0
    sum_x: int = 0
    sum_x2: int = 0
    valid_windows: int = 0
    short_reads: int = 0

    def add(self, length: int, k: int, contributions: Mapping[str, float], valid_windows: int) -> None:
        if length <= 0 or k <= 0:
            raise ValueError('Read length and k must be positive')
        exposure = max(0, length-k+1)
        if not 0 <= valid_windows <= exposure or set(contributions)-self.families.keys():
            raise ValueError('Invalid windows or unknown family contribution')
        self.read_count += 1
        self.total_bases += length
        self.sum_x += exposure
        self.sum_x2 += exposure * exposure
        self.valid_windows += valid_windows
        self.short_reads += exposure == 0
        for family, contribution in contributions.items():
            self.families[family].add(exposure, contribution)

    def estimates(self, genome_size: int, diagnostic_counts: Mapping[str, int],
                  minimum_effective_reads: int = 20) -> list[ReadClusterEstimate]:
        if genome_size <= 0 or minimum_effective_reads < 2:
            raise ValueError('Positive genome size and at least two effective reads required')
        if self.read_count == 0 or self.sum_x == 0:
            raise ValueError('No usable k-mer exposure in input reads')
        if set(diagnostic_counts) != self.families.keys() or any(v < 0 for v in diagnostic_counts.values()):
            raise ValueError('Diagnostic counts must match every family')
        result = []
        for family, moment in self.families.items():
            warnings = ['experimental_read_sampling_model', 'catalogue_specificity_not_background_verified',
                        'sequence_error_and_library_bias_not_corrected', 'genome_size_treated_as_fixed']
            if self.sum_x != self.valid_windows:
                warnings.append('ambiguous_base_windows_in_exposure')
            effective = moment.sum_y**2/moment.sum_y2 if moment.sum_y2 else 0.0
            estimate = se = low = high = None
            if diagnostic_counts[family] == 0:
                status = 'no_diagnostic_kmers'
            else:
                ratio = moment.sum_y/self.sum_x
                estimate = genome_size*ratio
                if moment.positive_reads == 0:
                    status = 'no_observed_support_not_absence'
                elif self.read_count < 2 or effective < minimum_effective_reads:
                    status = 'insufficient_effective_reads'
                else:
                    terms = (moment.sum_y2, -2*ratio*moment.sum_xy, ratio*ratio*self.sum_x2)
                    residual_ss = math.fsum(terms)
                    # Detect material cancellation/invalid state; only roundoff is clamped.
                    if residual_ss < -1e-12*sum(abs(t) for t in terms):
                        raise ArithmeticError('Invalid read residual sum of squares')
                    se = genome_size*math.sqrt(max(0.0, residual_ss)*self.read_count/(self.read_count-1))/self.sum_x
                    if se == 0:
                        status = 'zero_observed_variance_uncalibrated'
                        se = None
                    else:
                        low, high = max(0.0, estimate-1.959963984540054*se), estimate+1.959963984540054*se
                        status = 'approximate_95pct_read_cluster_normal'
                        warnings.append('normal_coverage_requires_independent_calibration')
            result.append(ReadClusterEstimate(family, diagnostic_counts[family], estimate, se, low, high,
                                               moment.positive_reads, effective, status, ';'.join(warnings)))
        return result


def collect_read_moments(sequences: Iterable[str], diagnostic: Mapping[str, Mapping[str, int]],
                         k: int) -> ReadMoments:
    """Y_if is mean multiplicity-corrected target count in each individual read.

    X_i is max(0,L_i-k+1), including N-containing opportunities. They are missing
    evidence and explicitly flagged, not silently removed to inflate depth.
    Targets must already be canonical, family-exclusive, catalogue-filtered.
    """
    if not 1 <= k <= 31 or not diagnostic:
        raise ValueError('Provide diagnostic families and k in 1..31')
    targets: dict[int, tuple[str, float]] = {}
    for family, kmers in diagnostic.items():
        if not family:
            raise ValueError('Empty family identifier')
        for kmer, multiplicity in kmers.items():
            if (len(kmer) != k or set(kmer)-set('ACGT') or canonical_kmer(kmer) != kmer
                    or not isinstance(multiplicity, int) or multiplicity <= 0):
                raise ValueError('Diagnostic targets require canonical ACGT and positive integer multiplicity')
            code = canonical_kmer_code(kmer)
            if code in targets:
                raise ValueError('A diagnostic target cannot belong to multiple families')
            targets[code] = (family, 1/(len(kmers)*multiplicity))
    moments = ReadMoments(families={family: FamilyMoments() for family in diagnostic})
    for sequence in sequences:
        if not sequence or set(sequence.upper())-set('ACGTN'):
            raise ValueError('Reads must contain nonempty ACGTN sequences')
        contributions: dict[str, float] = defaultdict(float)
        valid = 0
        for _, code in iter_canonical_kmer_codes(sequence, k):
            valid += 1
            target = targets.get(code)
            if target is not None:
                family, weight = target
                contributions[family] += weight
        moments.add(len(sequence), k, contributions, valid)
    return moments
