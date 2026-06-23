use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};

const BIN_SIZE: usize = 5;
const REFINEMENT_RADIUS: isize = 2;
const MAX_SEED_GROUPS: usize = 128;
const MAX_IDENTITY_COMPARISONS: usize = 1024;
const ACCEPTANCE_SCORE: f64 = 0.75;
const SEED_WEIGHT: f64 = 0.20;
const IDENTITY_WEIGHT: f64 = 0.80;

#[pyclass(frozen)]
struct ScanResult {
    #[pyo3(get)]
    candidate_periods: Vec<usize>,
    #[pyo3(get)]
    spacing_support: Vec<(usize, usize)>,
    #[pyo3(get)]
    best_period: usize,
    #[pyo3(get)]
    periodicity_score: f64,
    #[pyo3(get)]
    overflow_count: usize,
    #[pyo3(get)]
    status: String,
}

#[pyclass]
struct DiagnosticKmerCounter {
    k: usize,
    target_indices: HashMap<u64, usize>,
    target_names: Vec<String>,
    counts: Vec<u64>,
}

fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn canonical_code(sequence: &[u8]) -> Option<u64> {
    let k = sequence.len();
    if k == 0 || k > 31 {
        return None;
    }
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    for (index, &base) in sequence.iter().enumerate() {
        let code = base_code(base)?;
        forward = (forward << 2) | code;
        reverse |= (3 - code) << (2 * index);
    }
    Some(forward.min(reverse))
}

#[pymethods]
impl DiagnosticKmerCounter {
    #[new]
    fn new(k: usize, targets: Vec<String>) -> PyResult<Self> {
        if !(1..=31).contains(&k) {
            return Err(PyValueError::new_err(
                "Rust k-mer counter requires k in 1..=31",
            ));
        }
        let mut target_indices = HashMap::new();
        let mut target_names = Vec::new();
        for target in targets {
            if target.len() != k {
                return Err(PyValueError::new_err(
                    "diagnostic k-mer length does not match k",
                ));
            }
            let code = canonical_code(target.as_bytes()).ok_or_else(|| {
                PyValueError::new_err("diagnostic k-mer contains an invalid base")
            })?;
            if !target_indices.contains_key(&code) {
                let index = target_indices.len();
                target_indices.insert(code, index);
                target_names.push(target);
            }
        }
        let counts = vec![0_u64; target_indices.len()];
        Ok(Self {
            k,
            target_indices,
            target_names,
            counts,
        })
    }

    fn count_sequence(&mut self, sequence: &str) {
        let mask = (1_u64 << (2 * self.k)) - 1;
        let reverse_shift = 2 * (self.k - 1);
        let mut forward = 0_u64;
        let mut reverse = 0_u64;
        let mut valid_length = 0_usize;
        for &base in sequence.as_bytes() {
            let Some(code) = base_code(base) else {
                forward = 0;
                reverse = 0;
                valid_length = 0;
                continue;
            };
            valid_length += 1;
            forward = ((forward << 2) | code) & mask;
            reverse = (reverse >> 2) | ((3 - code) << reverse_shift);
            if valid_length < self.k {
                continue;
            }
            let canonical = forward.min(reverse);
            if let Some(&index) = self.target_indices.get(&canonical) {
                self.counts[index] += 1;
            }
        }
    }

    fn counts(&self) -> HashMap<String, u64> {
        self.target_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), self.counts[index]))
            .collect()
    }
}

fn extract_repeated_positions(
    sequence: &[u8],
    k: usize,
    min_seed_occurrences: usize,
    max_pairs_per_kmer: usize,
) -> (HashMap<u64, Vec<usize>>, usize) {
    let mut first_positions: HashMap<u64, usize> = HashMap::new();
    let mut repeated_positions: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut overflowed: HashSet<u64> = HashSet::new();
    let position_cap = max_pairs_per_kmer + 1;
    let mask = (1_u64 << (2 * k)) - 1;
    let reverse_shift = 2 * (k - 1);
    let complexity_threshold = (4 * k).div_ceil(5);
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mut valid_length = 0_usize;
    let mut base_counts = [0_usize; 4];
    let mut distinct_bases = 0_usize;
    let mut window_codes = vec![0_usize; k];

    for (index, &base) in sequence.iter().enumerate() {
        let Some(code_u64) = base_code(base) else {
            forward = 0;
            reverse = 0;
            valid_length = 0;
            base_counts = [0; 4];
            distinct_bases = 0;
            continue;
        };
        let code = code_u64 as usize;
        let slot = index % k;
        if valid_length >= k {
            let outgoing = window_codes[slot];
            base_counts[outgoing] -= 1;
            if base_counts[outgoing] == 0 {
                distinct_bases -= 1;
            }
        }
        window_codes[slot] = code;
        if base_counts[code] == 0 {
            distinct_bases += 1;
        }
        base_counts[code] += 1;
        valid_length += 1;
        forward = ((forward << 2) | code_u64) & mask;
        reverse = (reverse >> 2) | ((3 - code_u64) << reverse_shift);
        if valid_length < k
            || distinct_bases <= 2
            || base_counts.iter().copied().max().unwrap_or(0) >= complexity_threshold
        {
            continue;
        }

        let canonical = forward.min(reverse);
        let position = index + 1 - k;
        if let Some(observed) = repeated_positions.get_mut(&canonical) {
            if observed.len() >= position_cap {
                overflowed.insert(canonical);
            } else {
                observed.push(position);
            }
            continue;
        }
        if let Some(first_position) = first_positions.remove(&canonical) {
            repeated_positions.insert(canonical, vec![first_position, position]);
        } else {
            first_positions.insert(canonical, position);
        }
    }

    if min_seed_occurrences > 2 {
        repeated_positions.retain(|_, positions| positions.len() >= min_seed_occurrences);
    }
    (repeated_positions, overflowed.len())
}

fn build_spacing_histogram(
    repeated_positions: &HashMap<u64, Vec<usize>>,
    min_period: usize,
    max_period: usize,
    max_pairs_per_kmer: usize,
) -> Vec<usize> {
    let mut histogram = vec![0_usize; max_period + 1];
    let inspection_cap = max_pairs_per_kmer.saturating_mul(20);
    for positions in repeated_positions.values() {
        let mut pair_count = 0_usize;
        let mut inspected_pairs = 0_usize;
        'gaps: for position_gap in 1..positions.len() {
            for left_index in 0..positions.len() - position_gap {
                let spacing = positions[left_index + position_gap] - positions[left_index];
                inspected_pairs += 1;
                if (min_period..=max_period).contains(&spacing) {
                    let raw_bin = ((spacing + BIN_SIZE / 2) / BIN_SIZE) * BIN_SIZE;
                    let binned = raw_bin.clamp(min_period, max_period);
                    histogram[binned] += 1;
                    pair_count += 1;
                }
                if pair_count >= max_pairs_per_kmer || inspected_pairs >= inspection_cap {
                    break 'gaps;
                }
            }
        }
    }
    histogram
}

fn select_candidate_periods(
    histogram: &[usize],
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_spacing_support: usize,
) -> Vec<(usize, usize)> {
    let mut peaks: Vec<(usize, usize)> = (min_period..=max_period)
        .filter_map(|period| {
            let support = histogram[period];
            (support >= min_spacing_support).then_some((period, support))
        })
        .collect();
    peaks.sort_unstable_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
    peaks.truncate(top_periods);
    peaks
}

fn circular_distance(left: usize, right: usize, period: usize) -> usize {
    let direct = left.abs_diff(right);
    direct.min(period - direct)
}

fn modulo_periodicity_score(positions: &[&Vec<usize>], period: usize, tolerance: usize) -> f64 {
    let mut supported = 0_usize;
    let mut total = 0_usize;
    for observed in positions {
        if observed.len() < 2 {
            continue;
        }
        let residues: Vec<usize> = observed.iter().map(|position| position % period).collect();
        let group_support = residues
            .iter()
            .map(|center| {
                residues
                    .iter()
                    .filter(|residue| circular_distance(**residue, *center, period) <= tolerance)
                    .count()
            })
            .max()
            .unwrap_or(0);
        supported += group_support;
        total += residues.len();
    }
    if total == 0 {
        0.0
    } else {
        supported as f64 / total as f64
    }
}

fn bounded_periodicity_score(sequence: &[u8], period: usize) -> f64 {
    let compared_span = sequence.len().saturating_sub(period);
    if compared_span == 0 {
        return 0.0;
    }
    let step = (compared_span / MAX_IDENTITY_COMPARISONS).max(1);
    let mut matches = 0_usize;
    let mut valid = 0_usize;
    for index in (0..compared_span).step_by(step) {
        let left = sequence[index].to_ascii_uppercase();
        let right = sequence[index + period].to_ascii_uppercase();
        if left == b'N' || right == b'N' {
            continue;
        }
        valid += 1;
        if left == right {
            matches += 1;
        }
        if valid >= MAX_IDENTITY_COMPARISONS {
            break;
        }
    }
    if valid == 0 {
        0.0
    } else {
        matches as f64 / valid as f64
    }
}

fn refine_candidate_period(
    sequence: &[u8],
    repeated_positions: &HashMap<u64, Vec<usize>>,
    candidate_periods: &[usize],
    min_period: usize,
    max_period: usize,
) -> (usize, f64) {
    let mut local_periods = BTreeSet::new();
    for &candidate in candidate_periods {
        for offset in -REFINEMENT_RADIUS..=REFINEMENT_RADIUS {
            let period = candidate as isize + offset;
            if period >= min_period as isize && period <= max_period as isize {
                local_periods.insert(period as usize);
            }
        }
    }
    let mut seed_groups: Vec<(u64, &Vec<usize>)> = repeated_positions
        .iter()
        .map(|(&seed, positions)| (seed, positions))
        .collect();
    seed_groups.sort_unstable_by(|left, right| {
        right.1.len().cmp(&left.1.len()).then(left.0.cmp(&right.0))
    });
    let refinement_positions: Vec<&Vec<usize>> = seed_groups
        .into_iter()
        .take(MAX_SEED_GROUPS)
        .map(|(_, positions)| positions)
        .collect();
    let minimum_viable_identity = (ACCEPTANCE_SCORE - SEED_WEIGHT) / IDENTITY_WEIGHT;
    let mut best_period = 0_usize;
    let mut best_score = 0.0_f64;
    for period in local_periods {
        let identity_score = bounded_periodicity_score(sequence, period);
        let score = if identity_score < minimum_viable_identity {
            IDENTITY_WEIGHT * identity_score
        } else {
            let seed_score = modulo_periodicity_score(&refinement_positions, period, 2);
            SEED_WEIGHT * seed_score + IDENTITY_WEIGHT * identity_score
        };
        if score > best_score || (score == best_score && (best_period == 0 || period < best_period))
        {
            best_period = period;
            best_score = score;
        }
    }
    (best_period, best_score)
}

fn validate_scan_parameters(
    k: usize,
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_seed_occurrences: usize,
    min_spacing_support: usize,
    max_pairs_per_kmer: usize,
) -> PyResult<()> {
    if !(1..=31).contains(&k) {
        return Err(PyValueError::new_err("Rust backend requires k in 1..=31"));
    }
    if min_period == 0 || max_period < min_period {
        return Err(PyValueError::new_err("invalid period range"));
    }
    if top_periods == 0
        || min_seed_occurrences < 2
        || min_spacing_support == 0
        || max_pairs_per_kmer == 0
    {
        return Err(PyValueError::new_err("invalid seed/period limits"));
    }
    Ok(())
}

fn scan_owned_sequence(
    owned_sequence: Vec<u8>,
    k: usize,
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_seed_occurrences: usize,
    min_spacing_support: usize,
    max_pairs_per_kmer: usize,
) -> ScanResult {
    let (repeated_positions, overflow_count) = extract_repeated_positions(
        &owned_sequence,
        k,
        min_seed_occurrences,
        max_pairs_per_kmer,
    );
    let histogram = build_spacing_histogram(
        &repeated_positions,
        min_period,
        max_period,
        max_pairs_per_kmer,
    );
    let peaks = select_candidate_periods(
        &histogram,
        min_period,
        max_period,
        top_periods,
        min_spacing_support,
    );
    let candidate_periods: Vec<usize> = peaks.iter().map(|(period, _)| *period).collect();
    if candidate_periods.is_empty() {
        return ScanResult {
            candidate_periods,
            spacing_support: peaks,
            best_period: 0,
            periodicity_score: 0.0,
            overflow_count,
            status: "no_spacing_peak".to_string(),
        };
    }
    let (best_period, periodicity_score) = refine_candidate_period(
        &owned_sequence,
        &repeated_positions,
        &candidate_periods,
        min_period,
        max_period,
    );
    ScanResult {
        candidate_periods,
        spacing_support: peaks,
        best_period,
        periodicity_score,
        overflow_count,
        status: if periodicity_score >= ACCEPTANCE_SCORE {
            "accepted".to_string()
        } else {
            "rejected_score".to_string()
        },
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (sequence, k, min_period, max_period, top_periods, min_seed_occurrences, min_spacing_support, max_pairs_per_kmer))]
fn scan_read_for_periods(
    py: Python<'_>,
    sequence: &str,
    k: usize,
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_seed_occurrences: usize,
    min_spacing_support: usize,
    max_pairs_per_kmer: usize,
) -> PyResult<ScanResult> {
    validate_scan_parameters(
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
    )?;
    let owned_sequence = sequence.as_bytes().to_vec();
    Ok(py.allow_threads(move || {
        scan_owned_sequence(
            owned_sequence,
            k,
            min_period,
            max_period,
            top_periods,
            min_seed_occurrences,
            min_spacing_support,
            max_pairs_per_kmer,
        )
    }))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (sequences, k, min_period, max_period, top_periods, min_seed_occurrences, min_spacing_support, max_pairs_per_kmer))]
fn scan_reads_for_periods(
    py: Python<'_>,
    sequences: Vec<String>,
    k: usize,
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_seed_occurrences: usize,
    min_spacing_support: usize,
    max_pairs_per_kmer: usize,
) -> PyResult<Vec<ScanResult>> {
    validate_scan_parameters(
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
    )?;
    let owned_sequences: Vec<Vec<u8>> = sequences
        .into_iter()
        .map(|sequence| sequence.into_bytes())
        .collect();
    Ok(py.allow_threads(move || {
        owned_sequences
            .into_iter()
            .map(|owned_sequence| {
                scan_owned_sequence(
                    owned_sequence,
                    k,
                    min_period,
                    max_period,
                    top_periods,
                    min_seed_occurrences,
                    min_spacing_support,
                    max_pairs_per_kmer,
                )
            })
            .collect()
    }))
}

#[pymodule]
fn _rust_core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ScanResult>()?;
    module.add_class::<DiagnosticKmerCounter>()?;
    module.add_function(wrap_pyfunction!(scan_read_for_periods, module)?)?;
    module.add_function(wrap_pyfunction!(scan_reads_for_periods, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_non_default_repeat_period() {
        let monomer = b"ACGTTCAGGACTAACCGTGATCGATCGATCG";
        let sequence = monomer.repeat(20);
        let (positions, _) = extract_repeated_positions(&sequence, 11, 2, 100);
        let histogram = build_spacing_histogram(&positions, 20, 100, 100);
        let peaks = select_candidate_periods(&histogram, 20, 100, 3, 2);
        assert!(peaks[0].0.abs_diff(monomer.len()) <= 2);
        let candidates: Vec<usize> = peaks.iter().map(|(period, _)| *period).collect();
        let (best_period, score) =
            refine_candidate_period(&sequence, &positions, &candidates, 20, 100);
        assert_eq!(best_period, monomer.len());
        assert!(score >= ACCEPTANCE_SCORE);
    }

    #[test]
    fn pair_cap_bounds_histogram_support() {
        let positions = HashMap::from([(1_u64, (0..1000).step_by(10).collect())]);
        let histogram = build_spacing_histogram(&positions, 10, 100, 3);
        assert!(histogram.iter().sum::<usize>() <= 3);
    }

    #[test]
    fn ambiguous_base_resets_rolling_seed() {
        let seed = b"ACGTTCAGGAC";
        let mut sequence = seed.to_vec();
        sequence.push(b'N');
        sequence.extend_from_slice(seed);
        let (positions, _) = extract_repeated_positions(&sequence, 11, 2, 100);
        assert_eq!(positions.values().next(), Some(&vec![0, 12]));
    }

    #[test]
    fn diagnostic_counter_counts_only_targets() {
        let target = "ACGTTCAGGAC".to_string();
        let mut counter = DiagnosticKmerCounter::new(11, vec![target.clone()]).unwrap();
        counter.count_sequence("ACGTTCAGGACNACGTTCAGGAC");
        assert_eq!(counter.counts().get(&target), Some(&2));
    }
}
