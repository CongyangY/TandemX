use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

const BIN_SIZE: usize = 5;
const REFINEMENT_RADIUS: isize = 2;
const MAX_SEED_GROUPS: usize = 128;
const ACCEPTANCE_SCORE: f64 = 0.75;
const SEED_WEIGHT: f64 = 0.20;
const IDENTITY_WEIGHT: f64 = 0.80;
const SHORT_PERIOD_SCAN_MAX: usize = 19;
const SHORT_PERIOD_ACCEPTANCE_SCORE: f64 = 0.80;

#[pyclass(frozen)]
struct SequenceStatsResult {
    #[pyo3(get)]
    record_count: usize,
    #[pyo3(get)]
    total_bases: usize,
    #[pyo3(get)]
    max_read_length: usize,
}

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
    repeat_start: usize,
    #[pyo3(get)]
    repeat_end: usize,
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

fn is_valid_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'N' | b'n'
    )
}

fn normalize_sequence_line<'a>(
    line: &'a str,
    path: &str,
    line_number: usize,
) -> PyResult<&'a [u8]> {
    let bytes = line.trim().as_bytes();
    if bytes.is_empty() {
        return Err(PyValueError::new_err(format!(
            "Empty sequence line in {} at line {}",
            path, line_number
        )));
    }
    if let Some(invalid) = bytes.iter().copied().find(|base| !is_valid_base(*base)) {
        return Err(PyValueError::new_err(format!(
            "Invalid base '{}' in {} at line {}",
            invalid as char, path, line_number
        )));
    }
    Ok(bytes)
}

fn detect_sequence_format(path: &str) -> PyResult<&'static str> {
    let lowercase = path.to_lowercase();
    let name = lowercase.strip_suffix(".gz").unwrap_or(&lowercase);
    if name.ends_with(".fa") || name.ends_with(".fasta") {
        return Ok("fasta");
    }
    if name.ends_with(".fq") || name.ends_with(".fastq") {
        return Ok("fastq");
    }
    Err(PyValueError::new_err(format!(
        "Unsupported sequence file extension for {}. Expected .fa, .fasta, .fq, .fastq, or .gz-compressed variants.",
        path
    )))
}

fn open_sequence_reader(path: &str) -> PyResult<Box<dyn BufRead>> {
    if path.to_lowercase().ends_with(".gz") {
        return Err(PyValueError::new_err(format!(
            "Rust sequence stats do not support gzip-compressed inputs: {}",
            path
        )));
    }
    let file = File::open(path)
        .map_err(|err| PyValueError::new_err(format!("Failed to open {}: {}", path, err)))?;
    Ok(Box::new(BufReader::new(file)))
}

fn count_fasta_stats(reader: &mut dyn BufRead, path: &str) -> PyResult<SequenceStatsResult> {
    let mut line = String::new();
    let mut line_number = 0_usize;
    let mut record_count = 0_usize;
    let mut total_bases = 0_usize;
    let mut max_read_length = 0_usize;
    let mut current_len = 0_usize;
    let mut current_header = false;
    let mut yielded = false;
    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|err| PyValueError::new_err(format!("Failed reading {}: {}", path, err)))?;
        if read == 0 {
            break;
        }
        line_number += 1;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('>') {
            if current_header {
                record_count += 1;
                total_bases += current_len;
                max_read_length = max_read_length.max(current_len);
                yielded = true;
            }
            if trimmed.len() == 1 {
                return Err(PyValueError::new_err(format!(
                    "Empty FASTA header in {} at line {}",
                    path, line_number
                )));
            }
            current_header = true;
            current_len = 0;
            continue;
        }
        if !current_header {
            return Err(PyValueError::new_err(format!(
                "Invalid FASTA in {}: sequence before header at line {}",
                path, line_number
            )));
        }
        let sequence = normalize_sequence_line(trimmed, path, line_number)?;
        current_len += sequence.len();
    }
    if current_header {
        record_count += 1;
        total_bases += current_len;
        max_read_length = max_read_length.max(current_len);
        yielded = true;
    }
    if !yielded {
        return Err(PyValueError::new_err(format!(
            "Sequence file is empty or contains no records: {}",
            path
        )));
    }
    Ok(SequenceStatsResult {
        record_count,
        total_bases,
        max_read_length,
    })
}

fn count_fastq_stats(reader: &mut dyn BufRead, path: &str) -> PyResult<SequenceStatsResult> {
    let mut line = String::new();
    let mut line_number = 0_usize;
    let mut record_count = 0_usize;
    let mut total_bases = 0_usize;
    let mut max_read_length = 0_usize;
    let mut yielded = false;
    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|err| PyValueError::new_err(format!("Failed reading {}: {}", path, err)))?;
        if read == 0 {
            break;
        }
        line_number += 1;
        let header = line.trim_end_matches(&['\n', '\r'][..]).trim();
        if header.is_empty() {
            continue;
        }
        if !header.starts_with('@') || header.len() == 1 {
            return Err(PyValueError::new_err(format!(
                "Invalid FASTQ header in {} at line {}",
                path, line_number
            )));
        }

        let mut sequence = String::new();
        let mut plus = String::new();
        let mut quality = String::new();
        if reader
            .read_line(&mut sequence)
            .map_err(|err| PyValueError::new_err(format!("Failed reading {}: {}", path, err)))?
            == 0
            || reader
                .read_line(&mut plus)
                .map_err(|err| PyValueError::new_err(format!("Failed reading {}: {}", path, err)))?
                == 0
            || reader
                .read_line(&mut quality)
                .map_err(|err| PyValueError::new_err(format!("Failed reading {}: {}", path, err)))?
                == 0
        {
            return Err(PyValueError::new_err(format!(
                "Truncated FASTQ record in {} starting at line {}",
                path, line_number
            )));
        }
        let sequence_line = line_number + 1;
        line_number += 3;
        let sequence_trimmed = sequence.trim_end_matches(&['\n', '\r'][..]).trim();
        let plus_trimmed = plus.trim_end_matches(&['\n', '\r'][..]).trim();
        let quality_trimmed = quality.trim_end_matches(&['\n', '\r'][..]).trim();
        if !plus_trimmed.starts_with('+') {
            return Err(PyValueError::new_err(format!(
                "Invalid FASTQ separator in {} at line {}",
                path,
                sequence_line + 1
            )));
        }
        let bases = normalize_sequence_line(sequence_trimmed, path, sequence_line)?;
        if bases.len() != quality_trimmed.len() {
            return Err(PyValueError::new_err(format!(
                "FASTQ sequence and quality lengths differ in {} at line {}",
                path, sequence_line
            )));
        }
        record_count += 1;
        total_bases += bases.len();
        max_read_length = max_read_length.max(bases.len());
        yielded = true;
    }
    if !yielded {
        return Err(PyValueError::new_err(format!(
            "Sequence file is empty or contains no records: {}",
            path
        )));
    }
    Ok(SequenceStatsResult {
        record_count,
        total_bases,
        max_read_length,
    })
}

fn count_sequence_file_stats_inner(path: &str) -> PyResult<SequenceStatsResult> {
    let mut reader = open_sequence_reader(path)?;
    match detect_sequence_format(path)? {
        "fasta" => count_fasta_stats(reader.as_mut(), path),
        "fastq" => count_fastq_stats(reader.as_mut(), path),
        _ => Err(PyValueError::new_err(format!(
            "Unsupported format for {}",
            path
        ))),
    }
}

#[pyfunction]
fn count_sequence_file_stats(py: Python<'_>, path: &str) -> PyResult<SequenceStatsResult> {
    let owned = path.to_string();
    py.allow_threads(move || count_sequence_file_stats_inner(&owned))
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

    fn count_sequence(&mut self, py: Python<'_>, sequence: String) {
        py.allow_threads(|| self.count_sequence_bytes(sequence.as_bytes()));
    }

    fn count_sequences(&mut self, py: Python<'_>, sequences: Vec<String>) {
        py.allow_threads(|| {
            for sequence in sequences {
                self.count_sequence_bytes(sequence.as_bytes());
            }
        });
    }

    fn counts(&self) -> HashMap<String, u64> {
        self.target_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), self.counts[index]))
            .collect()
    }
}

impl DiagnosticKmerCounter {
    fn count_sequence_bytes(&mut self, sequence: &[u8]) {
        let mask = (1_u64 << (2 * self.k)) - 1;
        let reverse_shift = 2 * (self.k - 1);
        let mut forward = 0_u64;
        let mut reverse = 0_u64;
        let mut valid_length = 0_usize;
        for &base in sequence {
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
            || base_counts.iter().copied().max().unwrap_or(0) >= complexity_threshold
            || (distinct_bases <= 2 && is_simple_periodic_window(&window_codes, index, k))
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

fn is_simple_periodic_window(window_codes: &[usize], index: usize, k: usize) -> bool {
    let ordered: Vec<usize> = (0..k)
        .map(|offset| window_codes[(index + 1 + offset) % k])
        .collect();
    if ordered.iter().all(|code| *code == ordered[0]) {
        return true;
    }
    k >= 4 && (2..k).all(|offset| ordered[offset] == ordered[offset % 2])
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

fn best_local_periodicity_score(
    sequence: &[u8],
    period: usize,
    min_repeat_span: usize,
    acceptance_score: f64,
) -> (f64, usize, usize) {
    let compared_span = sequence.len().saturating_sub(period);
    let min_compared = min_repeat_span.saturating_sub(period).max(period).max(1);
    if period == 0 || compared_span < min_compared {
        return (0.0, 0, 0);
    }
    let mismatch_penalty = acceptance_score / (1.0 - acceptance_score).max(1e-12);
    let mut prefix_score = Vec::with_capacity(compared_span + 1);
    let mut prefix_matches = Vec::with_capacity(compared_span + 1);
    let mut prefix_valid = Vec::with_capacity(compared_span + 1);
    prefix_score.push(0.0);
    prefix_matches.push(0_usize);
    prefix_valid.push(0_usize);
    for index in 0..compared_span {
        let left = sequence[index];
        let right = sequence[index + period];
        let ambiguous = matches!(left, b'N' | b'n') || matches!(right, b'N' | b'n');
        let is_match = !ambiguous && left.eq_ignore_ascii_case(&right);
        let value = if ambiguous {
            0.0
        } else if is_match {
            1.0
        } else {
            -mismatch_penalty
        };
        prefix_score.push(prefix_score.last().copied().unwrap_or(0.0) + value);
        prefix_matches.push(prefix_matches.last().copied().unwrap_or(0) + usize::from(is_match));
        prefix_valid.push(prefix_valid.last().copied().unwrap_or(0) + usize::from(!ambiguous));
    }

    let mut minimum_prefix_value = prefix_score[0];
    let mut minimum_prefix_index = 0_usize;
    let mut best_score = f64::NEG_INFINITY;
    let mut best_start = 0_usize;
    let mut best_end = 0_usize;
    for end in min_compared..=compared_span {
        let eligible = end - min_compared;
        let eligible_value = prefix_score[eligible];
        if eligible_value < minimum_prefix_value {
            minimum_prefix_value = eligible_value;
            minimum_prefix_index = eligible;
        }
        let interval_score = prefix_score[end] - minimum_prefix_value;
        let interval_length = end - minimum_prefix_index;
        let best_length = best_end - best_start;
        if interval_score > best_score
            || (interval_score == best_score && interval_length > best_length)
            || (interval_score == best_score
                && interval_length == best_length
                && minimum_prefix_index < best_start)
        {
            best_score = interval_score;
            best_start = minimum_prefix_index;
            best_end = end;
        }
    }
    let valid = prefix_valid[best_end] - prefix_valid[best_start];
    let matches = prefix_matches[best_end] - prefix_matches[best_start];
    let identity = if valid == 0 {
        0.0
    } else {
        matches as f64 / valid as f64
    };
    (
        identity,
        best_start,
        (best_end + period).min(sequence.len()),
    )
}

fn direct_period_scan(
    sequence: &[u8],
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_repeat_span: usize,
) -> ScanResult {
    let mut scored_periods: Vec<(usize, f64, usize, usize)> = (min_period..=max_period)
        .map(|period| {
            let (score, start, end) = best_local_periodicity_score(
                sequence,
                period,
                min_repeat_span,
                SHORT_PERIOD_ACCEPTANCE_SCORE,
            );
            (period, score, start, end)
        })
        .collect();
    scored_periods.sort_unstable_by(|left, right| {
        right
            .1
            .partial_cmp(&left.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then((right.3 - right.2).cmp(&(left.3 - left.2)))
            .then(left.0.cmp(&right.0))
    });
    let best = scored_periods.first().copied().unwrap_or((0, 0.0, 0, 0));
    let candidate_periods = scored_periods
        .iter()
        .take(top_periods)
        .filter(|(_, score, _, _)| *score >= SHORT_PERIOD_ACCEPTANCE_SCORE)
        .map(|(period, _, _, _)| *period)
        .collect();
    ScanResult {
        candidate_periods,
        spacing_support: Vec::new(),
        best_period: best.0,
        periodicity_score: best.1,
        repeat_start: best.2,
        repeat_end: best.3,
        overflow_count: 0,
        status: if best.1 >= SHORT_PERIOD_ACCEPTANCE_SCORE {
            "accepted".to_string()
        } else {
            "rejected_score".to_string()
        },
    }
}

fn refine_candidate_period(
    sequence: &[u8],
    repeated_positions: &HashMap<u64, Vec<usize>>,
    candidate_periods: &[usize],
    min_period: usize,
    max_period: usize,
    min_repeat_span: usize,
) -> (usize, f64, usize, usize) {
    let mut local_periods = BTreeSet::new();
    for &candidate in candidate_periods {
        for offset in -REFINEMENT_RADIUS..=REFINEMENT_RADIUS {
            let period = candidate as isize + offset;
            if period >= min_period as isize && period <= max_period as isize {
                local_periods.insert(period as usize);
                let local = period as usize;
                let mut divisor = 2_usize;
                while divisor.saturating_mul(divisor) <= local {
                    if local.is_multiple_of(divisor) {
                        for fundamental in [divisor, local / divisor] {
                            if (min_period..=max_period).contains(&fundamental) {
                                local_periods.insert(fundamental);
                            }
                        }
                    }
                    divisor += 1;
                }
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
    let mut best_start = 0_usize;
    let mut best_end = 0_usize;
    for period in local_periods {
        let (identity_score, start, end) = best_local_periodicity_score(
            sequence,
            period,
            min_repeat_span,
            minimum_viable_identity,
        );
        let score = if identity_score < minimum_viable_identity {
            IDENTITY_WEIGHT * identity_score
        } else {
            let seed_score = modulo_periodicity_score(&refinement_positions, period, 2);
            SEED_WEIGHT * seed_score + IDENTITY_WEIGHT * identity_score
        };
        if score > best_score
            || (score == best_score && end - start > best_end - best_start)
            || (score == best_score
                && end - start == best_end - best_start
                && (best_period == 0 || period < best_period))
        {
            best_period = period;
            best_score = score;
            best_start = start;
            best_end = end;
        }
    }
    (best_period, best_score, best_start, best_end)
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

#[derive(Clone, Copy)]
struct ScanParameters {
    k: usize,
    min_period: usize,
    max_period: usize,
    top_periods: usize,
    min_seed_occurrences: usize,
    min_spacing_support: usize,
    max_pairs_per_kmer: usize,
    min_repeat_span: usize,
}

fn scan_owned_sequence(owned_sequence: Vec<u8>, parameters: ScanParameters) -> ScanResult {
    let short_max = parameters.max_period.min(SHORT_PERIOD_SCAN_MAX);
    if parameters.min_period <= short_max {
        let short_result = direct_period_scan(
            &owned_sequence,
            parameters.min_period,
            short_max,
            parameters.top_periods,
            parameters.min_repeat_span,
        );
        if short_result.status == "accepted" {
            return short_result;
        }
    }
    let long_min_period = parameters.min_period.max(SHORT_PERIOD_SCAN_MAX + 1);
    if long_min_period > parameters.max_period {
        return ScanResult {
            candidate_periods: Vec::new(),
            spacing_support: Vec::new(),
            best_period: 0,
            periodicity_score: 0.0,
            repeat_start: 0,
            repeat_end: 0,
            overflow_count: 0,
            status: "no_spacing_peak".to_string(),
        };
    }
    let (repeated_positions, overflow_count) = extract_repeated_positions(
        &owned_sequence,
        parameters.k,
        parameters.min_seed_occurrences,
        parameters.max_pairs_per_kmer,
    );
    let histogram = build_spacing_histogram(
        &repeated_positions,
        long_min_period,
        parameters.max_period,
        parameters.max_pairs_per_kmer,
    );
    let peaks = select_candidate_periods(
        &histogram,
        long_min_period,
        parameters.max_period,
        parameters.top_periods,
        parameters.min_spacing_support,
    );
    let candidate_periods: Vec<usize> = peaks.iter().map(|(period, _)| *period).collect();
    if candidate_periods.is_empty() {
        return ScanResult {
            candidate_periods,
            spacing_support: peaks,
            best_period: 0,
            periodicity_score: 0.0,
            repeat_start: 0,
            repeat_end: 0,
            overflow_count,
            status: "no_spacing_peak".to_string(),
        };
    }
    let (best_period, periodicity_score, repeat_start, repeat_end) = refine_candidate_period(
        &owned_sequence,
        &repeated_positions,
        &candidate_periods,
        long_min_period,
        parameters.max_period,
        parameters.min_repeat_span,
    );
    ScanResult {
        candidate_periods,
        spacing_support: peaks,
        best_period,
        periodicity_score,
        repeat_start,
        repeat_end,
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
#[pyo3(signature = (sequence, k, min_period, max_period, top_periods, min_seed_occurrences, min_spacing_support, max_pairs_per_kmer, min_repeat_span=1))]
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
    min_repeat_span: usize,
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
    let parameters = ScanParameters {
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
        min_repeat_span,
    };
    Ok(py.allow_threads(move || scan_owned_sequence(owned_sequence, parameters)))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (sequences, k, min_period, max_period, top_periods, min_seed_occurrences, min_spacing_support, max_pairs_per_kmer, min_repeat_span=1))]
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
    min_repeat_span: usize,
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
    let parameters = ScanParameters {
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
        min_repeat_span,
    };
    Ok(py.allow_threads(move || {
        owned_sequences
            .into_iter()
            .map(|owned_sequence| scan_owned_sequence(owned_sequence, parameters))
            .collect()
    }))
}

#[pymodule]
fn _rust_core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<SequenceStatsResult>()?;
    module.add_class::<ScanResult>()?;
    module.add_class::<DiagnosticKmerCounter>()?;
    module.add_function(wrap_pyfunction!(count_sequence_file_stats, module)?)?;
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
        let (best_period, score, start, end) =
            refine_candidate_period(&sequence, &positions, &candidates, 20, 100, 100);
        assert_eq!(best_period, monomer.len());
        assert!(score >= ACCEPTANCE_SCORE);
        assert_eq!(start, 0);
        assert_eq!(end, sequence.len());
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
        counter.count_sequence_bytes(b"ACGTTCAGGACNACGTTCAGGAC");
        assert_eq!(counter.counts().get(&target), Some(&2));
    }
}
