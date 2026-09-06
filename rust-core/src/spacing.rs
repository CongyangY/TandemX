//! Expose the existing bounded seed histogram without legacy period refinement.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type HistogramResult = (Vec<(usize, usize)>, usize);

#[pyfunction]
pub fn seed_spacing_histogram(
    py: Python<'_>,
    sequence: &str,
    k: usize,
    min_period: usize,
    max_period: usize,
    min_seed_occurrences: usize,
    max_pairs_per_kmer: usize,
) -> PyResult<HistogramResult> {
    if !(1..=31).contains(&k)
        || min_period == 0
        || max_period == 0
        || min_seed_occurrences < 2
        || max_pairs_per_kmer == 0
        || max_pairs_per_kmer.checked_add(1).is_none()
        || !sequence.is_ascii()
    {
        return Err(PyValueError::new_err(
            "Invalid native seed histogram parameters",
        ));
    }
    let owned = sequence.as_bytes().to_vec();
    Ok(py.allow_threads(move || {
        let (positions, overflow) =
            super::extract_repeated_positions(&owned, k, min_seed_occurrences, max_pairs_per_kmer);
        // No observed spacing exceeds read length; retain enough extra bins for
        // rounding, without allocating according to a potentially enormous CLI cap.
        let upper = max_period.min(owned.len().saturating_add(super::BIN_SIZE));
        if min_period > upper {
            return (Vec::new(), overflow);
        }
        let histogram =
            super::build_spacing_histogram(&positions, min_period, upper, max_pairs_per_kmer);
        let nonzero = histogram
            .into_iter()
            .enumerate()
            .filter(|(_, count)| *count > 0)
            .collect();
        (nonzero, overflow)
    }))
}
