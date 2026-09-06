//! Sparse read-level weighted target counts; no reads or moments persist here.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;

use super::{base_code, canonical_code};

type SparseRead = Vec<(usize, f64)>;

#[pyclass(frozen)]
pub struct WeightedKmerCounter {
    k: usize,
    targets: HashMap<u64, (usize, f64)>,
}

#[pymethods]
impl WeightedKmerCounter {
    #[new]
    fn new(k: usize, family_count: usize, targets: Vec<(String, usize, f64)>) -> PyResult<Self> {
        if !(1..=31).contains(&k) || family_count == 0 {
            return Err(PyValueError::new_err(
                "Require k in 1..31 and positive family count",
            ));
        }
        let mut bank = HashMap::new();
        for (word, family, weight) in targets {
            if word.len() != k
                || family >= family_count
                || !weight.is_finite()
                || !(0.0 < weight && weight <= 1.0)
            {
                return Err(PyValueError::new_err("Invalid weighted diagnostic target"));
            }
            let code = canonical_code(word.as_bytes())
                .ok_or_else(|| PyValueError::new_err("Invalid target DNA"))?;
            if bank.insert(code, (family, weight)).is_some() {
                return Err(PyValueError::new_err(
                    "Duplicate canonical target; require exclusive words",
                ));
            }
        }
        Ok(Self { k, targets: bank })
    }

    fn count_sequences(&self, py: Python<'_>, sequences: Vec<String>) -> PyResult<Vec<SparseRead>> {
        if sequences.iter().any(|s| {
            s.is_empty()
                || s.bytes()
                    .any(|b| base_code(b).is_none() && !matches!(b, b'N' | b'n'))
        }) {
            return Err(PyValueError::new_err(
                "Reads require nonempty ACGTN sequences",
            ));
        }
        Ok(py.allow_threads(|| {
            sequences
                .iter()
                .map(|s| self.count_record(s.as_bytes()))
                .collect()
        }))
    }
}

impl WeightedKmerCounter {
    fn count_record(&self, sequence: &[u8]) -> SparseRead {
        let mask = (1_u64 << (2 * self.k)) - 1;
        let shift = 2 * (self.k - 1);
        let (mut forward, mut reverse, mut valid) = (0_u64, 0_u64, 0_usize);
        let mut values: HashMap<usize, f64> = HashMap::new();
        for &base in sequence {
            let Some(code) = base_code(base) else {
                forward = 0;
                reverse = 0;
                valid = 0;
                continue;
            };
            valid += 1;
            forward = ((forward << 2) | code) & mask;
            reverse = (reverse >> 2) | ((3 - code) << shift);
            if valid >= self.k {
                if let Some(&(family, weight)) = self.targets.get(&forward.min(reverse)) {
                    *values.entry(family).or_default() += weight;
                }
            }
        }
        let mut result: SparseRead = values.into_iter().collect();
        result.sort_unstable_by_key(|&(family, _)| family);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_counts_reset_at_n_and_merge_reverse_complements() {
        let counter =
            WeightedKmerCounter::new(3, 2, vec![("ACG".into(), 0, 0.5), ("AAA".into(), 1, 0.1)])
                .unwrap();
        assert_eq!(
            counter.count_record(b"ACGNNNCGTAAATTT"),
            vec![(0, 1.0), (1, 0.2)]
        );
        assert_eq!(counter.count_record(b"acgnnncgt"), vec![(0, 1.0)]);
        assert!(counter.count_record(b"NNNN").is_empty());
        assert!(counter.count_record(b"AC").is_empty());
    }

    #[test]
    fn malformed_and_ambiguous_target_banks_fail() {
        assert!(WeightedKmerCounter::new(
            3,
            2,
            vec![("ACG".into(), 0, 0.5), ("CGT".into(), 1, 0.1)]
        )
        .is_err());
        assert!(WeightedKmerCounter::new(3, 2, vec![("NNN".into(), 0, 0.5)]).is_err());
        assert!(WeightedKmerCounter::new(3, 2, vec![("ACG".into(), 2, 0.5)]).is_err());
        assert!(WeightedKmerCounter::new(3, 2, vec![("ACG".into(), 0, f64::NAN)]).is_err());
    }
}
