//! Exact native version of the legacy ungapped representative audit.
//! This is not the gapped/circular clustering criterion or biological identity.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn compare(a: &[u8], b: &[u8]) -> (f64, usize, &'static str) {
    let minimum_overlap = 50.min(a.len()).min(b.len());
    let reverse: Vec<u8> = b
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect();
    let mut best = (0.0, 0, "forward");
    for (orientation, oriented) in [("forward", b), ("reverse", reverse.as_slice())] {
        for offset in -(oriented.len() as isize) + 1..a.len() as isize {
            let start_a = offset.max(0) as usize;
            let start_b = (-offset).max(0) as usize;
            let overlap = (a.len() - start_a).min(oriented.len() - start_b);
            if overlap < minimum_overlap {
                continue;
            }
            let matches = a[start_a..start_a + overlap]
                .iter()
                .zip(&oriented[start_b..start_b + overlap])
                .filter(|(left, right)| left == right)
                .count();
            let identity = matches as f64 / overlap as f64;
            if identity > best.0 || (identity == best.0 && overlap > best.1) {
                best = (identity, overlap, orientation);
            }
        }
    }
    best
}

#[pyfunction]
pub fn ungapped_family_identity(
    py: Python<'_>,
    a: &str,
    b: &str,
) -> PyResult<(f64, usize, String)> {
    if a.is_empty() || b.is_empty() || !a.is_ascii() || !b.is_ascii() {
        return Err(PyValueError::new_err(
            "Need nonempty ASCII representative sequences",
        ));
    }
    let (identity, overlap, orientation) = py.allow_threads(|| compare(a.as_bytes(), b.as_bytes()));
    Ok((identity, overlap, orientation.to_string()))
}

#[cfg(test)]
mod tests {
    use super::compare;

    #[test]
    fn orientation_and_ties_match_reference_contract() {
        assert_eq!(compare(b"ACGACG", b"CGTCGT"), (1.0, 6, "reverse"));
        assert_eq!(compare(b"NNNN", b"NNNN"), (1.0, 4, "forward"));
        assert_eq!(compare(b"AAAA", b"CCCC"), (0.0, 4, "forward"));
        assert_eq!(compare(&[b'A'; 60], &[b'A'; 80]), (1.0, 60, "forward"));
    }
}
