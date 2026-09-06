use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn distance(a: &[u8], b: &[u8], limit: usize) -> usize {
    if a.len().abs_diff(b.len()) > limit {
        return limit + 1;
    }
    let mut start = 0;
    while start < a.len().min(b.len()) && a[start] == b[start] && b"ACGT".contains(&a[start]) {
        start += 1;
    }
    let (mut a, mut b) = (&a[start..], &b[start..]);
    while !a.is_empty()
        && !b.is_empty()
        && a.last() == b.last()
        && b"ACGT".contains(a.last().unwrap())
    {
        a = &a[..a.len() - 1];
        b = &b[..b.len() - 1];
    }
    if a.is_empty() || b.is_empty() {
        return (limit + 1).min(a.len().max(b.len()));
    }
    let width = 2 * limit + 1;
    let mut previous = vec![limit + 1; width];
    for j in 0..=b.len().min(limit) {
        previous[j + limit] = j;
    }
    for i in 1..=a.len() {
        let mut current = vec![limit + 1; width];
        for j in i.saturating_sub(limit)..=b.len().min(i + limit) {
            let slot = j + limit - i;
            if j == 0 {
                current[slot] = i;
                continue;
            }
            let mismatch = usize::from(a[i - 1] != b[j - 1] || !b"ACGT".contains(&a[i - 1]));
            let mut score = previous[slot] + mismatch;
            if slot + 1 < width {
                score = score.min(previous[slot + 1] + 1);
            }
            if slot > 0 {
                score = score.min(current[slot - 1] + 1);
            }
            current[slot] = score;
        }
        if *current.iter().min().unwrap() > limit {
            return limit + 1;
        }
        previous = current;
    }
    previous[b.len() + limit - a.len()].min(limit + 1)
}

#[pyfunction]
pub fn bounded_edit_distance(
    py: Python<'_>,
    reference: &str,
    query: &str,
    limit: usize,
) -> PyResult<usize> {
    if !reference.is_ascii() || !query.is_ascii() {
        return Err(PyValueError::new_err("ASCII sequences required"));
    }
    let cells = limit
        .checked_mul(2)
        .and_then(|n| n.checked_add(1))
        .and_then(|width| reference.len().max(1).checked_mul(width));
    if cells.is_none_or(|n| n > 32_000_000) {
        return Err(PyValueError::new_err(
            "Cluster comparison exceeds 32 million cells",
        ));
    }
    let a = reference.as_bytes().to_ascii_uppercase();
    let b = query.as_bytes().to_ascii_uppercase();
    Ok(py.allow_threads(move || distance(&a, &b, limit)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn edits_and_threshold() {
        assert_eq!(distance(b"ACGTA", b"ACCGTA", 1), 1);
        assert_eq!(distance(b"ACGTA", b"ATGTA", 0), 1);
        assert_eq!(distance(b"NN", b"NN", 2), 2);
    }
}
