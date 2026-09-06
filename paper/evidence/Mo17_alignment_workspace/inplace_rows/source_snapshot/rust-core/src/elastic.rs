//! Bounded, indel-aware read-local alignment. See Python reference alignment.py.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type Hit = (
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    i32,
    Vec<(usize, usize)>,
);

fn align(sequence: &[u8], period: usize, min_span: usize, band: usize, x_drop: i32) -> Vec<Hit> {
    let n = sequence.len();
    let width = 2 * band + 1;
    let lowest = period - band;
    if lowest >= n {
        return Vec::new();
    }
    let mut directions = vec![0_u8; (n + 1) * width];
    let mut previous = vec![0_i32; width];
    let mut previous_peak = vec![0_i32; width];
    // Traverse bands left to right in place: b and b+1 still hold the previous
    // row, while b-1 already holds this row. Every rejected cell is reset, and
    // the inactive tail is reset only after its old value can feed the last up
    // edge. Thus no full-row allocation, copy, or second pair of buffers is needed.
    let mut endpoints = Vec::new();
    let minimum_score = (40.min(min_span) as i32)
        .max((0.7 * period.max(min_span.saturating_sub(period)) as f64).ceil() as i32);
    for i in 1..=n - lowest {
        let mut row_best = (0, 0);
        let left = sequence[i - 1];
        let active_width = if b"ACGT".contains(&left) {
            width.min(n - i - lowest + 1)
        } else {
            0
        };
        for b in 0..active_width {
            let j = i + lowest + b;
            if !b"ACGT".contains(&sequence[j - 1]) {
                previous[b] = 0;
                previous_peak[b] = 0;
                continue;
            }
            let mut score = previous[b] + if left == sequence[j - 1] { 2 } else { -3 };
            let mut direction = 1;
            let mut peak = previous_peak[b];
            if b + 1 < width && previous[b + 1] - 4 > score {
                score = previous[b + 1] - 4;
                direction = 2;
                peak = previous_peak[b + 1];
            }
            if b > 0 && previous[b - 1] - 4 > score {
                score = previous[b - 1] - 4;
                direction = 3;
                peak = previous_peak[b - 1];
            }
            if score <= 0 || peak - score > x_drop {
                previous[b] = 0;
                previous_peak[b] = 0;
                continue;
            }
            previous[b] = score;
            previous_peak[b] = peak.max(score);
            directions[i * width + b] = direction;
            if score > row_best.0 {
                row_best = (score, b);
            }
        }
        if row_best.0 >= minimum_score {
            endpoints.push((row_best.0, i, row_best.1));
        }
        previous[active_width..].fill(0);
        previous_peak[active_width..].fill(0);
    }
    endpoints.sort_unstable_by_key(|&(score, i, b)| (-score, i, b));
    let mut hits: Vec<Hit> = Vec::new();
    let mut claimed = vec![false; n + 1];
    for (score, endpoint, end_band) in endpoints {
        if claimed[endpoint] {
            continue;
        }
        let mut i = endpoint;
        let mut b = end_band as isize;
        let end = i + lowest + end_band;
        let mut pairs = Vec::new();
        let (mut matches, mut columns, mut gaps) = (0, 0, 0);
        while i > 0 && b >= 0 && b < width as isize {
            let direction = directions[i * width + b as usize];
            if direction == 0 {
                break;
            }
            let j = i + lowest + b as usize;
            columns += 1;
            match direction {
                1 => {
                    pairs.push((i - 1, j - 1));
                    matches += usize::from(sequence[i - 1] == sequence[j - 1]);
                    i -= 1;
                }
                2 => {
                    gaps += 1;
                    i -= 1;
                    b += 1;
                }
                _ => {
                    gaps += 1;
                    b -= 1;
                }
            }
        }
        if pairs.is_empty() || (matches as f64 / columns as f64) < 0.75 {
            continue;
        }
        pairs.reverse();
        let mut offsets: Vec<usize> = pairs.iter().map(|&(l, r)| r - l).collect();
        offsets.sort_unstable();
        let measured = offsets[(offsets.len() - 1) / 2];
        if end - i < min_span || ((endpoint - i) as f64) < 0.8 * measured as f64 {
            continue;
        }
        if hits
            .iter()
            .any(|h| end.min(h.1).saturating_sub(i.max(h.0)) * 2 >= (end - i).min(h.1 - h.0))
        {
            continue;
        }
        hits.push((i, end, measured, matches, columns, gaps, score, pairs));
        claimed[i..=end].fill(true);
    }
    hits.sort_unstable_by_key(|h| (h.0, h.1, h.2));
    hits
}

#[pyfunction]
pub fn banded_self_align(
    py: Python<'_>,
    sequence: &str,
    period: usize,
    min_span: usize,
    band: usize,
    x_drop: i32,
) -> PyResult<Vec<Hit>> {
    if period == 0 || min_span == 0 || x_drop <= 0 || band >= period {
        return Err(PyValueError::new_err(
            "invalid elastic alignment parameters",
        ));
    }
    let cells = sequence.len().checked_add(1).and_then(|n| {
        band.checked_mul(2)
            .and_then(|b| b.checked_add(1))
            .and_then(|width| n.checked_mul(width))
    });
    if cells.is_none_or(|cells| cells > 32_000_000) {
        return Err(PyValueError::new_err(
            "Self-alignment trace exceeds 32 million cells",
        ));
    }
    if !sequence.is_ascii() {
        return Err(PyValueError::new_err("sequence must be ASCII"));
    }
    let sequence = sequence.as_bytes().to_ascii_uppercase();
    Ok(py.allow_threads(move || align(&sequence, period, min_span, band, x_drop)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gap_free_alignment_has_complete_boundaries() {
        let unit = b"ACGTTCAGGACTAACCGTGATCGATCGATCG";
        let sequence = unit.repeat(10);
        let hits = align(&sequence, unit.len(), 100, 3, 40);
        assert_eq!(hits.len(), 1);
        assert_eq!(
            (hits[0].0, hits[0].1, hits[0].2),
            (0, sequence.len(), unit.len())
        );
        assert_eq!(hits[0].3, hits[0].4);
    }
    #[test]
    fn ambiguity_breaks_alignment() {
        assert!(align(&vec![b'N'; 500], 50, 100, 4, 40).is_empty());
    }
}

#[pyfunction]
pub fn global_align_ops(
    py: Python<'_>,
    reference: &str,
    query: &str,
    band: usize,
) -> PyResult<String> {
    let (n, m) = (reference.len(), query.len());
    let cells = n.checked_add(1).and_then(|n| {
        band.checked_mul(2)
            .and_then(|b| b.checked_add(1))
            .and_then(|w| n.checked_mul(w))
    });
    if cells.is_none_or(|cells| cells > 32_000_000) || n.abs_diff(m) > band {
        return Err(PyValueError::new_err(
            "Invalid or excessive unit alignment band",
        ));
    }
    if !reference.is_ascii() || !query.is_ascii() {
        return Err(PyValueError::new_err("sequence must be ASCII"));
    }
    let reference = reference.as_bytes().to_vec();
    let query = query.as_bytes().to_vec();
    Ok(py.allow_threads(move || global_ops(&reference, &query, band)))
}

fn global_ops(reference: &[u8], query: &[u8], band: usize) -> String {
    let (n, m) = (reference.len(), query.len());
    let width = 2 * band + 1;
    let mut directions = vec![0_u8; (n + 1) * width];
    let mut previous = vec![-1_000_000_000_i32; width];
    for j in 0..=m.min(band) {
        previous[j + band] = -4 * j as i32;
        directions[j + band] = if j > 0 { 3 } else { 0 };
    }
    for i in 1..=n {
        let mut current = vec![-1_000_000_000_i32; width];
        for j in i.saturating_sub(band)..=m.min(i + band) {
            let b = j + band - i;
            if j == 0 {
                current[b] = -4 * i as i32;
                directions[i * width + b] = 2;
                continue;
            }
            let mut score = previous[b]
                + if reference[i - 1] == query[j - 1] {
                    2
                } else {
                    -3
                };
            let mut direction = 1;
            if b + 1 < width && previous[b + 1] - 4 > score {
                score = previous[b + 1] - 4;
                direction = 2;
            }
            if b > 0 && current[b - 1] - 4 > score {
                score = current[b - 1] - 4;
                direction = 3;
            }
            current[b] = score;
            directions[i * width + b] = direction;
        }
        previous = current;
    }
    let (mut i, mut j) = (n, m);
    let mut ops = Vec::new();
    while i > 0 || j > 0 {
        match directions[i * width + j + band - i] {
            1 => {
                ops.push(b'M');
                i -= 1;
                j -= 1;
            }
            2 => {
                ops.push(b'D');
                i -= 1;
            }
            _ => {
                ops.push(b'I');
                j -= 1;
            }
        }
    }
    ops.reverse();
    String::from_utf8(ops).expect("ASCII alignment operations")
}
