//! Exact canonical circular-word gate. Alignment and clustering order stay in Python.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

// A leading sentinel and three bits/base distinguish all ACGTN words of length <=9.
fn encode(word: &str) -> Option<u32> {
    if word.is_empty() || word.len() > 9 {
        return None;
    }
    word.bytes().try_fold(1_u32, |value, base| {
        let digit = match base {
            b'A' => 1,
            b'C' => 2,
            b'G' => 3,
            b'T' => 4,
            b'N' => 5,
            _ => return None,
        };
        Some((value << 3) | digit)
    })
}

fn normalized_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'C' | b'G' | b'N' | b'T' => Some(base),
        _ => None,
    }
}

fn complement(base: u8) -> Option<u8> {
    match normalized_base(base)? {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'N' => Some(b'N'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

fn encoded_digit(base: u8) -> Option<u32> {
    match normalized_base(base)? {
        b'A' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        b'T' => Some(4),
        b'N' => Some(5),
        _ => None,
    }
}

fn canonical_circular_words(sequence: &str) -> Result<Vec<(u32, u32)>, &'static str> {
    let bytes = sequence.as_bytes();
    if bytes.is_empty() {
        return Err("Representative sequence must be nonempty");
    }
    if bytes.iter().any(|base| normalized_base(*base).is_none()) {
        return Err("Representative sequence must contain uppercase ACGTN bases");
    }
    let k = bytes.len().min(9);
    let mut counts: HashMap<u32, u32> = HashMap::with_capacity(bytes.len().min(262_144));
    for start in 0..bytes.len() {
        let mut orientation = Ordering::Equal;
        for offset in 0..k {
            let forward = bytes[(start + offset) % bytes.len()];
            let reverse = complement(bytes[(start + k - 1 - offset) % bytes.len()]).unwrap();
            orientation = forward.cmp(&reverse);
            if orientation != Ordering::Equal {
                break;
            }
        }
        let use_reverse = orientation == Ordering::Greater;
        let mut encoded = 1_u32;
        for offset in 0..k {
            let base = if use_reverse {
                complement(bytes[(start + k - 1 - offset) % bytes.len()]).unwrap()
            } else {
                bytes[(start + offset) % bytes.len()]
            };
            encoded = (encoded << 3) | encoded_digit(base).unwrap();
        }
        *counts.entry(encoded).or_default() += 1;
    }
    Ok(counts.into_iter().collect())
}

fn validate_words(length: u32, words: Vec<(String, u32)>) -> Result<Vec<(u32, u32)>, &'static str> {
    if length == 0 {
        return Err("Representative length must be positive");
    }
    let mut seen = HashSet::with_capacity(words.len());
    let mut encoded = Vec::with_capacity(words.len());
    let mut total = 0_u64;
    for (word, count) in words {
        let key = encode(&word).ok_or("Index words must contain 1-9 uppercase ACGTN bases")?;
        if word.len() != length.min(9) as usize || count == 0 || !seen.insert(key) {
            return Err("Index words must be unique, with expected length and positive counts");
        }
        total += u64::from(count);
        encoded.push((key, count));
    }
    if total != u64::from(length) {
        return Err("Circular-word multiplicities must sum to sequence length");
    }
    Ok(encoded)
}

#[pyclass]
#[derive(Default)]
pub struct RepresentativeIndex {
    postings: HashMap<u32, Vec<(u32, u32)>>,
    lengths: Vec<u32>,
    length_ids: HashMap<u32, usize>,
    representative_buckets: Vec<usize>,
    overlaps: Vec<u32>,
}

#[pymethods]
impl RepresentativeIndex {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn append(&mut self, length: u32, words: Vec<(String, u32)>) -> PyResult<u32> {
        // Validate the complete append before changing any index state.
        let encoded = validate_words(length, words).map_err(PyValueError::new_err)?;
        self.append_encoded(length, encoded)
    }

    pub fn append_sequence(&mut self, sequence: String) -> PyResult<u32> {
        let length = u32::try_from(sequence.len())
            .map_err(|_| PyValueError::new_err("Representative length exceeds unsigned 32 bits"))?;
        let encoded = canonical_circular_words(&sequence).map_err(PyValueError::new_err)?;
        self.append_encoded(length, encoded)
    }

    pub fn candidates_sequence(
        &mut self,
        sequence: String,
        minimum_identity: f64,
    ) -> PyResult<Vec<u32>> {
        let length = u32::try_from(sequence.len())
            .map_err(|_| PyValueError::new_err("Representative length exceeds unsigned 32 bits"))?;
        let encoded = canonical_circular_words(&sequence).map_err(PyValueError::new_err)?;
        self.candidates_encoded(length, encoded, minimum_identity)
    }

    pub fn candidates(
        &mut self,
        length: u32,
        words: Vec<(String, u32)>,
        minimum_identity: f64,
    ) -> PyResult<Vec<u32>> {
        let encoded = validate_words(length, words).map_err(PyValueError::new_err)?;
        self.candidates_encoded(length, encoded, minimum_identity)
    }
}

impl RepresentativeIndex {
    fn append_encoded(&mut self, length: u32, encoded: Vec<(u32, u32)>) -> PyResult<u32> {
        let id = u32::try_from(self.overlaps.len())
            .map_err(|_| PyValueError::new_err("Representative ID exceeds unsigned 32 bits"))?;
        let bucket = *self.length_ids.entry(length).or_insert_with(|| {
            let next = self.lengths.len();
            self.lengths.push(length);
            next
        });
        for (word, count) in encoded {
            self.postings.entry(word).or_default().push((id, count));
        }
        self.representative_buckets.push(bucket);
        self.overlaps.push(0);
        Ok(id)
    }

    fn candidates_encoded(
        &mut self,
        length: u32,
        encoded: Vec<(u32, u32)>,
        minimum_identity: f64,
    ) -> PyResult<Vec<u32>> {
        if !minimum_identity.is_finite() || minimum_identity <= 0.0 || minimum_identity > 1.0 {
            return Err(PyValueError::new_err("Cluster identity must be in (0,1]"));
        }
        if minimum_identity <= 0.9 || length < 20 {
            return Ok((0..self.overlaps.len()).map(|i| i as u32).collect());
        }
        let thresholds: Vec<Option<u32>> = self
            .lengths
            .iter()
            .map(|&other| {
                let size = length.max(other);
                // Keep Python's operation order, tolerance and floor exactly.
                let limit = ((1.0 - minimum_identity) * f64::from(size) + 1e-9).floor() as u32;
                (length.abs_diff(other) <= limit)
                    .then(|| size.saturating_sub(length.min(other).min(9) * limit))
            })
            .collect();
        let mut touched = Vec::new();
        for (word, count) in encoded {
            if let Some(postings) = self.postings.get(&word) {
                for &(id, other_count) in postings {
                    let j = id as usize;
                    if thresholds[self.representative_buckets[j]].is_some() {
                        if self.overlaps[j] == 0 {
                            touched.push(id);
                        }
                        // Unique validated words bound each sum by the query length.
                        self.overlaps[j] += count.min(other_count);
                    }
                }
            }
        }
        let mut result = Vec::new();
        for id in touched {
            let j = id as usize;
            if self.overlaps[j] >= thresholds[self.representative_buckets[j]].unwrap() {
                result.push(id);
            }
            self.overlaps[j] = 0;
        }
        result.sort_unstable();
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_is_injective_including_length_and_ambiguity() {
        let mut words = vec![String::new()];
        let mut codes = HashSet::new();
        for _ in 0..6 {
            words = words
                .iter()
                .flat_map(|w| "ACGTN".chars().map(move |b| format!("{w}{b}")))
                .collect();
            for word in &words {
                assert!(codes.insert(encode(word).unwrap()));
            }
        }
        assert_eq!(encode("a"), None);
        assert_eq!(encode("ACGTACGTAC"), None);
    }

    #[test]
    fn gate_preserves_multiplicities_length_exclusion_and_repeated_queries() {
        let words = vec![("AAAAAAAAA".into(), 20)];
        let mut index = RepresentativeIndex::new();
        index.append(20, words.clone()).unwrap();
        index.append(40, vec![("AAAAAAAAA".into(), 40)]).unwrap();
        index.append(20, vec![("CCCCCCCCC".into(), 20)]).unwrap();
        for _ in 0..2 {
            assert_eq!(index.candidates(20, words.clone(), 0.95).unwrap(), vec![0]);
            assert_eq!(
                index.candidates(20, words.clone(), 0.9).unwrap(),
                vec![0, 1, 2]
            );
        }
    }

    #[test]
    fn malformed_word_banks_are_rejected_before_append() {
        assert!(validate_words(0, vec![]).is_err());
        assert!(validate_words(20, vec![("AAAA".into(), 20)]).is_err());
        assert!(
            validate_words(20, vec![("AAAAAAAAA".into(), 10), ("AAAAAAAAA".into(), 10)]).is_err()
        );
        assert!(validate_words(20, vec![("AAAAAAAAA".into(), 19)]).is_err());
    }

    #[test]
    fn native_circular_words_match_python_orientation_and_multiplicity_rules() {
        let expected = validate_words(
            10,
            vec![
                ("AACGTNACG".into(), 1),
                ("ACGTAACGT".into(), 1),
                ("ACGTNACGT".into(), 1),
                ("CGTAACGTN".into(), 1),
                ("CGTNACGTA".into(), 1),
                ("CGTTACGTN".into(), 1),
                ("GTAACGTNA".into(), 1),
                ("GTNACGTAA".into(), 1),
                ("GTNACGTTA".into(), 1),
                ("GTTACGTNA".into(), 1),
            ],
        )
        .unwrap();
        let observed = canonical_circular_words("AACGTNACGT").unwrap();
        assert_eq!(
            observed.into_iter().collect::<HashMap<_, _>>(),
            expected.into_iter().collect()
        );
        assert_eq!(
            canonical_circular_words("AAAA").unwrap(),
            vec![(encode("AAAA").unwrap(), 4)]
        );
        assert!(canonical_circular_words("").is_err());
        assert!(canonical_circular_words("acgt").is_err());
        assert!(canonical_circular_words("ACG?").is_err());
    }

    #[test]
    fn sequence_methods_match_validated_word_methods() {
        let sequence = "ACGTNACGTACGTNACGTAC";
        let encoded = canonical_circular_words(sequence).unwrap();
        let words: Vec<(String, u32)> = encoded
            .iter()
            .map(|(code, count)| {
                let word = (0..9)
                    .rev()
                    .map(|offset| match (code >> (3 * offset)) & 7 {
                        1 => 'A',
                        2 => 'C',
                        3 => 'G',
                        4 => 'T',
                        5 => 'N',
                        _ => unreachable!(),
                    })
                    .collect();
                (word, *count)
            })
            .collect();
        let mut sequence_index = RepresentativeIndex::new();
        let mut word_index = RepresentativeIndex::new();
        assert_eq!(sequence_index.append_sequence(sequence.into()).unwrap(), 0);
        assert_eq!(
            word_index
                .append(sequence.len() as u32, words.clone())
                .unwrap(),
            0
        );
        assert_eq!(
            sequence_index
                .candidates_sequence(sequence.into(), 0.95)
                .unwrap(),
            word_index
                .candidates(sequence.len() as u32, words, 0.95)
                .unwrap()
        );
    }
}
