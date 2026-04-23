use std::fmt;
use yui_core::{Elem, lc::Gen};
use crate::grid_diagram::{GridDiagram, GridPoint, point::even_points};

#[derive(Clone, Debug, Default)]
pub struct GridComplexGen {
    pub code: usize,
    pub grid_number: u8,
    pub maslov_deg: isize,
    pub alexander_deg: isize,
}

impl GridComplexGen {
    pub fn new(diagram: &GridDiagram, seq: &[u8]) -> Self {
        let pts = even_points(seq);
        let maslov_deg = diagram.maslov_degree(&pts);
        let alexander_deg = diagram.alexander_degree(&pts);
        Self::from_raw(
            Self::encode(seq),
            seq.len() as u8,
            maslov_deg,
            alexander_deg,
        )
    }

    pub fn from_raw(code: usize, grid_number: u8, maslov_deg: isize, alexander_deg: isize) -> Self {
        GridComplexGen { code, grid_number, maslov_deg, alexander_deg }
    }

    pub fn sequence(&self) -> Vec<u8> {
        Self::decode(self.code, self.grid_number)
    }

    pub fn points(&self) -> Vec<GridPoint> {
        even_points(&self.sequence())
    }

    pub fn is_adjacent(&self, other: &Self) -> bool {
        let s = self.sequence();
        let t = other.sequence();
        s.iter().zip(t.iter()).filter(|(a, b)| a != b).count() == 2
    }

    // Knuth, Volume 2, Section 3.3.2, Algorithm P
    pub fn encode(seq: &[u8]) -> usize {
        let n = seq.len();
        let mut tmp = seq.to_vec();
        let mut code = 0usize;

        for k in (0..n).rev() {
            let i = (0..=k).find(|&i| tmp[i] == k as u8).unwrap();
            code = code * (k + 1) + i;
            tmp.swap(i, k);
        }

        code
    }

    pub fn decode(code: usize, n: u8) -> Vec<u8> {
        let n = n as usize;
        let mut tmp = code;
        let mut seq: Vec<u8> = (0..n as u8).collect();

        for k in 1..n {
            let i = tmp % (k + 1);
            tmp /= k + 1;
            seq.swap(i, k);
        }

        seq
    }
}

impl PartialEq for GridComplexGen {
    fn eq(&self, other: &Self) -> bool {
        self.code == other.code && self.grid_number == other.grid_number
    }
}

impl Eq for GridComplexGen {}

impl std::hash::Hash for GridComplexGen {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.code.hash(state);
    }
}

impl PartialOrd for GridComplexGen {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GridComplexGen {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.code.cmp(&other.code)
    }
}

impl fmt::Display for GridComplexGen {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GridComplexGen({:?})", self.sequence())
    }
}

impl Elem for GridComplexGen {
    fn math_symbol() -> String {
        String::from("GridComplexGen")
    }
}

impl Gen for GridComplexGen {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let seq: Vec<u8> = vec![2, 0, 3, 1];
        let code = GridComplexGen::encode(&seq);
        assert_eq!(GridComplexGen::decode(code, 4), seq);
    }

    #[test]
    fn test_encode_decode() {
        // encode/decode are inverses; identity does not map to code 0
        for seq in [[0u8,1,2,3], [1,0,2,3], [2,1,0,3], [3,2,1,0]] {
            let code = GridComplexGen::encode(&seq);
            assert_eq!(GridComplexGen::decode(code, 4), seq);
        }
        // code 0 decodes to the last permutation, not the identity
        assert_ne!(GridComplexGen::decode(0, 4), vec![0u8, 1, 2, 3]);
    }

    #[test]
    fn test_is_adjacent() {
        let seq1: Vec<u8> = vec![0, 1, 2, 3];
        let seq2: Vec<u8> = vec![1, 0, 2, 3]; // differs at 2 positions → adjacent
        let seq3: Vec<u8> = vec![1, 2, 0, 3]; // differs at 3 positions → not adjacent
        let g1 = GridComplexGen::from_raw(GridComplexGen::encode(&seq1), 4, 0, 0);
        let g2 = GridComplexGen::from_raw(GridComplexGen::encode(&seq2), 4, 0, 0);
        let g3 = GridComplexGen::from_raw(GridComplexGen::encode(&seq3), 4, 0, 0);
        assert!(g1.is_adjacent(&g2));
        assert!(!g1.is_adjacent(&g3));
    }

    #[test]
    fn test_ord_via_code() {
        let g0 = GridComplexGen::from_raw(0, 3, 0, 0);
        let g1 = GridComplexGen::from_raw(1, 3, 0, 0);
        assert!(g0 < g1);
    }

    #[test]
    fn test_encode_identity() {
        // The identity permutation [0,1,...,n-1] encodes to n! - 1
        for n in 1u8..=5 {
            let seq: Vec<u8> = (0..n).collect();
            let n_fact: usize = (1..=n as usize).product();
            assert_eq!(GridComplexGen::encode(&seq), n_fact - 1);
            assert_eq!(GridComplexGen::decode(n_fact - 1, n), seq);
        }
    }
}
