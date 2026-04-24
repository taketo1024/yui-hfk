use std::fmt;
use yui_core::poly::MultiVar;
use yui_core::{Elem, lc::Gen};
use crate::grid_diagram::{GridDiagram, GridPoint, point::even_points};

#[derive(Clone, Debug, Default)]
pub struct GCGen {
    pub code: usize,
    pub grid_number: u8,
    pub maslov_deg: isize,
    pub alex_deg: isize,
}

impl GCGen {
    pub fn new(diagram: &GridDiagram, seq: &[u8]) -> Self {
        let pts = even_points(seq);
        let maslov_deg = diagram.maslov_degree(&pts);
        let alex_deg = diagram.alexander_degree(&pts);
        Self::from_raw(
            Self::encode(seq),
            seq.len() as u8,
            maslov_deg,
            alex_deg,
        )
    }

    pub fn from_raw(code: usize, grid_number: u8, maslov_deg: isize, alex_deg: isize) -> Self {
        GCGen { code, grid_number, maslov_deg, alex_deg }
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

impl PartialEq for GCGen {
    fn eq(&self, other: &Self) -> bool {
        self.code == other.code && self.grid_number == other.grid_number
    }
}

impl Eq for GCGen {}

impl std::hash::Hash for GCGen {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.code.hash(state);
    }
}

impl PartialOrd for GCGen {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GCGen {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.code.cmp(&other.code)
    }
}

impl fmt::Display for GCGen {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GridComplexGen({:?})", self.sequence())
    }
}

impl Elem for GCGen {
    fn math_symbol() -> String {
        String::from("GridComplexGen")
    }
}

impl Gen for GCGen {}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct GCInfGen {
    pub u: MultiVar<'U', usize>,
    pub gen: GCGen,
}

impl GCInfGen {
    pub fn new(gen: GCGen, u: MultiVar<'U', usize>) -> Self {
        GCInfGen { u, gen }
    }

    pub fn maslov_deg(&self) -> isize {
        self.gen.maslov_deg - 2 * self.u.total_deg() as isize
    }

    pub fn alex_deg(&self) -> isize {
        self.gen.alex_deg - self.u.total_deg() as isize
    }
}

impl PartialOrd for GCInfGen {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GCInfGen {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.gen.cmp(&other.gen).then_with(|| self.u.cmp(&other.u))
    }
}

impl fmt::Display for GCInfGen {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ⊗ {}", self.u, self.gen)
    }
}

impl Elem for GCInfGen {
    fn math_symbol() -> String {
        String::from("GCInfGen")
    }
}

impl Gen for GCInfGen {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let seq: Vec<u8> = vec![2, 0, 3, 1];
        let code = GCGen::encode(&seq);
        assert_eq!(GCGen::decode(code, 4), seq);
    }

    #[test]
    fn test_encode_decode() {
        // encode/decode are inverses; identity does not map to code 0
        for seq in [[0u8,1,2,3], [1,0,2,3], [2,1,0,3], [3,2,1,0]] {
            let code = GCGen::encode(&seq);
            assert_eq!(GCGen::decode(code, 4), seq);
        }
        // code 0 decodes to the last permutation, not the identity
        assert_ne!(GCGen::decode(0, 4), vec![0u8, 1, 2, 3]);
    }

    #[test]
    fn test_is_adjacent() {
        let seq1: Vec<u8> = vec![0, 1, 2, 3];
        let seq2: Vec<u8> = vec![1, 0, 2, 3]; // differs at 2 positions → adjacent
        let seq3: Vec<u8> = vec![1, 2, 0, 3]; // differs at 3 positions → not adjacent
        let g1 = GCGen::from_raw(GCGen::encode(&seq1), 4, 0, 0);
        let g2 = GCGen::from_raw(GCGen::encode(&seq2), 4, 0, 0);
        let g3 = GCGen::from_raw(GCGen::encode(&seq3), 4, 0, 0);
        assert!(g1.is_adjacent(&g2));
        assert!(!g1.is_adjacent(&g3));
    }

    #[test]
    fn test_ord_via_code() {
        let g0 = GCGen::from_raw(0, 3, 0, 0);
        let g1 = GCGen::from_raw(1, 3, 0, 0);
        assert!(g0 < g1);
    }

    #[test]
    fn test_inf_gen_degrees() {
        let gen = GCGen::from_raw(0, 3, 4, 2);
        // U_0^2 * U_1^1 has total_deg = 3
        let u = MultiVar::<'U', usize>::from([2usize, 1]);
        let t = GCInfGen::new(gen, u);
        assert_eq!(t.maslov_deg(), 4 - 2 * 3);  // -2
        assert_eq!(t.alex_deg(), 2 - 3);        // -1
    }

    #[test]
    fn test_encode_identity() {
        // The identity permutation [0,1,...,n-1] encodes to n! - 1
        for n in 1u8..=5 {
            let seq: Vec<u8> = (0..n).collect();
            let n_fact: usize = (1..=n as usize).product();
            assert_eq!(GCGen::encode(&seq), n_fact - 1);
            assert_eq!(GCGen::decode(n_fact - 1, n), seq);
        }
    }
}
