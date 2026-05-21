use num_integer::Integer;
use num_traits::Zero;
use yui_homology::{Grid2, isize2};
use crate::HFK;

use super::Region;

/// The k-th screen Sₖ of HFK(K):
///
///   Sₖ = { Uᵃx | a ∈ Z, x ∈ HFK(K), M(Uᵃx) = k }
///
/// Sₖ contains all the 'particles' projected from HFK(K) that hit M = k.
/// Each Uᵃx ∈ Sₖ has M(Uᵃx) = -2a + M(x) = k, hence a = (M(x) - k) / 2.
/// (See arXiv:2002.09210, p. 23, definition before Prop. 5.3)

pub struct Screen {
    data: Grid2<usize>,
}

impl Screen {
    pub fn from_hfk(hfk: &HFK, k: isize) -> Screen {
        let data = hfk.iter().filter_map(|(isize2(m, j), &r)| {
            if r.is_zero() || !(m - k).is_even() { return None }
            let a = (m - k) / 2;
            Some((isize2(-a, j - a), r))
        }).collect();
        Screen { data }
    }

    /// Returns true if the screen contains an element of bidegree (i, j).
    pub fn has_element_at(&self, i: isize, j: isize) -> bool {
        self.data[(i, j)] > 0
    }

    /// Returns true if some element of the screen lies inside the region r.
    pub fn has_element_in(&self, region: &Region) -> bool {
        self.data.iter().any(|(&isize2(i, j), r)| *r > 0 && region.contains((i, j)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make(triples: &[(isize, isize, usize)]) -> HFK {
        HFK::from_iter(triples.iter().map(|&(m, a, r)| (isize2(m, a), r)))
    }

    #[test]
    fn test_screen_parity_filter() {
        let h = make(&[(1, 0, 1)]);
        assert!(!Screen::from_hfk(&h, 0).has_element_at(0, 0));
        assert!(Screen::from_hfk(&h, 1).has_element_at(0, 0));
    }

    #[test]
    fn test_has_element_at() {
        let h = make(&[(0, 0, 1), (1, 1, 1), (2, 2, 1), (3, 3, 1)]);
        let s = Screen::from_hfk(&h, 0);
        assert!(s.has_element_at(0, 0));
        assert!(s.has_element_at(-1, 1));
        assert!(!s.has_element_at(0, 1));
        assert!(!s.has_element_at(-1, 2));

        let s = Screen::from_hfk(&h, 1);
        assert!(s.has_element_at(0, 1));
        assert!(s.has_element_at(-1, 2));
        assert!(!s.has_element_at(0, 0));
        assert!(!s.has_element_at(-1, 1));
    }

    #[test]
    fn test_has_element_in() {
        let h = make(&[(0, 0, 1)]);
        let s = Screen::from_hfk(&h, 0);
        assert!(s.has_element_in(&Region::simple(0, 0)));
        assert!(s.has_element_in(&Region::simple(1, 0)));
        assert!(s.has_element_in(&Region::simple(0, 1)));
        assert!(!s.has_element_in(&Region::simple(0, -1)));
    }
}
