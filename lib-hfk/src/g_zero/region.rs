use std::fmt;
use std::str::FromStr;
use num_traits::Zero;
use yui_core::lc::Lc;
use yui_core::num::FF2;
use crate::grid_complex::GCInfGen;

/// A single corner (i, j) representing the rectangle (-∞, i] × (-∞, j].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Corner {
    pub i: isize,
    pub j: isize,
}

impl Corner {
    pub fn new(i: isize, j: isize) -> Self {
        Corner { i, j }
    }

    /// Returns true if point p = (x, y) satisfies x ≤ i and y ≤ j.
    pub fn contains(&self, p: (isize, isize)) -> bool {
        p.0 <= self.i && p.1 <= self.j
    }

    /// Returns true if other.i ≤ self.i and other.j ≤ self.j (self dominates other).
    pub fn includes(&self, other: Corner) -> bool {
        other.i <= self.i && other.j <= self.j
    }

    pub fn transposed(&self) -> Self {
        Corner::new(self.j, self.i)
    }

    pub fn shift(&self, di: isize, dj: isize) -> Self {
        Corner::new(self.i + di, self.j + dj)
    }
}

impl PartialOrd for Corner {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        use std::cmp::Ordering::*;
        match (self.includes(*other), other.includes(*self)) {
            (true,  true)  => Some(Equal),
            (true,  false) => Some(Greater),
            (false, true)  => Some(Less),
            (false, false) => None,
        }
    }
}

impl fmt::Display for Corner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", coord(self.i), coord(self.j))
    }
}

// note: isize::MAX and isize::MIN stand for ±∞, as in the regions of tau.
fn coord(x: isize) -> String {
    match x {
        isize::MAX => "∞".to_string(),
        isize::MIN => "-∞".to_string(),
        _ => x.to_string(),
    }
}

/// A closed region in the grid, represented as a minimal antichain of corners.
///
/// The region is the union of rectangles (-∞, i] × (-∞, j] for each corner (i, j).
/// The corners are stored in reduced (antichain) form, sorted by (i, j).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Region {
    corners: Vec<Corner>, // reduced antichain, sorted ascending by (i, j)
}

impl Region {
    /// Build from any iterable of corners, reducing to canonical antichain form.
    pub fn new(corners: impl IntoIterator<Item = Corner>) -> Self {
        Region { corners: Self::reduce(corners) }
    }

    pub fn empty() -> Self {
        Region { corners: vec![] }
    }

    pub fn simple(i: isize, j: isize) -> Self {
        Region { corners: vec![Corner::new(i, j)] }
    }

    pub fn corners(&self) -> &[Corner] {
        &self.corners
    }

    pub fn is_empty(&self) -> bool {
        self.corners.is_empty()
    }

    pub fn is_simple(&self) -> bool {
        self.corners.len() == 1
    }

    /// Maximum i-coordinate among corners (used as the "shift number").
    pub fn right_end(&self) -> isize {
        self.corners.iter().map(|c| c.i).max().unwrap_or(0)
    }

    /// Returns true if p is in the region (some corner dominates p).
    pub fn contains(&self, p: (isize, isize)) -> bool {
        self.corners.iter().any(|c| c.contains(p))
    }

    pub fn contains_gen(&self, x: &GCInfGen) -> bool {
        self.contains((x.alg_deg(), x.alex_deg()))
    }

    pub fn contains_chain(&self, z: &Lc<GCInfGen, FF2>) -> bool {
        z.iter().all(|(g, c)| c.is_zero() || self.contains_gen(g))
    }

    /// Returns true if every corner of `other` is contained in self.
    /// i.e., self is a superset of the region defined by `other`.
    pub fn includes(&self, other: &Region) -> bool {
        other.corners.iter().all(|c| self.contains((c.i, c.j)))
    }

    pub fn transposed(&self) -> Region {
        Region::new(self.corners.iter().map(|c| c.transposed()))
    }

    pub fn is_symmetric(&self) -> bool {
        self.corners.iter().all(|c| {
            let t = c.transposed();
            self.corners.binary_search_by_key(&(t.i, t.j), |x| (x.i, x.j)).is_ok()
        })
    }

    pub fn shift(&self, di: isize, dj: isize) -> Region {
        // Already in canonical form after shifting (shift preserves the antichain structure).
        Region { corners: self.corners.iter().map(|c| c.shift(di, dj)).collect() }
    }

    pub fn union(&self, other: &Region) -> Region {
        Region::new(self.corners.iter().chain(other.corners.iter()).copied())
    }

    pub fn add_corner(&self, c: Corner) -> Region {
        Region::new(self.corners.iter().copied().chain(std::iter::once(c)))
    }

    /// The corner with the smallest i > c.i, if any.
    pub fn corner_next_to(&self, c: Corner) -> Option<&Corner> {
        self.corners.iter()
            .filter(|&&x| x.i > c.i)
            .min_by_key(|x| x.i)
    }

    // Drop terms of `z` whose (alg_deg, alex_deg) is inside self.
    pub fn cutoff(&self, z: &Lc<GCInfGen, FF2>) -> Lc<GCInfGen, FF2> {
        z.filtered(|g| !self.contains_gen(g))
    }

    // ---- internal ----

    fn reduce(corners: impl IntoIterator<Item = Corner>) -> Vec<Corner> {
        let mut result: Vec<Corner> = Vec::new();
        for c in corners {
            if result.iter().any(|r| r.includes(c)) {
                continue;
            }
            result.retain(|r| !c.includes(*r));
            result.push(c);
        }
        result.sort_by_key(|c| (c.i, c.j));
        result
    }
}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        use std::cmp::Ordering::*;
        match (other.includes(self), self.includes(other)) {
            (true,  true)  => Some(Equal),
            (true,  false) => Some(Less),
            (false, true)  => Some(Greater),
            (false, false) => None,
        }
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.corners.len() == 1 {
            write!(f, "R{}", self.corners[0])
        } else {
            let cs: Vec<_> = self.corners.iter().map(|c| c.to_string()).collect();
            write!(f, "R[{}]", cs.join(", "))
        }
    }
}

/// Parses the `Display` format, e.g. `R(1, 1)` or `R[(0, 1), (1, 0)]`. The leading `R` is optional.
impl FromStr for Region {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let json = s.trim().trim_start_matches('R').replace('(', "[").replace(')', "]");
        let corners: Vec<[isize; 2]> = serde_json::from_str(&json)
            .or_else(|_| serde_json::from_str::<[isize; 2]>(&json).map(|c| vec![c]))
            .map_err(|_| format!("invalid region: {s}"))?;
        Ok(Region::new(corners.into_iter().map(|[i, j]| Corner::new(i, j))))
    }
}

#[cfg(test)]
mod tests {
    use crate::grid_complex::chain::tests::{lc_of, make_gen};

    use super::*;

    #[test]
    fn test_corner_includes() {
        let c = Corner::new(2, 3);
        assert!(c.includes(Corner::new(1, 2)));
        assert!(c.includes(Corner::new(2, 3)));
        assert!(!c.includes(Corner::new(3, 2)));
    }

    #[test]
    fn test_region_reduce() {
        // (2, 3) dominates (1, 2), so only (2, 3) survives.
        let r = Region::new([Corner::new(1, 2), Corner::new(2, 3)]);
        assert_eq!(r.corners(), &[Corner::new(2, 3)]);
    }

    #[test]
    fn test_region_includes() {
        let r1 = Region::simple(1, 2);
        let r2 = Region::simple(2, 3);
        assert!(r2.includes(&r1));
        assert!(!r1.includes(&r2));
    }

    #[test]
    fn test_region_from_str() {
        let r1 = Region::simple(1, 1);
        let r2 = Region::new([Corner::new(0, 1), Corner::new(1, 0)]);

        assert_eq!(r1.to_string().parse::<Region>(), Ok(r1.clone()));
        assert_eq!(r2.to_string().parse::<Region>(), Ok(r2.clone()));
        assert_eq!("(1,1)".parse::<Region>(), Ok(r1));
        assert_eq!("[(0, 1), (1, 0)]".parse::<Region>(), Ok(r2));
        assert_eq!("R[(-2, 1), (-1, 0), (1, -1)]".parse::<Region>().map(|r| r.corners().len()), Ok(3));
        assert!("(1, 1, 1)".parse::<Region>().is_err());
        assert!("foo".parse::<Region>().is_err());
    }

    #[test]
    fn test_region_display() {
        let r1 = Region::new([Corner::new(-1, isize::MAX), Corner::new(0, 0)]);
        let r2 = Region::new([Corner::new(isize::MIN, 2)]);
        assert_eq!(r1.to_string(), "R[(-1, ∞), (0, 0)]");
        assert_eq!(r2.to_string(), "R(-∞, 2)");
    }

    #[test]
    fn test_partial_ord() {
        let r1 = Region::simple(0, 1);
        let r2 = Region::simple(1, 2);
        let r3 = Region::simple(1, 0); // incomparable to r1: (1,0) vs (0,1)
        assert!(r1 < r2);
        assert!(r3 < r2);
        // r1=(0,1) and r3=(1,0): neither dominates the other
        assert_eq!(r1.partial_cmp(&r3), None);
    }

    // ---- cutoff tests ----

    #[test]
    fn test_cutoff_empty() {
        let r = Region::simple(0, 0);
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert!(r.cutoff(&z).is_zero());
    }

    #[test]
    fn test_cutoff_drops_inside_keeps_outside() {
        // g_in:  (alg=0, alex=0)  → inside  R(0,0)  → dropped
        // g_out: (alg=0, alex=1)  → outside R(0,0) (alex 1 > 0) → kept
        let g_in  = make_gen(0, 2, 0, 0, &[]);  // alg=0, alex=0
        let g_out = make_gen(1, 2, 0, 1, &[]);  // alg=0, alex=1
        let lc = lc_of(&[g_in, g_out.clone()]);
        let r = Region::simple(0, 0);
        let result = r.cutoff(&lc);
        // Only g_out should remain
        assert_eq!(result.nterms(), 1);
        assert!(!result.coeff(&g_out).is_zero());
    }

    // A term with no U₀ factor sits at alg = 0 however large its other
    // exponents are. Under the old (total-degree) convention it would have been
    // pushed to alg = -3 and wrongly counted as inside R(-1, ∞).
    #[test]
    fn test_contains_gen_ignores_higher_vars() {
        let g = make_gen(0, 4, 0, 3, &[0, 3]); // U₁³: alg = 0, alex = 3 - 3 = 0

        assert_eq!(g.alg_deg(), 0);
        assert_eq!(g.alex_deg(), 0);

        // alg = 0 is outside any region whose corners all have i ≤ -1
        let r = Region::simple(-1, isize::MAX);
        assert!(!r.contains_gen(&g));

        // but inside one that reaches alg = 0
        assert!(Region::simple(0, 0).contains_gen(&g));
    }

    #[test]
    fn test_cutoff_keeps_mixed_monomial_outside() {
        let g_u0 = make_gen(0, 4, 0, 3, &[3]);    // U₀³: alg = -3 → inside
        let g_u1 = make_gen(1, 4, 0, 3, &[0, 3]); // U₁³: alg =  0 → outside
        let lc = lc_of(&[g_u0, g_u1.clone()]);

        let r = Region::simple(-1, isize::MAX);
        let result = r.cutoff(&lc);

        assert_eq!(result.nterms(), 1);
        assert!(!result.coeff(&g_u1).is_zero());
    }

    #[test]
    fn test_cutoff_all_outside() {
        // U^1 generator: alg = -1, alex = alex_raw - 1
        // With alex_raw=2: alex = 1, coord = (-1, 1) → outside R(0,0)? -1 ≤ 0 and 1 > 0 → outside
        let g = make_gen(0, 2, 0, 2, &[1]); // alg=-1, alex=2-1=1
        let lc = lc_of(&[g.clone()]);
        let r = Region::simple(0, 0);
        let result = r.cutoff(&lc);
        assert_eq!(result.nterms(), 1);
        assert!(!result.coeff(&g).is_zero());
    }
}
