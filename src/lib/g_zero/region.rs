use std::fmt;

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
        write!(f, "({}, {})", self.i, self.j)
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

    /// Remove corner `c` and replace it with (c.i-1, c.j) and (c.i, c.j-1).
    pub fn cut_off(&self, c: Corner) -> Region {
        let mut cs: Vec<Corner> = self.corners.iter().filter(|&&x| x != c).copied().collect();
        cs.push(Corner::new(c.i - 1, c.j));
        cs.push(Corner::new(c.i, c.j - 1));
        Region::new(cs)
    }

    /// The corner with the smallest i > c.i, if any.
    pub fn corner_next_to(&self, c: Corner) -> Option<&Corner> {
        self.corners.iter()
            .filter(|&&x| x.i > c.i)
            .min_by_key(|x| x.i)
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

#[cfg(test)]
mod tests {
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
    fn test_partial_ord() {
        let r1 = Region::simple(0, 1);
        let r2 = Region::simple(1, 2);
        let r3 = Region::simple(1, 0); // incomparable to r1: (1,0) vs (0,1)
        assert!(r1 < r2);
        assert!(r3 < r2);
        // r1=(0,1) and r3=(1,0): neither dominates the other
        assert_eq!(r1.partial_cmp(&r3), None);
    }

    #[test]
    fn test_cut_off() {
        let r = Region::new([Corner::new(1, 2), Corner::new(0, 3)]);
        let cut = r.cut_off(Corner::new(1, 2));
        // Replaces (1,2) with (0,2) and (1,1); then reduces with (0,3)
        assert!(cut.corners().contains(&Corner::new(0, 3)));
    }

}
