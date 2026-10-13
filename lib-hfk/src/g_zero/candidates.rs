use crate::HFK;
use super::{Corner, Region, Screen};

pub struct Candidates {
    pub regions: Vec<Region>,
    pub accepted: Vec<Region>,
    pub resigned: Vec<Region>,
}

impl Candidates {
    pub fn new() -> Self { 
        Candidates { regions: vec![], accepted: vec![], resigned: vec![] }
    }

    /// Generates all combinatorial candidates for genus g.
    ///
    /// A candidate is a closed region R such that R includes the minimal region
    ///   M_g = { (-g+k, -k) | k = 0, ..., g }
    /// and every corner of R lies in the diamond |i+j| ≤ g, |i-j| ≤ g.
    ///
    /// ```text
    ///                    ^ j
    ///                    |
    ///                   /|\
    ///        j - i = g / | \ i + j = g
    ///                 /  |  \
    ///                /   |   \          i
    ///      --------------+-------------->
    ///                \   |   /
    ///                 \  |  /
    ///       i + j = -g \ | / i - j = g
    ///                   \|/
    ///                    |
    /// ```
    pub fn generate_from(g: usize) -> Self {
        let g = g as isize;
        let minimal = Region::new((0..=g).map(|k| Corner::new(-g + k, -k)));

        fn extend(r: Region, minimal: &Region, g: isize, out: &mut Vec<Region>) {
            if r.includes(minimal) {
                out.push(r.clone());
            }

            let last = *r.corners().last().unwrap();

            // Each new corner must be strictly right (larger i) and strictly below
            // (smaller j) the last corner, and must lie inside the diamond:
            //   |i+j| ≤ g  and  |i-j| ≤ g,  but NOT on the upper-right boundary i+j = g.
            for j1 in -g..last.j {
                for i1 in last.i + 1..=g {
                    if i1 + j1 >= -g && i1 + j1 < g && (i1 - j1).abs() <= g {
                        extend(r.add_corner(Corner::new(i1, j1)), minimal, g, out);
                    }
                }
            }
        }

        let mut regions = Vec::new();

        // The seed starts at j = 0, not j = -g: the minimal region has corner (-g, 0), 
        // so any valid candidate needs a corner with j ≥ 0. Since extend only adds corners
        // with strictly smaller j, the seed always has the maximum j in the region.
        // A seed with j < 0 therefore cannot include the minimal region.
        for j in 0..=g {
            // run i over the horizontal segment intersecting the diamond. 
            for i in (j - g)..=(g - j) { 
                extend(Region::simple(i, j), &minimal, g, &mut regions);
            }
        }

        Candidates { regions, accepted: vec![], resigned: vec![] }
    }

    /// Removes symmetric pairs {R, Rᵀ} in place, keeping the one with the
    /// smaller shift number (max i-coordinate).
    pub fn desymmetrize(&mut self) {
        let mut result: Vec<Region> = Vec::new();
        for r in self.regions.drain(..) {
            let rt = r.transposed();
            if let Some(pos) = result.iter().position(|s| *s == rt) {
                if r.right_end() < result[pos].right_end() {
                    result[pos] = r;
                }
            } else {
                result.push(r);
            }
        }
        self.regions = result;
    }

    pub fn add(&mut self, r: Region) { 
        self.regions.push(r);
    }

    /// Marks `r` as realizable: moves `r` and `rᵀ` into `accepted`, removes any
    /// previously accepted regions larger than them (non-minimal), and removes all
    /// candidates larger than or equal to `r` or `rᵀ` from `regions`.
    pub fn accept(&mut self, r: &Region) {
        let rt = r.transposed();

        // Keep only minimal accepted regions: drop those that r or rᵀ now dominate.
        self.accepted.retain(|s| !s.includes(r) && !s.includes(&rt));
        self.accepted.push(r.clone());
        if rt != *r {
            self.accepted.push(rt.clone());
        }

        // Remove r, rᵀ, and all larger regions from the remaining candidates.
        self.regions.retain(|s| !s.includes(r) && !s.includes(&rt));
        
        // Sort accepted regions.
        Self::sort_regions(&mut self.accepted);
    }

    /// Marks `r` as unrealizable: removes all candidates equal to or smaller than
    /// `r` or `rᵀ` from `regions` (they are automatically unrealizable).
    pub fn reject(&mut self, r: &Region) {
        let rt = r.transposed();
        self.regions.retain(|s| !r.includes(s) && !rt.includes(s));
    }

    /// Marks `r` as undecided: moves `r` and `rᵀ` into `resigned` and removes
    /// them from `regions`.
    pub fn resign(&mut self, r: &Region) {
        let rt = r.transposed();

        self.resigned.push(r.clone());
        if rt != *r {
            self.resigned.push(rt.clone());
        }

        self.regions.retain(|s| s != r && s != &rt);
        Self::sort_regions(&mut self.resigned);
    }

    /// Sorts candidates primarily by `right_end` (ascending), then by `Region::partial_cmp`
    /// (smaller regions first; incomparable regions are treated as equal).
    pub fn sort(&mut self) {
        Self::sort_regions(&mut self.regions);
    }

    fn sort_regions(v: &mut Vec<Region>) {
        v.sort_by(|a, b| {
            a.right_end().cmp(&b.right_end())
                .then_with(|| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        });
    }

    /// Filters candidates in place using the 4-genus g₄.
    ///
    /// Two conditions must hold for each candidate R:
    ///   1. Every corner (i, j) satisfies |i + j| ≤ g₄.
    ///   2. If R has more than one corner, no corner lies on the boundary i + j = g₄.
    pub fn filter_by_g4(&mut self, g4: usize) {
        let g4 = g4 as isize;
        self.regions.retain(|r| {
            r.corners().iter().all(|c| (c.i + c.j).abs() <= g4)
            && (r.is_simple() || r.corners().iter().all(|c| c.i + c.j < g4))
        });
    }

    pub fn filter_by_tau(&mut self, tau: isize) {
        let accept = Region::new([Corner::new(-1, isize::MAX), Corner::new(0, tau)]);
        let reject = Region::new([Corner::new(-1, isize::MAX), Corner::new(0, tau - 1)]);
        self.accept(&accept);
        self.reject(&reject);
    }

    /// Filters candidates in place by the HFK structure via Proposition 5.3
    /// (arXiv:2002.09210).
    pub fn filter_by_hfk(&mut self, hfk: &HFK) {
        let s0 = Screen::from_hfk(hfk, 0);
        let s1 = Screen::from_hfk(hfk, -1);

        // Prop 5.3, Condition 1: every corner must appear in S₀.
        self.regions.retain(|r| 
            r.corners().iter().all(|c| s0.has_element_at(c.i, c.j))
        );

        // Prop 5.3, Condition 2: for each adjacent corner pair (c0, c1) with
        // c0.i < c1.i, S₋₁ must have element in the simple region R(c0.i, c1.j).
        self.regions.retain(|r| 
            r.corners().iter().all(|c0| {
                r.corner_next_to(*c0)
                    .map(|c1| s1.has_element_in(&Region::simple(c0.i, c1.j)))
                    .unwrap_or(true)
            })
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use yui_homology::isize2;
    use crate::GridDiagram;

    fn hfk_from(triples: &[(isize, isize, usize)]) -> HFK {
        HFK::from_iter(triples.iter().map(|&(m, a, r)| (isize2(m, a), r)))
    }

    fn make_candidates(regions: &[Region]) -> Candidates {
        Candidates { regions: regions.to_vec(), accepted: vec![], resigned: vec![] }
    }

    #[test]
    fn test_accept_moves_to_accepted_and_prunes_larger() {
        // regions = [R(0,0), R(0,1), R(-1,0)∪R(0,-1)]
        // R(0,1) ⊇ R(0,0), so accepting R(0,0) removes R(0,1) from regions.
        // R(-1,0)∪R(0,-1) does not include R(0,0), so it stays.
        let multi = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);
        let mut cs = make_candidates(&[Region::simple(0, 0), Region::simple(0, 1), multi.clone()]);
        cs.accept(&Region::simple(0, 0));
        assert_eq!(cs.accepted, vec![Region::simple(0, 0)]);
        assert!(cs.regions.contains(&multi));
        assert!(!cs.regions.contains(&Region::simple(0, 0)));
        assert!(!cs.regions.contains(&Region::simple(0, 1)));
    }

    #[test]
    fn test_accept_prunes_non_minimal_from_accepted() {
        // If R(0,1) was already accepted and we now accept R(0,0),
        // R(0,1) must be removed from accepted (it is larger than R(0,0)).
        let mut cs = make_candidates(&[]);
        cs.accepted = vec![Region::simple(0, 1)];
        cs.accept(&Region::simple(0, 0));
        assert!(!cs.accepted.contains(&Region::simple(0, 1)));
        assert!(cs.accepted.contains(&Region::simple(0, 0)));
    }

    #[test]
    fn test_reject_removes_equal_and_smaller() {
        // Rejecting R(-1,0)∪R(0,-1) removes it and all subregions.
        // R(0,1) and R(0,0) are NOT subregions, so they stay.
        let multi = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);
        let mut cs = make_candidates(&[Region::simple(0, 1), Region::simple(0, 0), multi.clone()]);
        cs.reject(&multi);
        assert!(!cs.regions.contains(&multi));
        assert!(cs.regions.contains(&Region::simple(0, 1)));
        assert!(cs.regions.contains(&Region::simple(0, 0)));
    }

    #[test]
    fn test_filter_by_g4_removes_corners_outside_slice() {
        // genus=1 candidates after dedup: R(-1,0)∪R(0,-1), R(0,0), R(0,1)
        // With g4=0: only corners with |i+j| ≤ 0 (i.e. i+j=0) survive.
        // R(0,0): |0|=0 ✓. R(-1,0)∪R(0,-1): |-1|=1 > 0. R(0,1): |1|=1 > 0.
        let mut cs = Candidates::generate_from(1);
        cs.desymmetrize();
        cs.filter_by_g4(0);
        assert_eq!(cs.regions, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_filter_by_g4_removes_multi_corner_on_boundary() {
        // Manually construct a multi-corner region with one corner on i+j = g4.
        // R[(0,1),(1,0)] has corners with i+j = 1 and i+j = 1 (both on boundary g4=1).
        // As a multi-corner region it must be removed; R(0,1) alone survives.
        let r_multi = Region::new([Corner::new(0, 1), Corner::new(1, 0)]);
        let mut cs = Candidates { regions: vec![Region::simple(0, 1), r_multi], accepted: vec![], resigned: vec![] };
        cs.filter_by_g4(1);
        assert_eq!(cs.regions, vec![Region::simple(0, 1)]);
    }

    #[test]
    fn test_generate_genus0() {
        let cs = Candidates::generate_from(0);
        assert_eq!(cs.regions, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_generate_genus1_before_dedup() {
        let cs = Candidates::generate_from(1);
        // Before dedup both R(0,1) and R(1,0) are present.
        assert!(cs.regions.contains(&Region::simple(0, 1)));
        assert!(cs.regions.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_deduplicate_genus1() {
        let mut cs = Candidates::generate_from(1);
        cs.desymmetrize();
        // Three candidates after dedup: R(-1,0)∪R(0,-1), R(0,0), R(0,1)
        assert_eq!(cs.regions.len(), 3);
        assert!(cs.regions.contains(&Region::new([Corner::new(-1, 0), Corner::new(0, -1)])));
        assert!(cs.regions.contains(&Region::simple(0, 0)));
        assert!(cs.regions.contains(&Region::simple(0, 1)));
        // R(1,0) replaced by R(0,1) (smaller shift number)
        assert!(!cs.regions.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_candidates_unknot() {
        let gd = GridDiagram::unknot();
        let hfk = HFK::from(&gd);
        let mut cs = Candidates::generate_from(hfk.genus() as usize);
        cs.desymmetrize();
        cs.filter_by_hfk(&hfk);
        assert_eq!(cs.regions, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_candidates_trefoil() {
        // ĤFK(trefoil) = {(0,1):1, (-1,0):1, (-2,-1):1}
        // S₀ = {(0,1),(1,0)}, S₋₁ = {(0,0)}
        // Only R(0,1) passes both conditions.
        let hfk = hfk_from(&[(0, 1, 1), (-1, 0, 1), (-2, -1, 1)]);
        let mut cs = Candidates::generate_from(hfk.genus() as usize);
        cs.desymmetrize();
        cs.filter_by_hfk(&hfk);
        assert_eq!(cs.regions, vec![Region::simple(0, 1)]);
    }

    #[test]
    fn test_candidates_figure8() {
        // ĤFK(4₁) = {(1,1):1, (0,0):3, (-1,-1):1}
        // S₀ = {(0,0)}, so only corners at (0,0) survive.
        let hfk = hfk_from(&[(1, 1, 1), (0, 0, 3), (-1, -1, 1)]);
        let mut cs = Candidates::generate_from(hfk.genus() as usize);
        cs.desymmetrize();
        cs.filter_by_hfk(&hfk);
        assert_eq!(cs.regions, vec![Region::simple(0, 0)]);
    }
}
