use num_integer::binomial;
use yui_core::num::FF2;
use yui_homology::Summand;

use crate::grid_complex::{GCFiltered, GCGen, GCInfGen, GridComplex};
use super::Region;

impl GridComplex<GCFiltered> {
    // Returns the (alg_deg, alex_deg) coordinates outside r that can host inf_gens at maslov_deg k.
    pub(crate) fn complementary_points(&self, k: isize, r: &Region) -> Vec<(isize, isize)> {
        let a_range = self.gens.alex_deg_range();
        let j_max = *a_range.end();

        a_range.flat_map(|j| {
            let i_min = -(j_max - j); // <= 0
            (i_min..=0).filter_map(move |i| {
                (!r.contains((i, j)) && self.gens.has_gen_at(k - 2 * i, j - i))
                    .then_some((i, j))
            })
        }).collect()
    }

    // Count GCInfGens of maslov_deg = k at the given (alg, alex) points without allocating them. 
    pub(crate) fn count_inf_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> usize {
        let n = self.grid_number() as usize;
        pts.iter().map(|&(i, j)| {
            let e = (-i) as usize;
            let n_bases = self.gens.raw_gens(k - 2 * i, j - i).len();
            let n_multi = binomial(n + e - 1, e);
            n_bases * n_multi
        }).sum()
    }

    // Collect all GCInfGens of maslov_deg = k at the given (alg, alex) points.
    pub(crate) fn inf_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> Summand<GCInfGen, FF2> {
        let gens = pts.iter().flat_map(|&(i, j)| {
            self.inf_gens_e(k, j, (-i) as usize)
        });
        Summand::from_raw_generators(gens)
    }

    // Count base GCGens at the given points without allocating them.
    pub(crate) fn count_raw_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> usize {
        pts.iter().map(|&(i, j)| {
            self.gens.raw_gens(k - 2 * i, j - i).len()
        }).sum()
    }

    // Collect all base GCGens whose (maslov, alex) corresponds to one of the given
    // inflated (alg, alex) points at maslov_deg = k. 
    pub(crate) fn raw_gens_in(&self, k: isize, pts: &[(isize, isize)]) -> Summand<GCGen, FF2> {
        let gens = pts.iter().flat_map(|&(i, j)| {
            self.gens.raw_gens(k - 2 * i, j - i).iter().cloned()
        });
        Summand::from_raw_generators(gens)
    }
}

#[cfg(test)]
mod tests {
    use yui_homology::SummandTrait;
    use super::*;
    use crate::{GridComplex, GridDiagram};
    use super::super::Region;

    // ---- complementary_points tests ----

    #[test]
    fn test_complementary_points_unknot() {
        // unknot: raw gens at (maslov=0, alex=0) and (maslov=-1, alex=-1), so j_max=0.
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);

        // k=0: only raw gen at (0,0). R(-1,-1) does not contain (0,0) → point (0,0).
        let pts = cx.complementary_points(0, &Region::simple(-1, -1));
        assert_eq!(pts, vec![(0, 0)]);

        // k=-1: only raw gen at (-1,-1). R(-1,-1) does not contain (0,-1) → point (0,-1).
        let pts = cx.complementary_points(-1, &Region::simple(-1, -1));
        assert_eq!(pts, vec![(0, -1)]);

        // R(0,0): contains all (i,j) with i ≤ 0 and j ≤ 0 → every point is inside r.
        let pts = cx.complementary_points(0, &Region::simple(0, 0));
        assert!(pts.is_empty());
        let pts = cx.complementary_points(-1, &Region::simple(0, 0));
        assert!(pts.is_empty());
    }

    // ---- inf_gens_in tests ----

    #[test]
    fn test_inf_gens_in_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);

        // The unknot has a single raw gen at (maslov=0, alex=0).
        // Its only inf_gen at (k=0, j=0) is the trivial one: alg=0, alex=0.
        // R(0,0) contains (0,0), so the points are empty.
        let pts = cx.complementary_points(0, &Region::simple(0, 0));
        let gen0 = cx.inf_gens_in(0, &pts);
        assert_eq!(gen0.dim(), 0);

        // R(-1,-1) does not contain (0,0), so we get 1 gen.
        let pts = cx.complementary_points(0, &Region::simple(-1, -1));
        let gen1 = cx.inf_gens_in(0, &pts);
        assert_eq!(gen1.dim(), 1);
    }

}
