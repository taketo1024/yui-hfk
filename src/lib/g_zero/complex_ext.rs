use yui_core::num::FF2;
use yui_homology::Summand;

use crate::grid_complex::{GCFiltered, GCInfGen, GridComplex};
use super::Region;

impl GridComplex<GCFiltered> {
    // Returns the (alg_deg, alex_deg) coordinates outside r that can host inf_gens at maslov_deg k.
    //
    // For each inflated alex_deg j and alg_deg i ≤ 0, keeps (i, j) only when:
    //   - (i, j) ∉ r, and
    //   - a raw gen exists at the corresponding raw bidegree (k − 2i, j − i).
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

    // Collect all GCInfGens of maslov_deg = k whose (alg_deg, alex_deg) lies outside r.
    //
    // Uses complementary_points to enumerate only the relevant (alg, alex) coordinates,
    // then calls inf_gens_e for the exact U-degree at each point — avoiding a full sweep
    // of inf_gens over all e at every alex degree.
    pub(crate) fn complementary_gens(&self, k: isize, r: &Region) -> Summand<GCInfGen, FF2> {
        let gens = self.complementary_points(k, r).into_iter().flat_map(|(i, j)| {
            self.inf_gens_e(k, j, (-i) as usize)
        });
        Summand::from_raw_gens(gens)
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

    // ---- complementary_gens tests ----

    #[test]
    fn test_complementary_gens_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);

        // The unknot has a single raw gen at (maslov=0, alex=0).
        // Its only inf_gen at (k=0, j=0) is the trivial one: alg=0, alex=0.
        // R(0,0) contains (0,0), so complementary_gens(0, R(0,0)) = empty.
        let gen0 = cx.complementary_gens(0, &Region::simple(0, 0));
        assert_eq!(gen0.dim(), 0);

        // R(-1,-1) does not contain (0,0), so complementary_gens(0, R(-1,-1)) has 1 gen.
        let gen1 = cx.complementary_gens(0, &Region::simple(-1, -1));
        assert_eq!(gen1.dim(), 1);
    }

}
