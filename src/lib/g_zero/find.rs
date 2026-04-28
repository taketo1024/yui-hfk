use log::debug;
use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_homology::{Summand, SummandTrait};
use yui_matrix::MatTrait;
use yui_matrix::sparse::SpVec;
use yui_matrix::sparse::pluq::solve_pluq;

use crate::grid_complex::{GCFiltered, GCInfGen, GridComplex};
use crate::grid_complex::chain::GCInfChainExt;
use super::Region;

impl GridComplex<GCFiltered> {
    /// Tests whether `z` has a homologous cycle contained in the region `r`.
    ///
    /// Finds `w` of maslov degree `k+1` (where `k = maslov_degree(z)`) with
    /// generators outside `r` such that `d(w) ≡ z (mod r)`.
    /// Returns `z - d(w)`, which is homologous to `z` and lies inside `r`.
    pub fn find_homologous(
        &self,
        z: &Lc<GCInfGen, FF2>,
        r: &Region,
    ) -> Option<Lc<GCInfGen, FF2>> {
        assert!(z.is_maslov_homogeneous());
        let k = z.maslov_degree();

        debug!("inflate generators ...");
        let gen_w = self.complementary_gens(k + 1, r); // space where w lives
        debug!("  in-gens: {}", gen_w.rank());

        let gen_z = self.complementary_gens(k, r);     // space where z lives mod r
        debug!("  out-gens: {}", gen_z.rank());

        let z0 = cutoff(z, r);
        let y  = gen_z.vectorize(&z0);

        debug!("create matrix ...");
        let a  = gen_w.make_matrix(&gen_z, |x| cutoff(&self.differentiate(x), r));
        debug!("  a: {} x {}", a.nrows(), a.ncols());

        debug!("solve linear system ...");
        let x_sol = solve_pluq(&a, &y.to_dense())?;
        let x_vec = SpVec::from(x_sol);
        let w = gen_w.devectorize(&x_vec);

        // z - d(w) is homologous to z and lies inside r.
        let mut z1 = z.clone();
        z1 -= &self.differentiate(&w);
        Some(z1)
    }

    // Collect all GCInfGens of maslov_deg = k whose (alg_deg, alex_deg) lies outside r.
    //
    // Uses complementary_points to enumerate only the relevant (alg, alex) coordinates,
    // then calls inf_gens_e for the exact U-degree at each point — avoiding a full sweep
    // of inf_gens over all e at every alex degree.
    fn complementary_gens(&self, k: isize, r: &Region) -> Summand<GCInfGen, FF2> {
        let gens = self.complementary_points(k, r).into_iter().flat_map(|(i, j)| {
            self.inf_gens_e(k, j, (-i) as usize)
        });
        Summand::from_raw_gens(gens)
    }

    // Returns the (alg_deg, alex_deg) coordinates outside r that can host inf_gens at maslov_deg k.
    //
    // For each inflated alex_deg j and alg_deg i ≤ 0, keeps (i, j) only when:
    //   - (i, j) ∉ r, and
    //   - a raw gen exists at the corresponding raw bidegree (k − 2i, j − i).
    fn complementary_points(&self, k: isize, r: &Region) -> Vec<(isize, isize)> {
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
}

// Drop terms of `z` whose (alg_deg, alex_deg) is inside r.
fn cutoff(z: &Lc<GCInfGen, FF2>, r: &Region) -> Lc<GCInfGen, FF2> {
    z.filter_gens(|g| !r.contains_gen(g))
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_homology::SummandTrait;
    use super::*;
    use crate::grid_complex::{GCInfGen};
    use crate::{GridComplex, GridDiagram};
    use crate::grid_complex::chain::tests::{make_gen, lc_of};
    use super::super::{Corner, Region};

    // ---- cutoff tests ----

    #[test]
    fn test_cutoff_empty() {
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert!(cutoff(&z, &Region::simple(0, 0)).is_zero());
    }

    #[test]
    fn test_cutoff_drops_inside_keeps_outside() {
        // g_in:  (alg=0, alex=0)  → inside  R(0,0)  → dropped
        // g_out: (alg=0, alex=1)  → outside R(0,0) (alex 1 > 0) → kept
        let g_in  = make_gen(0, 2, 0, 0, &[]);  // alg=0, alex=0
        let g_out = make_gen(1, 2, 0, 1, &[]);  // alg=0, alex=1
        let lc = lc_of(&[g_in, g_out.clone()]);
        let r = Region::simple(0, 0);
        let result = cutoff(&lc, &r);
        // Only g_out should remain
        assert_eq!(result.nterms(), 1);
        assert!(!result.coeff(&g_out).is_zero());
    }

    #[test]
    fn test_cutoff_all_outside() {
        // U^1 generator: alg = -1, alex = alex_raw - 1
        // With alex_raw=2: alex = 1, coord = (-1, 1) → outside R(0,0)? -1 ≤ 0 and 1 > 0 → outside
        let g = make_gen(0, 2, 0, 2, &[1]); // alg=-1, alex=2-1=1
        let lc = lc_of(&[g.clone()]);
        let r = Region::simple(0, 0);
        let result = cutoff(&lc, &r);
        assert_eq!(result.nterms(), 1);
        assert!(!result.coeff(&g).is_zero());
    }

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

    // ---- find_homologous tests ----

    #[test]
    fn test_find_homologous_unknot_already_inside() {
        // For the unknot the hom_gen lies at (alg=0, alex=0) ∈ R(0,0),
        // so z0 = cutoff(z, R(0,0)) = 0.  The system is trivially satisfied
        // with w = 0, and the function returns z unchanged.
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        let r = Region::simple(0, 0);

        let z1 = cx.find_homologous(&z, &r).expect("expected Some");
        assert!(cx.differentiate(&z1).is_zero());
        assert!(r.contains_chain(&z1));
    }

    #[test]
    fn test_find_homologous_unknot_none() {
        // R(-1,-1) lies strictly below the hom_gen at (0,0).
        // The unknot has no generators of maslov degree 1, so no w can be found.
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        let r = Region::new([Corner::new(-1, 0), Corner::new(0, -1)]);

        assert!(cx.find_homologous(&z, &r).is_none());
    }

    #[test]
    fn test_find_homologous_trefoil() {
        // For the trefoil (genus 1), R(0,1) is the expected realizable region.
        // find_homologous should return Some(z1) where z1 is a cycle inside R(0,1).
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCFiltered>::new(gd);
        let z = cx.hom_gen();
        let r = Region::simple(0, 1);

        let z1 = cx.find_homologous(&z, &r).expect("expected Some for trefoil R(0,1)");
        assert!(cx.differentiate(&z1).is_zero(), "z1 must be a cycle");
        assert!(r.contains_chain(&z1), "z1 must lie inside r");
    }
}
