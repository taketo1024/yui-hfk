use itertools::Itertools;
use num_traits::Zero;
use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_core::poly::MultiVar;
use yui_homology::{Grid2, isize2};

use crate::grid_complex::GCInfGen;

pub trait GCInfChainExt {
    fn is_maslov_homogeneous(&self) -> bool;
    fn maslov_degree(&self) -> isize;
    fn alex_degree(&self) -> isize;
    fn alg_degree(&self) -> isize;
    fn mul(&self, u: &MultiVar<'U', usize>) -> Self;
    fn support_grid(&self) -> Grid2<usize>;
}

impl GCInfChainExt for Lc<GCInfGen, FF2> {
    fn is_maslov_homogeneous(&self) -> bool {
        let mut degs = self.iter()
            .filter(|(_, c)| !c.is_zero())
            .map(|(g, _)| g.maslov_deg());
        match degs.next() {
            None => true,
            Some(first) => degs.all(|d| d == first),
        }
    }

    fn maslov_degree(&self) -> isize { max_deg(self, GCInfGen::maslov_deg) }
    fn alex_degree(&self)   -> isize { max_deg(self, GCInfGen::alex_deg)   }
    fn alg_degree(&self)    -> isize { max_deg(self, GCInfGen::alg_deg)    }

    fn mul(&self, u: &MultiVar<'U', usize>) -> Self {
        self.map_ref(|g, r| (GCInfGen::new(&g.u * u, g.x.clone()), *r))
    }

    fn support_grid(&self) -> Grid2<usize> {
        self.iter()
            .filter(|(_, c)| !c.is_zero())
            .counts_by(|(g, _)| isize2(g.alg_deg(), g.alex_deg()))
            .into_iter()
            .collect()
    }
}

fn max_deg(lc: &Lc<GCInfGen, FF2>, f: impl Fn(&GCInfGen) -> isize) -> isize {
    lc.iter()
        .filter(|(_, c)| !c.is_zero())
        .map(|(g, _)| f(g))
        .max()
        .unwrap_or(0)
}

#[cfg(test)]
pub(crate) mod tests {
    use num_traits::One;
    use yui_core::poly::MultiVar;
    use super::*;
    use crate::grid_complex::{GCGen, GCInfGen};

    // code distinguishes generators; grid_number must be >= u_degs.len()
    pub fn make_gen(code: usize, grid_number: u8, maslov: isize, alex: isize, u_degs: &[usize]) -> GCInfGen {
        let base = GCGen::from_raw(code, grid_number, maslov, alex);
        let u = MultiVar::<'U', usize>::from_iter(
            u_degs.iter().copied().enumerate().filter(|(_, d)| *d > 0)
        );
        GCInfGen::new(u, base)
    }

    pub fn lc_of(gens: &[GCInfGen]) -> Lc<GCInfGen, FF2> {
        gens.iter().cloned().map(|g| (g, FF2::one())).collect()
    }

    #[test]
    fn test_is_maslov_homogeneous_empty() {
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert!(z.is_maslov_homogeneous());
    }

    #[test]
    fn test_is_maslov_homogeneous_single() {
        let lc = lc_of(&[make_gen(0, 2, 3, 1, &[])]);
        assert!(lc.is_maslov_homogeneous());
    }

    #[test]
    fn test_is_maslov_homogeneous_same() {
        // both have maslov_deg() = 2 (no U-factors)
        let g1 = make_gen(0, 2, 2, 0, &[]);
        let g2 = make_gen(1, 2, 2, 1, &[]);
        let lc = lc_of(&[g1, g2]);
        assert!(lc.is_maslov_homogeneous());
    }

    #[test]
    fn test_is_maslov_homogeneous_different() {
        // distinct codes → don't cancel; maslov_deg 2 vs 4
        let g1 = make_gen(0, 2, 2, 0, &[]);
        let g2 = make_gen(1, 2, 4, 0, &[]);
        let lc = lc_of(&[g1, g2]);
        assert!(!lc.is_maslov_homogeneous());
    }

    #[test]
    fn test_maslov_degree_empty() {
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert_eq!(z.maslov_degree(), 0);
    }

    #[test]
    fn test_maslov_degree() {
        // maslov_deg() = gen.maslov_deg - 2 * u.total_deg
        let g1 = make_gen(0, 2, 4, 0, &[]);    // 4 - 0 = 4
        let g2 = make_gen(1, 2, 6, 0, &[1]);   // 6 - 2 = 4
        let g3 = make_gen(2, 2, 2, 0, &[]);    // 2 - 0 = 2
        let lc = lc_of(&[g1, g2, g3]);
        assert_eq!(lc.maslov_degree(), 4);
    }

    #[test]
    fn test_alex_degree_empty() {
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert_eq!(z.alex_degree(), 0);
    }

    #[test]
    fn test_alex_degree() {
        // alex_deg() = gen.alex_deg - u.total_deg
        let g1 = make_gen(0, 2, 0, 3, &[]);    // 3 - 0 = 3
        let g2 = make_gen(1, 2, 0, 5, &[2]);   // 5 - 2 = 3
        let g3 = make_gen(2, 2, 0, 2, &[]);    // 2 - 0 = 2
        let lc = lc_of(&[g1, g2, g3]);
        assert_eq!(lc.alex_degree(), 3);
    }

    #[test]
    fn test_alg_degree_empty() {
        let z: Lc<GCInfGen, FF2> = Lc::new();
        assert_eq!(z.alg_degree(), 0);
    }

    #[test]
    fn test_alg_degree() {
        // alg_degree of a term = -(u.total_deg); max = least-negative
        let g1 = make_gen(0, 2, 0, 0, &[1]);    // -(1) = -1
        let g2 = make_gen(1, 2, 0, 0, &[0, 3]); // -(3) = -3
        let g3 = make_gen(2, 2, 0, 0, &[]);     // -(0) = 0
        let lc = lc_of(&[g1, g2, g3]);
        assert_eq!(lc.alg_degree(), 0);
    }

    #[test]
    fn test_alg_degree_no_unit() {
        let g1 = make_gen(0, 2, 0, 0, &[2]);    // -(2) = -2
        let g2 = make_gen(1, 2, 0, 0, &[0, 1]); // -(1) = -1
        let lc = lc_of(&[g1, g2]);
        assert_eq!(lc.alg_degree(), -1);
    }

    #[test]
    fn test_support_grid() {
        let g1 = make_gen(0, 2, 0, 0, &[]);     // alg=0,  alex=0
        let g2 = make_gen(1, 2, 0, 3, &[1]);    // alg=-1, alex=3-1=2
        let g3 = make_gen(2, 2, 0, 3, &[0, 1]); // alg=-1, alex=3-1=2 (same bidegree as g2)
        let lc = lc_of(&[g1, g2, g3]);
        let g = lc.support_grid();
        assert_eq!(g[(0, 0)], 1);
        assert_eq!(g[(-1, 2)], 2);
        assert_eq!(g[(0, 1)], 0);
    }
}
