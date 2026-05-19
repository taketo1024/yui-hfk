use std::collections::HashSet;
use std::marker::PhantomData;
use std::ops::RangeInclusive;
use itertools::Itertools;
use num_traits::{One, Zero};
use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_homology::{ChainComplex, ChainComplex2, Grid, Homology, Homology2, Summand, isize2, range_of};
use yui_core::poly::MultiVar;
use crate::grid_complex::gc_gen::GCInfGen;
use crate::grid_diagram::{GridDiagram, GridRect};
use super::gc_gen::GCGen;
use super::gen_table::GCGenSet;

pub trait GCType {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool;
    fn nvars(grid_number: u8) -> usize;
}

pub struct GCTilde;     // [Book] p.72,  Def 4.4.1
pub struct GCHat;       // [Book] p.80,  Def 4.6.12.
pub struct GCMinus;     // [Book] p.75,  Def 4.6.1
pub struct GCFiltered;  // [Book] p.252, Def 13.2.1

impl GCType for GCTilde {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        !rect.intersects(&diagram.os, grid_size, false) &&
        !rect.intersects(&diagram.xs, grid_size, false)
    }
    fn nvars(_: u8) -> usize { 0 }
}

impl GCType for GCHat {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        let n = diagram.grid_number() as usize;
        !rect.intersects(&diagram.xs, grid_size, false) &&
        !rect.contains(diagram.os[n - 1], grid_size, false)
    }
    fn nvars(n: u8) -> usize { (n - 1) as usize }
}

impl GCType for GCMinus {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        !rect.intersects(&diagram.xs, grid_size, false)
    }
    fn nvars(n: u8) -> usize { n as usize }
}

impl GCType for GCFiltered {
    fn is_admissible(_rect: &GridRect, _diagram: &GridDiagram, _grid_size: u8) -> bool {
        true
    }
    fn nvars(n: u8) -> usize { n as usize }
}

pub trait GCBigradedType: GCType {}
impl GCBigradedType for GCTilde {}
impl GCBigradedType for GCHat {}
impl GCBigradedType for GCMinus {}

pub trait GCInfType: GCType {}
impl GCInfType for GCHat {}
impl GCInfType for GCMinus {}
impl GCInfType for GCFiltered {}

pub struct GridComplex<V: GCType> {
    pub diagram: GridDiagram,
    pub gens: GCGenSet,
    _phantom: PhantomData<V>,
}

impl<V: GCType> GridComplex<V> {
    pub fn new(diagram: GridDiagram) -> Self {
        Self::new_filtered(diagram, |_| true)
    }

    pub fn new_filtered(diagram: GridDiagram, pred: impl Fn(&GCGen) -> bool) -> Self {
        let generators = GCGenSet::generate_filtered(&diagram, pred);
        GridComplex { diagram, gens: generators, _phantom: PhantomData }
    }

    pub fn grid_number(&self) -> u8 {
        self.diagram.grid_number()
    }

    pub fn maslov_deg_range(&self) -> RangeInclusive<isize> { 
        self.gens.maslov_deg_range()
    }

    pub fn alex_deg_range(&self) -> RangeInclusive<isize> { 
        self.gens.alex_deg_range()
    }

    fn raw_gens(&self, i: isize, j: isize) -> &[GCGen] {
        self.gens.raw_gens(i, j)
    }

    pub(crate) fn inf_gens(&self, i: isize, j: isize) -> Vec<GCInfGen> {
        let i_max = *self.gens.maslov_deg_range().end();
        let max_e = if i <= i_max { (i_max - i) as usize / 2 } else { 0 };
        (0..=max_e).flat_map(|e| self.inf_gens_e(i, j, e)).collect_vec()
    }

    // Generators at (maslov=i, alex=j) with total U-degree exactly e.
    // Raw gens sit at (i + 2e, j + e); each is tensored with all monomials of degree e.
    //
    //              * x ← (i + 2e, j + e)
    //              |
    //              | e
    //        2e    |
    //  * ----------
    //  ↑
    //  U^e ⊗ x ... (i, j)
    //
    pub(crate) fn inf_gens_e(&self, i: isize, j: isize, e: usize) -> Vec<GCInfGen> {
        let n = V::nvars(self.grid_number());
        self.raw_gens(i + 2 * e as isize, j + e as isize)
            .iter()
            .flat_map(move |x| {
                (0..n).combinations_with_replacement(e).map(move |combo| {
                    let u = MultiVar::from_iter(
                        (0..n).map(|k| (k, combo.iter().filter(|&&v| v == k).count()))
                    );
                    GCInfGen::new(u, x.clone())
                })
            })
            .collect_vec()
    }

    fn differentiate_gen(&self, x: &GCGen) -> Vec<(GCGen, GridRect)> {
        let n = self.grid_number();
        let grid_size = 2 * n;
        let pts = x.points();
        let seq = x.sequence();
        let mut tmp = vec![0u8; n as usize];

        (0..n as usize).tuple_combinations()
            .flat_map(|(i, j)| [(pts[i], pts[j]), (pts[j], pts[i])])
            .filter_map(|(p, q)| {
                let rect = GridRect::between(p, q, grid_size);
                if rect.intersects(&pts, grid_size, true) { return None; }
                if !V::is_admissible(&rect, &self.diagram, grid_size) { return None; }

                let ri = ((rect.origin.x / 2) % n) as usize;
                let rj = ((rect.end().x  / 2) % n) as usize;

                tmp.copy_from_slice(&seq);
                tmp.swap(ri, rj);

                let code = GCGen::encode_in_place(&mut tmp);
                self.gens.gen_from_code(code).cloned().map(|g| (g, rect))
            })
            .collect()
    }
}

impl GridComplex<GCTilde> {
    fn differentiate(&self, z: &Lc<GCGen, FF2>) -> Lc<GCGen, FF2> {
        z.iter().filter_map(|(x, r)| {
            if r.is_zero() { return None }
            let dx = self.differentiate_gen(x).into_iter().map(|(y, _)| 
                (y, FF2::one())
            );
            Some(dx)
        }).flatten().collect()
    }

    pub fn into_complex(self) -> ChainComplex2<GCGen, FF2> {
        let support = self.gens.support().map(|&(i, j)| isize2(i, j)).collect::<HashSet<_>>();

        let this = &self;
        let summands = {
            Grid::generate(support, move |idx| {
                Summand::from_raw_generators(this.raw_gens(idx.0, idx.1).iter().cloned())
            })
        };

        ChainComplex2::new(summands, isize2(-1, 0), move |_, x| {
            self.differentiate(x)
        })
    }

    pub fn into_homology(self) -> Homology2<GCGen, FF2> {
        self.into_complex().reduced().homology()
    }
}

impl<V: GCInfType> GridComplex<V> {
    pub fn differentiate(&self, z: &Lc<GCInfGen, FF2>) -> Lc<GCInfGen, FF2> {
        let n = V::nvars(self.grid_number());
        let grid_size = 2 * self.grid_number();
        let os = &self.diagram.os;

        z.iter().filter_map(|(t, r)| {
            if r.is_zero() { return None }
            let dx = self.differentiate_gen(&t.x).into_iter().map(|(y, rect)| {
                // U_k appears in the differential if Os[k] is inside the rect.
                let v = MultiVar::from_iter(
                    (0..n).filter_map(|k| 
                        rect.contains(os[k], grid_size, false).then_some((k, 1))
                    )
                );
                (GCInfGen::new(&t.u * v, y), FF2::one())
            });
            Some(dx)
        }).flatten().collect()
    }
}

// Implementation for GCHat and GCMinus.
impl<V: GCBigradedType + GCInfType + Send + Sync + 'static> GridComplex<V> {
    pub fn into_complex(self) -> ChainComplex2<GCInfGen, FF2> {
        let a_range = self.gens.alex_deg_range();
        let a_min = *a_range.start();
        let m_range = range_of(self.gens.support().flat_map(|&(i, j)| 
            [i - 2 * (j - a_min), i]
        ));

        let support: Vec<isize2> = a_range.rev().flat_map(|j|
            m_range.clone().rev().map(move |i| isize2(i, j))
        ).collect();

        let this = &self;
        let summands = Grid::generate(support, move |idx| {
            Summand::from_raw_generators(this.inf_gens(idx.0, idx.1))
        });

        ChainComplex2::new(summands, isize2(-1, 0), move |_, z| {
            self.differentiate(z)
        })
    }

    pub fn into_homology(self) -> Homology2<GCInfGen, FF2> {
        self.into_complex().reduced().homology()
    }
}

// Implementation for GCFiltered.
impl GridComplex<GCFiltered> {
    pub fn into_complex(self, maslov_bound: isize) -> ChainComplex<GCInfGen, FF2> {
        let m_max = *self.gens.maslov_deg_range().end();
        let m_range = maslov_bound ..= m_max;
        let a_range = range_of(self.gens.support().flat_map(|&(i, j)| { 
            if i >= maslov_bound { 
                vec![j - (i - maslov_bound) / 2, j]
            } else {
                vec![]
            }
        }));

        let this = &self;
        let summands = Grid::generate(m_range.rev(), move |i| {
            let gens = a_range.clone().flat_map(|j| this.inf_gens(i, j));
            Summand::from_raw_generators(gens)
        });

        ChainComplex::new(summands, -1, move |_, z| {
            self.differentiate(z)
        })
    }

    pub fn into_homology(self, maslov_bound: isize) -> Homology<GCInfGen, FF2> {
        self.into_complex(maslov_bound - 2).reduced().homology().truncated(maslov_bound ..= 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid_diagram::GridDiagram;
    use yui_homology::SummandTrait;

    fn unknot_twist() -> GridDiagram { 
        GridDiagram::from_arc_presentation(&[1,2,3,1,2,3])
    }

    #[test]
    fn test_tilde_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCTilde>::new(gd);
        assert_eq!(cx.grid_number(), 2);
        assert_eq!(cx.gens.len(), 2);

        let h = cx.into_homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(-1, -1)].rank(), 1);
    }

    #[test]
    fn test_tilde_unknot_twist() {
        let gd = unknot_twist();
        let cx = GridComplex::<GCTilde>::new(gd);
        assert_eq!(cx.grid_number(), 3);
        assert_eq!(cx.gens.len(), 6);

        // note: not invariant
        let h = cx.into_homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(-1, -1)].rank(), 2);
        assert_eq!(h[(-2, -2)].rank(), 1);
    }

    #[test]
    fn test_tilde_trefoil() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCTilde>::new(gd);
        assert_eq!(cx.grid_number(), 5);
        assert_eq!(cx.gens.len(), 120);
        assert_eq!(cx.maslov_deg_range(), -6..=0);
        assert_eq!(cx.alex_deg_range(), -5..=1);

        let h = cx.into_homology();
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(-1, 0)].rank(), 5);
        assert_eq!(h[(-2, -1)].rank(), 11);
        assert_eq!(h[(-3, -2)].rank(), 14);
        assert_eq!(h[(-4, -3)].rank(), 11);
        assert_eq!(h[(-5, -4)].rank(), 5);
        assert_eq!(h[(-6, -5)].rank(), 1);
    }

    #[test]
    fn test_tilde_trefoil_mirror() {
        let gd = GridDiagram::trefoil().rotate90();
        let cx = GridComplex::<GCTilde>::new(gd);

        let h = cx.into_homology();
        assert_eq!(h[(2, 1)].rank(), 1);
        assert_eq!(h[(1, 0)].rank(), 5);
        assert_eq!(h[(0, -1)].rank(), 11);
        assert_eq!(h[(-1, -2)].rank(), 14);
        assert_eq!(h[(-2, -3)].rank(), 11);
        assert_eq!(h[(-3, -4)].rank(), 5);
        assert_eq!(h[(-4, -5)].rank(), 1);
    }

    #[test]
    fn test_tilde_trefoil_filtered() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCTilde>::new_filtered(gd, |x| x.alex_deg >= 0);

        assert_eq!(cx.grid_number(), 5);
        assert_eq!(cx.gens.len(), 6);
        assert_eq!(cx.maslov_deg_range(), -1..=0);
        assert_eq!(cx.alex_deg_range(), 0..=1);

        let h = cx.into_homology();
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(-1, 0)].rank(), 5);
        assert_eq!(h[(-2, -1)].rank(), 0);
        assert_eq!(h[(-3, -2)].rank(), 0);
        assert_eq!(h[(-4, -3)].rank(), 0);
        assert_eq!(h[(-5, -4)].rank(), 0);
        assert_eq!(h[(-6, -5)].rank(), 0);
    }

    #[test]
    fn test_inf_gens_tilde() {
        // GCTilde has nvars = 0, so inf_gens(i, j) wraps raw_gens(i, j) with the unit monomial.
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCTilde>::new(gd);

        let unit = MultiVar::<'U', usize>::default();

        let gens_00 = cx.inf_gens(0, 0);
        assert_eq!(gens_00.len(), 1);
        assert_eq!(gens_00[0].u, unit);

        let gens_m1 = cx.inf_gens(-1, -1);
        assert_eq!(gens_m1.len(), 1);
        assert_eq!(gens_m1[0].u, unit);

        assert_eq!(cx.inf_gens(0, 1).len(), 0);
    }

    #[test]
    fn test_inf_gens_minus() {
        // GCMinus on unknot (n=2): nvars = 2, raw gens at (0,0) and (-1,-1).
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCMinus>::new(gd);

        // e=0 only: 1 raw gen × 1 monomial
        assert_eq!(cx.inf_gens(0, 0).len(), 1);

        // e=0: nothing at (-2,-1); e=1: 1 raw gen at (0,0) × 2 monomials (U_0, U_1)
        let v = cx.inf_gens(-2, -1);
        assert_eq!(v.len(), 2);
        assert!(v.iter().all(|t| t.maslov_deg() == -2 && t.alex_deg() == -1));

        // e=2: 1 raw gen at (0,0) × 3 monomials (U_0^2, U_0*U_1, U_1^2)
        let v = cx.inf_gens(-4, -2);
        assert_eq!(v.len(), 3);
        assert!(v.iter().all(|t| t.maslov_deg() == -4 && t.alex_deg() == -2));
    }

    #[test]
    fn test_differentiate_hat_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCHat>::new(gd);

        for (i, j) in [(0isize, 0isize), (-1, -1), (-2, -1)] {
            for t in cx.inf_gens(i, j) {
                let dz = cx.differentiate(&Lc::from(t));
                for (s, _) in dz.iter() {
                    assert_eq!(s.maslov_deg(), i - 1);
                    assert_eq!(s.alex_deg(), j);
                }
                let ddz = cx.differentiate(&dz);
                assert!(ddz.is_zero());
            }
        }
    }

    #[test]
    fn test_differentiate_minus_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<GCMinus>::new(gd);

        // Check degree preservation and d² = 0 at several bidegrees.
        // Includes generators with non-trivial U-monomials (e.g. at (-2,-1), (-3,-2)).
        for (i, j) in [(0isize, 0isize), (-1, -1), (-2, -1), (-3, -2)] {
            for t in cx.inf_gens(i, j) {
                let dz = cx.differentiate(&Lc::from(t));
                for (s, _) in dz.iter() {
                    assert_eq!(s.maslov_deg(), i - 1,);
                    assert_eq!(s.alex_deg(), j);
                }
                let ddz = cx.differentiate(&dz);
                assert!(ddz.is_zero());
            }
        }
    }

    // -- GCHat tests --

    #[test]
    fn test_hat_unknot() {
        let gd = GridDiagram::unknot();
        let h = GridComplex::<GCHat>::new(gd).into_homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h.total_rank(), 1);
    }

    #[test]
    fn test_hat_unknot_twist() {
        let gd = GridDiagram::from_arc_presentation(&[1,2,3,1,2,3]);
        let h = GridComplex::<GCHat>::new(gd).into_homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h.total_rank(), 1);
    }

    #[test]
    fn test_hat_trefoil() {
        let gd = GridDiagram::trefoil();
        let h = GridComplex::<GCHat>::new(gd).into_homology();
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(-1, 0)].rank(), 1);
        assert_eq!(h[(-2, -1)].rank(), 1);
        assert_eq!(h.total_rank(), 3);
    }

    #[test]
    fn test_hat_trefoil_mirror() {
        let gd = GridDiagram::trefoil().rotate90();
        let h = GridComplex::<GCHat>::new(gd).into_homology();
        assert_eq!(h[(2, 1)].rank(), 1);
        assert_eq!(h[(1, 0)].rank(), 1);
        assert_eq!(h[(0, -1)].rank(), 1);
        assert_eq!(h.total_rank(), 3);
    }

    // -- GCMinus tests --

    #[test]
    fn test_minus_unknot() {
        let gd = GridDiagram::unknot();
        let h = GridComplex::<GCMinus>::new_filtered(gd, |x| x.alex_deg >= -2).into_homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(-2, -1)].rank(), 1);
        assert_eq!(h.total_rank(), 2);
    }

    #[test]
    fn test_minus_trefoil() {
        let gd = GridDiagram::trefoil();
        let h = GridComplex::<GCMinus>::new_filtered(gd, |x| x.alex_deg >= -2).into_homology();
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(-2, -1)].rank(), 1);
        assert_eq!(h[(-4, -2)].rank(), 1);
        assert_eq!(h.total_rank(), 3);
    }

    // -- GCFiltered tests --

    #[test]
    fn test_filt_unknot() {
        let gd = GridDiagram::from_arc_presentation(&[1,2,3,1,2,3]);
        let gc = GridComplex::<GCFiltered>::new(gd);
        let gh = gc.into_homology(-4);
        assert_eq!(gh[0].rank(), 1);
        assert_eq!(gh[-2].rank(), 1);
        assert_eq!(gh[-4].rank(), 1);
        assert_eq!(gh.total_rank(), 3);
    }

    #[test]
    fn test_filt_trefoil() {
        let gd = GridDiagram::trefoil();
        let gc = GridComplex::<GCFiltered>::new(gd);
        let gh = gc.into_homology(-4);
        assert_eq!(gh[0].rank(), 1);
        assert_eq!(gh[-2].rank(), 1);
        assert_eq!(gh[-4].rank(), 1);
        assert_eq!(gh.total_rank(), 3);
    }
}
