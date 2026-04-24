use std::collections::HashMap;
use std::marker::PhantomData;
use std::ops::RangeInclusive;
use num_integer::binomial;
use num_traits::Zero;
use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_homology::{ChainComplex2, Grid, Grid2, Homology2, Summand, SummandTrait, isize2, range_of};
use crate::grid_diagram::{GridDiagram, GridRect};
use super::gen::GCGen;
use super::gen_table::GCGenSet;

pub trait GCType {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool;
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
}

impl GCType for GCHat {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        let n = diagram.grid_number() as usize;
        !rect.intersects(&diagram.xs, grid_size, false) &&
        !rect.contains(diagram.os[n - 1], grid_size, false)
    }
}

impl GCType for GCMinus {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        !rect.intersects(&diagram.xs, grid_size, false)
    }
}

impl GCType for GCFiltered {
    fn is_admissible(_rect: &GridRect, _diagram: &GridDiagram, _grid_size: u8) -> bool {
        true
    }
}

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
        range_of(self.gens.support().map(|(i, _)| i))
    }

    pub fn alex_deg_range(&self) -> RangeInclusive<isize> { 
        range_of(self.gens.support().map(|(_, j)| j))
    }

    fn differentiate_gen(&self, x: &GCGen) -> Lc<GCGen, FF2> {
        use itertools::Itertools;

        let n = self.grid_number();
        let grid_size = 2 * n;
        let pts = x.points();

        let mut seq = x.sequence();
        let mut result = Lc::new();

        for (i, j) in (0..n as usize).tuple_combinations() {
            for &(p, q) in &[(pts[i], pts[j]), (pts[j], pts[i])] {
                let rect = GridRect::between(p, q, grid_size);

                if rect.intersects(&pts, grid_size, true) {
                    continue;
                }
                if !V::is_admissible(&rect, &self.diagram, grid_size) {
                    continue;
                }

                let ri = (rect.origin.x / 2) % n;
                let rj = (rect.end().x / 2) % n;

                seq.swap(ri as usize, rj as usize);
                let y = self.gens.gen_from_seq(&seq);
                seq.swap(ri as usize, rj as usize);

                if let Some(y) = y {
                    result += Lc::from(y.clone());
                }
            }
        }

        result
    }
}

impl GridComplex<GCTilde> {
    fn differentiate(&self, x: &Lc<GCGen, FF2>) -> Lc<GCGen, FF2> {
        let mut result = Lc::new();
        for (gen, _) in x.iter() {
            result += self.differentiate_gen(gen);
        }
        result
    }

    pub fn into_complex(self) -> ChainComplex2<GCGen, FF2> {
        let this = &self;

        let support: Vec<isize2> = this.gens.support()
            .map(|(m, a)| isize2(m, a))
            .collect();

        let summands = {
            Grid::generate(support, move |idx| {
                let gens = this.gens.gens_of_degree(idx.0, idx.1).cloned();
                Summand::from_raw_gens(gens)
            })
        };

        ChainComplex2::new(summands, isize2(-1, 0), move |_, x| {
            self.differentiate(x)
        })
    }

    pub fn into_homology(self) -> Homology2<GCGen, FF2> {
        self.into_complex().reduced().homology()
    }

    pub fn genus(gd: &GridDiagram) -> usize { 
        let cx = GridComplex::<GCTilde>::new_filtered(gd.clone(), |x| x.alexander_deg > 0);
        let m_range = cx.maslov_deg_range();
        let a_range = cx.alex_deg_range();

        let h = cx.into_homology();
        for j in a_range.rev() { 
            for i in m_range.clone().rev() { 
                if h[(i, j)].rank() > 0 { 
                    return j as usize
                }
            }
        }
        0
    }

    /// Computed from tilde-GC. 
    /// H(tilde-GC(K)) = HFK(K) ⊗ W^(n - 1)
    pub fn hfk(gd: &GridDiagram) -> Grid2<usize> { 
        let n = gd.grid_number();

        let cx = GridComplex::<GCTilde>::new_filtered(gd.clone(), |x| x.alexander_deg >= 0);
        let m_range = cx.maslov_deg_range();
        let a_range = cx.alex_deg_range();
        let h = cx.into_homology();

        let mut ranks: HashMap<(isize, isize), usize> = HashMap::new();

        for j in a_range.clone() { 
            for i in m_range.clone() { 
                if !h[(i, j)].is_zero() { 
                    ranks.insert((i, j), h[(i, j)].rank());
                }
            }
        }

        // divide by W^(n - 1).
        for j in a_range.clone().rev() { 
            for i in m_range.clone() { 
                let Some(&r) = ranks.get(&(i, j)) else { continue };
                if r.is_zero() { continue }

                for k in 1 .. n { 
                    let k = k as isize;
                    if k > j { break }
                    
                    let t = ranks.get_mut(&(i - k, j - k)).unwrap();
                    *t -= r * binomial(n - 1, k as u8) as usize;
                }
            }
        }

        ranks.retain(|_, r| *r > 0);

        // use symmetry: HFK_{i, j}(K) = HFK_{i - 2j, -j}(K)
        for j in a_range.clone().rev() { 
            if j == 0 { break }

            for i in m_range.clone() { 
                let Some(r) = ranks.get(&(i, j)) else { continue };
                ranks.insert((i - 2 * j, -j), *r);
            }
        }

        Grid2::generate(
            ranks.keys().map(|&(i, j)| isize2(i, j)),
            |isize2(i, j)| ranks[&(i, j)]
        )
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
    fn test_tilde_trefoil_filtered() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<GCTilde>::new_filtered(gd, |x| x.alexander_deg >= 0);

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
    fn test_genus_unknot() {
        let gd = GridDiagram::unknot();
        let g = GridComplex::<GCTilde>::genus(&gd);
        assert_eq!(g, 0);
    }

    #[test]
    fn test_genus_trefoil() {
        let gd = GridDiagram::trefoil();
        let g = GridComplex::<GCTilde>::genus(&gd);
        assert_eq!(g, 1);
    }

    #[test]
    fn test_genus_5_1() {
        let gd = GridDiagram::from_sequences(&[5,4,3,2,1,0,6], &[0,6,5,4,3,2,1]);
        let g = GridComplex::<GCTilde>::genus(&gd);
        assert_eq!(g, 2);
    }

    #[test]
    fn test_hfk_unknot() {
        let gd = GridDiagram::unknot();
        let h = GridComplex::<GCTilde>::hfk(&gd);
        assert_eq!(h[(0, 0)], 1);
    }

    #[test]
    fn test_hfk_trefoil() {
        let gd = GridDiagram::trefoil();
        let h = GridComplex::<GCTilde>::hfk(&gd);
        assert_eq!(h[(0, 1)], 1);
        assert_eq!(h[(-1, 0)], 1);
        assert_eq!(h[(-2, -1)], 1);
    }

    #[test]
    fn test_hfk_4_1() {
        let gd = GridDiagram::figure8();
        let h = GridComplex::<GCTilde>::hfk(&gd);
        assert_eq!(h[(1, 1)], 1);
        assert_eq!(h[(0, 0)], 3);
        assert_eq!(h[(-1, -1)], 1);
    }
}
