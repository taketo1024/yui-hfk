use std::marker::PhantomData;
use std::sync::Arc;
use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_homology::{ChainComplex2, Grid, isize2, Summand};
use crate::grid_diagram::{GridDiagram, GridRect};
use super::gen::GridComplexGen;
use super::gen_table::GridComplexGenTable;

pub trait GridComplexVariant {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool;
}

pub struct Tilde;     // [Book] p.72,  Def 4.4.1
pub struct Hat;       // [Book] p.80,  Def 4.6.12.
pub struct Minus;     // [Book] p.75,  Def 4.6.1
pub struct Filtered;  // [Book] p.252, Def 13.2.1

impl GridComplexVariant for Tilde {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        !rect.intersects(&diagram.os, grid_size, false) &&
        !rect.intersects(&diagram.xs, grid_size, false)
    }
}

impl GridComplexVariant for Hat {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        let n = diagram.grid_number() as usize;
        !rect.intersects(&diagram.xs, grid_size, false) &&
        !rect.contains(diagram.os[n - 1], grid_size, false)
    }
}

impl GridComplexVariant for Minus {
    fn is_admissible(rect: &GridRect, diagram: &GridDiagram, grid_size: u8) -> bool {
        !rect.intersects(&diagram.xs, grid_size, false)
    }
}

impl GridComplexVariant for Filtered {
    fn is_admissible(_rect: &GridRect, _diagram: &GridDiagram, _grid_size: u8) -> bool {
        true
    }
}

pub struct GridComplex<V: GridComplexVariant> {
    pub diagram: GridDiagram,
    pub gens: GridComplexGenTable,
    _phantom: PhantomData<V>,
}

impl<V: GridComplexVariant> GridComplex<V> {
    pub fn new(diagram: GridDiagram) -> Self {
        let generators = GridComplexGenTable::from(&diagram);
        GridComplex { diagram, gens: generators, _phantom: PhantomData }
    }

    pub fn grid_number(&self) -> u8 {
        self.diagram.grid_number()
    }

    pub fn into_complex(self) -> ChainComplex2<GridComplexGen, FF2>
    where V: Send + Sync + 'static
    {
        let this = Arc::new(self);

        let support: Vec<isize2> = this.gens.support()
            .map(|(m, a)| isize2(m, a))
            .collect();

        let summands = {
            let this = this.clone();
            Grid::generate(support, move |idx| {
                let gens = this.gens.gens_of_degree(idx.0, idx.1).cloned();
                Summand::from_raw_gens(gens)
            })
        };

        ChainComplex2::new(summands, isize2(-1, 0), move |_, x| {
            this.differentiate(x)
        })
    }

    fn differentiate(&self, x: &Lc<GridComplexGen, FF2>) -> Lc<GridComplexGen, FF2> {
        let mut result = Lc::new();
        for (gen, _) in x.iter() {
            result += self.differentiate_gen(gen);
        }
        result
    }

    fn differentiate_gen(&self, x: &GridComplexGen) -> Lc<GridComplexGen, FF2> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid_diagram::GridDiagram;
    use yui_homology::{ChainComplexTrait, SummandTrait};

    fn unknot_twist() -> GridDiagram { 
        GridDiagram::from_arc_presentation(&[1,2,3,1,2,3])
    }

    #[test]
    fn test_tilde_unknot() {
        let gd = GridDiagram::unknot();
        let cx = GridComplex::<Tilde>::new(gd);
        assert_eq!(cx.grid_number(), 2);
        assert_eq!(cx.gens.len(), 2);

        let cx = cx.into_complex();
        cx.check_d_all();

        let h = cx.homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(-1, -1)].rank(), 1);
    }

    #[test]
    fn test_tilde_unknot_twist() {
        let gd = unknot_twist();
        let cx = GridComplex::<Tilde>::new(gd);
        assert_eq!(cx.grid_number(), 3);
        assert_eq!(cx.gens.len(), 6);

        let cx = cx.into_complex();
        cx.check_d_all();

        // note: not invariant
        let h = cx.homology();
        assert_eq!(h[(0, 0)].rank(), 1);
        assert_eq!(h[(-1, -1)].rank(), 2);
        assert_eq!(h[(-2, -2)].rank(), 1);
    }

    #[test]
    fn test_tilde_trefoil_d() {
        let gd = GridDiagram::trefoil();
        let cx = GridComplex::<Tilde>::new(gd);
        assert_eq!(cx.grid_number(), 5);
        assert_eq!(cx.gens.len(), 120);

        let cx = cx.into_complex();
        cx.check_d_all();

        let h = cx.homology();
        assert_eq!(h[(0, 1)].rank(), 1);
        assert_eq!(h[(-1, 0)].rank(), 5);
        assert_eq!(h[(-2, -1)].rank(), 11);
        assert_eq!(h[(-3, -2)].rank(), 14);
        assert_eq!(h[(-4, -3)].rank(), 11);
        assert_eq!(h[(-5, -4)].rank(), 5);
        assert_eq!(h[(-6, -5)].rank(), 1);
    }
}
