use yui_core::lc::Lc;
use yui_core::num::FF2;
use yui_core::poly::MultiVar;

use crate::{GridComplex, GridDiagram, HFK};
use crate::grid_complex::{GCFiltered, GCInfGen};
use crate::grid_complex::chain::GCInfChainExt;
use super::{Candidates, Region};

pub struct G0Calc {
    gc: GridComplex<GCFiltered>,
    z: Lc<GCInfGen, FF2>,
    candidates: Candidates,
}

impl G0Calc {
    pub fn new(gd: GridDiagram) -> Self {
        let hfk = HFK::from(&gd);
        let gc = GridComplex::<GCFiltered>::new(gd);
        let z = gc.hom_gen();

        let g = hfk.genus() as usize;
        let mut candidates = Candidates::generate_from(g);
        candidates.deduplicate();
        candidates.filter_by_hfk(&hfk);
        candidates.sort();

        G0Calc { gc, z, candidates }
    }

    /// Processes all remaining candidates and populates `candidates.accepted`.
    pub fn run(&mut self) {
        while !self.candidates.regions.is_empty() {
            let r = self.candidates.regions.remove(0);
            self.process(r);
        }
    }

    pub fn accepted(&self) -> &[Region] {
        &self.candidates.accepted
    }

    // Tests whether `r` is realizable: tries to find a cycle homologous to z
    // inside r. Shifts both r and z by s = right_end(r) so that the shifted
    // region has right_end = 0, which is required by find_homologous.
    fn process(&mut self, r: Region) {
        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);
        match self.gc.find_homologous(&z_shifted, &r_shifted) {
            Some(_) => self.candidates.accept(&r),
            None    => self.candidates.reject(&r),
        }
    }

    // Returns z multiplied by U_0^s, shifting all generator coordinates by
    // (-s, -s) to match r.shift(-s, -s).
    fn shift_z(&self, s: usize) -> Lc<GCInfGen, FF2> {
        let u = MultiVar::from((0, s));
        self.z.mul(&u)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GridDiagram;
    use crate::g_zero::Corner;

    fn accepted(gd: GridDiagram) -> Vec<Region> {
        let mut calc = G0Calc::new(gd);
        calc.run();
        let mut rs = calc.accepted().to_vec();
        rs.sort_by_key(|r| r.right_end());
        rs
    }

    #[test]
    fn test_unknot() {
        let rs = accepted(GridDiagram::load("0_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_trefoil() {
        let rs = accepted(GridDiagram::load("3_1").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_trefoil_mirror() {
        let rs = accepted(GridDiagram::load("3_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }

    #[test]
    fn test_figure8() {
        let rs = accepted(GridDiagram::load("4_1").unwrap());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_figure8_mirror() {
        let rs = accepted(GridDiagram::load("4_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::simple(0, 0)]);
    }

    #[test]
    fn test_5_1() {
        let rs = accepted(GridDiagram::load("5_1").unwrap());
        assert_eq!(rs.len(), 3);
        assert!(rs.contains(&Region::simple(0, 2)));
        assert!(rs.contains(&Region::simple(1, 1)));
        assert!(rs.contains(&Region::simple(2, 0)));
    }

    #[test]
    fn test_5_1_mirror() {
        let rs = accepted(GridDiagram::load("5_1").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([
            Corner::new(-2, 0),
            Corner::new(-1, -1),
            Corner::new(0, -2),
        ])]);
    }

    #[test]
    fn test_5_2() {
        let rs = accepted(GridDiagram::load("5_2").unwrap());
        assert_eq!(rs.len(), 2);
        assert!(rs.contains(&Region::simple(0, 1)));
        assert!(rs.contains(&Region::simple(1, 0)));
    }

    #[test]
    fn test_5_2_mirror() {
        let rs = accepted(GridDiagram::load("5_2").unwrap().rotate90());
        assert_eq!(rs, vec![Region::new([Corner::new(-1, 0), Corner::new(0, -1)])]);
    }
}
