use log::info;
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
        info!("compute HFK ...");
        let hfk = HFK::from(&gd);
        info!("HFK:\n{}", hfk.display_table("M", "A"));

        let g = hfk.genus() as usize;
        info!("g = {g}");

        let candidates = Self::make_candidates(&hfk, g);
        info!("candidate regions: {}", candidates.regions.len());
        for r in &candidates.regions {
            info!("  {r}");
        }

        info!("prepare generators ...");
        let gc = Self::build_complex(gd, &candidates);
        info!("  #generators: {}", gc.gens.len());

        info!("compute homological generator");
        let z = gc.hom_gen();

        G0Calc { gc, z, candidates }
    }

    fn make_candidates(hfk: &HFK, g: usize) -> Candidates { 
        let mut candidates = Candidates::generate_from(g);
        candidates.deduplicate();
        candidates.filter_by_hfk(&hfk);
        candidates.sort();
        candidates
    }

    fn build_complex(gd: GridDiagram, candidates: &Candidates) -> GridComplex<GCFiltered> {
        let s_max = candidates.regions.iter()
            .map(|r| r.right_end().max(0))
            .max()
            .unwrap_or(0);
        let s_range = -2 * s_max ..= 1;

        GridComplex::new_filtered(gd, |x| {
            if !s_range.contains(&x.maslov_deg) { return false; }
            let s = -(x.maslov_deg.div_euclid(2)); // >= 0
            candidates.regions.iter()
                .filter(|r| r.right_end() == s)
                .any(|r| !r.contains((s, x.alex_deg + s))) // equivalent to: !r.shift(-s, -s).contains(0, x.alex_deg)
        })
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
        info!("try: {r} (remain: {})", self.candidates.regions.len());

        let s = r.right_end().max(0);
        let r_shifted = r.shift(-s, -s);
        let z_shifted = self.shift_z(s as usize);
        match self.gc.find_homologous(&z_shifted, &r_shifted) {
            Some(_) => { 
                info!("OK: {r}"); 
                self.candidates.accept(&r); 
            },
            None => { 
                info!("NG: {r}"); 
                self.candidates.reject(&r); 
            }
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
