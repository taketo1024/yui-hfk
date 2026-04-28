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
