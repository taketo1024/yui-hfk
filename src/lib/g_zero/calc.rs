use yui_core::lc::Lc;
use yui_core::num::FF2;

use crate::{GridComplex, GridDiagram, HFK};
use crate::grid_complex::{GCFiltered, GCInfGen};
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
        todo!()
    }

    pub fn accepted(&self) -> &[Region] {
        &self.candidates.accepted
    }

    // Tests whether `r` is realizable: tries to find a cycle homologous to `z`
    // that lies entirely inside `r`.
    fn process(&mut self, r: Region) {
        todo!()
    }
}
