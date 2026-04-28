use yui_core::lc::Lc;
use yui_core::num::FF2;

use crate::{GridComplex, GridDiagram, HFK};
use crate::grid_complex::{GCFiltered, GCInfGen};
use super::{Candidates, Region};

pub struct G0Calc {
    cx: GridComplex<GCFiltered>,
    z: Lc<GCInfGen, FF2>,
    candidates: Candidates,
}

impl G0Calc {
    pub fn new(gd: GridDiagram, hfk: &HFK) -> Self {
        todo!()
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
