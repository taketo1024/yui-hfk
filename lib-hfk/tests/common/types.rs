//! The data types of the `G0(K)` regression tests: `G0` as a set of regions, the
//! `G0`-types of [SS], and the knots.
#![allow(dead_code)] // note: each test crate uses only some of these

use std::collections::BTreeSet;

use yui_hfk::g_zero::{self, Corner, Region};
use yui_hfk::{G0CalcConfig, GridDiagram};

use super::check_g0;

/// `G0(K)` as a set of regions, each a set of corners.
pub type Regions = BTreeSet<BTreeSet<(isize, isize)>>;

pub fn regions_of(rs: &[Region]) -> Regions {
    rs.iter()
        .map(|r| r.corners().iter().map(|c| (c.i, c.j)).collect())
        .collect()
}

pub fn to_regions(g0: &Regions) -> Vec<Region> {
    g0.iter()
        .map(|r| Region::new(r.iter().map(|&(i, j)| Corner::new(i, j))))
        .collect()
}

/// `G0(K*)` from `G0(K)`, by `yui_hfk::g_zero::dual` ([SS] Cor. A.7).
pub fn dual(g0: &Regions) -> Regions {
    regions_of(&g_zero::dual(&to_regions(g0)))
}

/// A `G0`-type, as a list of regions, each a list of corners.
pub type G0Type = &'static [&'static [(isize, isize)]];

// The G0-types of Table 2 in [SS].
pub const O: G0Type = &[&[(0, 0)]];
pub const T2_3: G0Type = &[&[(0, 1)], &[(1, 0)]];
pub const T2_5: G0Type = &[&[(0, 2)], &[(1, 1)], &[(2, 0)]];
pub const T2_3_2_1: G0Type = &[&[(0, 2)], &[(0, 1), (1, 0)], &[(2, 0)]];
pub const T3_4: G0Type = &[&[(0, 3)], &[(1, 1)], &[(3, 0)]];
pub const T2_7: G0Type = &[&[(0, 3)], &[(1, 2)], &[(2, 1)], &[(3, 0)]];
pub const T3_5: G0Type = &[&[(0, 4)], &[(1, 2)], &[(2, 1)], &[(4, 0)]];
// note: not in Table 2. Computed in the 12-crossing tests for all seven knots of this type.
pub const T2_3_T3_5: G0Type = &[&[(0, 5)], &[(1, 3)], &[(2, 2)], &[(3, 1)], &[(5, 0)]];

pub struct K {
    pub name: &'static str,
    pub os: &'static [u8],
    pub xs: &'static [u8],
    /// The test uses the mirror of the diagram above.
    pub mirror: bool,
    pub g: usize,
    pub g4: usize,
    /// `tau` of the diagram tested.
    pub tau: isize,
    /// The `G0`-type of the side with `tau ≥ 0`.
    pub g0_type: G0Type,
}

impl K {
    pub fn grid_number(&self) -> usize {
        self.os.len()
    }

    /// The diagram the test uses.
    pub fn diagram(&self) -> GridDiagram {
        let gd = GridDiagram::from_sequences(self.os, self.xs);
        if self.mirror { gd.rotate90() } else { gd }
    }

    /// `G0` of the diagram tested: the dual type ([SS] Cor. A.7) where `tau < 0`.
    pub fn expected(&self) -> Regions {
        let g0 = self.g0_type.iter().map(|r| r.iter().copied().collect()).collect();
        if self.tau < 0 { dual(&g0) } else { g0 }
    }

    /// Runs `check_g0` on the diagram tested, passing `g`, `g4` and `tau`.
    pub fn check(&self, config: G0CalcConfig) -> usize {
        let (g4, tau) = (Some(self.g4), Some(self.tau));
        check_g0(self.name, &self.diagram(), self.g, g4, tau, &self.expected(), config)
    }
}
