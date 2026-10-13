//! Regression fixture for `G0(K)`, checked against results computed by the
//! original Swift implementation (which uses the correct `Alg = -a_0`).
//!
//! Both the expected values and the grid diagrams are inlined, so these tests
//! need no external data. The Swift results came from a different diagram
//! source, so some knots are mirrored relative to the diagram recorded here;
//! `mirror` says which, and `tau` refers to the oriented diagram.
//!
//! Cost is governed by the grid number (generator enumeration is O(n!)), so
//! only the smallest grids run by default. For the rest:
//!
//!     cargo test --release --test oracle -- --ignored --nocapture

use std::collections::BTreeSet;

use yui_hfk::g_zero::Region;
use yui_hfk::g_zero::calc::G0CalcError;
use yui_hfk::{G0Calc, G0CalcConfig, GridDiagram};

struct K {
    name: &'static str,
    os: &'static [u8],
    xs: &'static [u8],
    /// The Swift result describes the mirror of the diagram above.
    mirror: bool,
    genus: usize,
    tau: isize,
    /// `G0(K)` as a set of regions, each a set of corners.
    regions: &'static [&'static [(isize, isize)]],
}

impl K {
    fn grid_number(&self) -> usize {
        self.os.len()
    }

    /// The diagram in the orientation the expected values refer to.
    fn diagram(&self) -> GridDiagram {
        let gd = GridDiagram::from_sequences(self.os, self.xs);
        if self.mirror { gd.rotate90() } else { gd }
    }

    fn expected(&self) -> Regions {
        self.regions.iter().map(|r| r.iter().copied().collect()).collect()
    }
}

type Regions = BTreeSet<BTreeSet<(isize, isize)>>;

fn regions_of(rs: &[Region]) -> Regions {
    rs.iter()
        .map(|r| r.corners().iter().map(|c| (c.i, c.j)).collect())
        .collect()
}

/// Checks one knot, returning the number of regions left undecided.
///
/// A large system can exhaust `max_size` and resign without being wrong, so
/// undecided regions are tolerated. What must hold is that nothing spurious
/// was accepted and nothing real was rejected.
fn check(k: &K) -> usize {
    let mut calc = G0Calc::from_gd(&k.diagram(), None, None, G0CalcConfig::default());
    calc.run();

    let accepted = regions_of(calc.accepted());
    let resigned = regions_of(calc.resigned());
    let expected = k.expected();

    assert_eq!(calc.genus(), k.genus, "{}: genus", k.name);

    let spurious: Vec<_> = accepted.difference(&expected).collect();
    assert!(spurious.is_empty(), "{}: accepted non-realizable {spurious:?}", k.name);

    let missed: Vec<_> = expected.difference(&accepted)
        .filter(|r| !resigned.contains(*r))
        .collect();
    assert!(missed.is_empty(), "{}: rejected realizable {missed:?}", k.name);

    resigned.len()
}

fn run_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    assert!(!selected.is_empty(), "{label}: no knots selected");

    let mut undecided = 0;
    for k in &selected {
        let n = check(k);
        if n > 0 {
            undecided += 1;
            eprintln!("  {}: agrees, but {n} region(s) undecided", k.name);
        }
    }

    eprintln!("{label}: {} knots, {undecided} with undecided regions", selected.len());
}

#[test]
fn g0_small_grids() {
    run_over("g0_small_grids", |k| k.grid_number() <= 8);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_all() {
    run_over("g0_all", |_| true);
}

/// tau comes from a separate search, so check it against the same fixture.
/// As in `check`, resigning is "undecided" rather than wrong.
fn run_tau_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    assert!(!selected.is_empty(), "{label}: no knots selected");

    let mut undecided = 0;
    for k in &selected {
        match yui_hfk::tau(&k.diagram()) {
            Ok(t) => assert_eq!(t, k.tau, "{}: tau", k.name),
            Err(G0CalcError::Resigned) => {
                undecided += 1;
                eprintln!("  {}: tau undecided", k.name);
            }
        }
    }

    eprintln!("{label}: {} knots, {undecided} undecided", selected.len());
}

#[test]
fn tau_small_grids() {
    run_tau_over("tau_small_grids", |k| k.grid_number() <= 8);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn tau_all() {
    run_tau_over("tau_all", |_| true);
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "8_19", os: &[3,2,1,0,6,5,4], xs: &[6,5,4,3,2,1,0], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "10_124", os: &[4,3,2,1,0,7,6,5], xs: &[7,6,5,4,3,2,1,0], mirror: false, genus: 4, tau: 4, regions: &[&[(0,4)], &[(1,2)], &[(2,1)], &[(4,0)]] },
    K { name: "9_42", os: &[5,6,0,1,7,3,4,2], xs: &[7,3,4,5,2,0,1,6], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "10_128", os: &[2,1,0,8,5,4,3,7,6], xs: &[7,4,3,2,1,0,6,5,8], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "10_132m", os: &[1,0,2,5,4,3,6,7,8], xs: &[4,3,7,1,0,8,2,5,6], mirror: true, genus: 2, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "10_136", os: &[8,4,5,0,6,2,3,1,7], xs: &[6,7,3,4,1,8,0,5,2], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "10_139", os: &[6,4,2,1,0,8,7,5,3], xs: &[8,7,6,5,4,3,1,2,0], mirror: false, genus: 4, tau: 4, regions: &[&[(0,4)], &[(1,2)], &[(2,1)], &[(4,0)]] },
    K { name: "10_145", os: &[0,2,1,4,3,5,6,7,8], xs: &[3,7,5,0,8,2,1,4,6], mirror: true, genus: 2, tau: 2, regions: &[&[(0,2)], &[(1,1)], &[(2,0)]] },
    K { name: "10_161", os: &[6,4,2,1,8,0,7,5,3], xs: &[8,7,6,5,3,4,1,2,0], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n19", os: &[6,0,1,2,4,5,3,7,8], xs: &[2,3,7,8,0,1,6,4,5], mirror: true, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n38", os: &[3,2,4,0,8,1,5,7,6], xs: &[7,5,1,3,2,6,0,4,8], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "10_152", os: &[7,5,3,2,9,1,8,0,6,4], xs: &[9,8,7,6,5,4,2,3,1,0], mirror: false, genus: 4, tau: 4, regions: &[&[(0,4)], &[(1,2)], &[(2,1)], &[(4,0)]] },
    K { name: "10_154", os: &[7,5,3,2,9,1,0,6,8,4], xs: &[9,8,7,6,4,5,2,1,3,0], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n102", os: &[5,0,8,2,1,3,4,7,6,9], xs: &[1,4,3,5,6,0,9,2,8,7], mirror: true, genus: 2, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n104", os: &[7,8,2,1,9,0,6,5,4,3], xs: &[9,6,7,5,4,3,1,2,0,8], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n107", os: &[2,8,1,0,9,5,4,3,6,7], xs: &[5,4,6,3,2,1,0,7,8,9], mirror: false, genus: 4, tau: 2, regions: &[&[(0,2)], &[(1,1)], &[(2,0)]] },
    K { name: "11n111", os: &[7,5,3,2,0,9,4,8,6,1], xs: &[9,8,7,6,5,1,0,2,3,4], mirror: false, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n12", os: &[9,8,7,1,2,0,4,6,3,5], xs: &[7,6,4,5,9,3,1,2,8,0], mirror: false, genus: 2, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n126", os: &[7,5,3,9,2,1,0,8,6,4], xs: &[9,8,6,4,7,5,3,2,1,0], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n133", os: &[7,8,2,1,0,9,5,4,6,3], xs: &[9,5,6,4,3,2,1,7,0,8], mirror: false, genus: 4, tau: 2, regions: &[&[(0,1), (1,0)], &[(0,2)], &[(2,0)]] },
    K { name: "11n135", os: &[7,8,2,9,0,6,5,4,3,1], xs: &[9,6,7,4,3,1,2,0,8,5], mirror: false, genus: 3, tau: 2, regions: &[&[(0,2)], &[(1,1)], &[(2,0)]] },
    K { name: "11n138", os: &[9,8,7,6,0,2,1,4,3,5], xs: &[7,6,4,1,3,9,5,0,8,2], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n143", os: &[7,5,3,6,9,0,2,1,8,4], xs: &[9,8,7,1,4,5,6,3,2,0], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n145", os: &[9,8,7,6,0,1,5,2,3,4], xs: &[7,5,3,2,4,9,0,8,6,1], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n183", os: &[7,5,2,3,1,9,0,8,6,4], xs: &[9,8,6,7,4,2,5,3,1,0], mirror: false, genus: 3, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n20", os: &[7,8,4,5,0,6,9,3,1,2], xs: &[9,6,7,3,4,1,2,0,5,8], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n24", os: &[9,8,7,1,0,2,6,4,3,5], xs: &[7,6,4,5,3,9,1,0,8,2], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n27", os: &[7,5,2,1,0,6,9,8,4,3], xs: &[9,8,7,4,3,2,5,1,0,6], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n57", os: &[7,8,4,2,1,0,9,6,5,3], xs: &[9,6,7,5,4,3,2,1,0,8], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n61", os: &[2,1,8,0,9,5,4,3,6,7], xs: &[5,4,3,6,2,1,0,7,8,9], mirror: false, genus: 4, tau: 2, regions: &[&[(0,1), (1,0)], &[(0,2)], &[(2,0)]] },
    K { name: "11n70", os: &[9,8,7,6,1,0,2,5,4,3], xs: &[7,6,5,2,4,3,9,1,0,8], mirror: false, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n79", os: &[8,1,0,2,5,4,3,7,6,9], xs: &[2,4,3,9,1,0,6,5,8,7], mirror: false, genus: 2, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n81", os: &[9,8,1,0,7,6,3,2,5,4], xs: &[6,3,5,4,2,1,9,8,0,7], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n88", os: &[7,5,4,2,0,9,8,1,3,6], xs: &[9,8,7,6,5,3,2,4,0,1], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
    K { name: "11n92", os: &[9,8,7,1,2,0,6,3,4,5], xs: &[7,6,3,4,9,5,1,8,0,2], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n96", os: &[1,8,6,0,9,4,3,2,5,7], xs: &[4,3,2,5,1,0,6,7,8,9], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n151", os: &[10,9,8,3,6,1,0,7,4,2,5], xs: &[8,6,4,7,10,5,3,2,1,9,0], mirror: false, genus: 2, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n152", os: &[10,9,8,3,6,5,0,7,4,2,1], xs: &[8,6,4,7,10,2,3,1,0,9,5], mirror: false, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n31", os: &[8,9,4,1,2,0,10,5,7,3,6], xs: &[10,7,8,6,5,4,3,1,2,9,0], mirror: false, genus: 3, tau: 2, regions: &[&[(0,2)], &[(1,1)], &[(2,0)]] },
    K { name: "11n34", os: &[2,0,1,4,7,5,6,3,8,9,10], xs: &[5,3,9,0,2,1,10,7,4,6,8], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n45", os: &[8,6,4,3,5,0,1,10,7,9,2], xs: &[10,9,8,7,1,2,6,3,0,4,5], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n46", os: &[2,0,3,1,9,7,10,5,4,6,8], xs: &[7,5,6,4,2,0,3,1,8,9,10], mirror: false, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n47", os: &[8,6,4,2,0,3,10,9,5,7,1], xs: &[10,9,8,5,4,7,6,1,0,2,3], mirror: false, genus: 4, tau: 2, regions: &[&[(0,2)], &[(1,1)], &[(2,0)]] },
    K { name: "11n6", os: &[10,7,8,6,0,3,2,1,9,4,5], xs: &[8,9,4,2,5,10,7,6,3,0,1], mirror: false, genus: 3, tau: 0, regions: &[&[(0,0)]] },
    K { name: "11n77", os: &[8,6,4,3,10,2,1,7,0,9,5], xs: &[10,9,8,7,6,5,3,2,4,1,0], mirror: false, genus: 4, tau: 4, regions: &[&[(0,4)], &[(1,2)], &[(2,1)], &[(4,0)]] },
    K { name: "11n80", os: &[10,7,8,6,0,3,1,2,9,4,5], xs: &[8,9,4,2,5,10,6,7,3,0,1], mirror: true, genus: 3, tau: 1, regions: &[&[(0,1)], &[(1,0)]] },
    K { name: "11n9", os: &[8,9,4,2,1,0,10,7,5,3,6], xs: &[10,7,8,6,5,4,3,1,2,9,0], mirror: false, genus: 4, tau: 3, regions: &[&[(0,3)], &[(1,1)], &[(3,0)]] },
];
