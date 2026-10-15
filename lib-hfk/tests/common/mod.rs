//! Shared helpers for the `G0(K)` regression tests.

use std::collections::BTreeSet;
use std::time::Instant;

use yui_hfk::g_zero::{self, Corner, Region};
use yui_hfk::{G0Calc, G0CalcConfig, GridDiagram};

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
#[allow(dead_code)] // note: not every test crate uses it
pub fn dual(g0: &Regions) -> Regions {
    regions_of(&g_zero::dual(&to_regions(g0)))
}

/// Checks one knot, returning the number of regions left undecided. The known
/// `genus`, `g4` and `tau` are passed to `G0Calc::init`, and `expected` is passed
/// to `G0Calc::set_expected` as a hint for the order of the tests.
///
/// A large system can exhaust `max_size` and resign without being wrong, so
/// undecided regions are tolerated. What must hold is that nothing spurious
/// was accepted and nothing real was rejected.
pub fn check_g0(
    name: &str,
    gd: &GridDiagram,
    genus: usize,
    g4: Option<usize>,
    tau: Option<isize>,
    expected: &Regions,
    config: G0CalcConfig,
) -> usize {
    let mut calc = G0Calc::init(gd, Some(genus), g4, tau, config);
    calc.set_expected(to_regions(expected));
    calc.run();

    let accepted = regions_of(calc.accepted());
    let resigned = regions_of(calc.resigned());

    // note: an accepted region is realizable, but it is known to be minimal (in G0)
    // only when no resigned region lies inside it.
    let minimal: Vec<Region> = calc.accepted().iter()
        .filter(|a| !calc.resigned().iter().any(|r| a.includes(r)))
        .cloned()
        .collect();
    let spurious: Vec<_> = regions_of(&minimal).difference(expected).cloned().collect();
    assert!(spurious.is_empty(), "{name}: accepted non-realizable {spurious:?}");

    let missed: Vec<_> = expected.difference(&accepted)
        .filter(|r| !resigned.contains(*r))
        .collect();
    assert!(missed.is_empty(), "{name}: rejected realizable {missed:?}");

    resigned.len()
}

/// Runs `check` on each knot in turn, with one line per knot and a summary.
/// Set `RUST_LOG` (e.g. `yui_hfk=debug`) to log the calculation as well.
pub fn run_all<K>(label: &str, knots: &[&K], name: impl Fn(&K) -> &str, check: impl Fn(&K) -> usize) {
    assert!(!knots.is_empty(), "{label}: no knots selected");
    let _ = env_logger::builder().is_test(true).try_init();

    let undecided = knots.iter().filter(|&&k| {
        eprint!("  {}: ", name(k));
        let start = Instant::now();
        let n = check(k);
        let secs = start.elapsed().as_secs_f64();
        if n == 0 {
            eprintln!("agrees ({secs:.0}s)");
        } else {
            eprintln!("agrees, but {n} region(s) undecided ({secs:.0}s)");
        }
        n > 0
    }).count();

    eprintln!("{label}: {} knots, {undecided} with undecided regions", knots.len());
}
