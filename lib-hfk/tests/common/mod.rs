//! Shared helpers for the `G0(K)` regression tests.

use std::collections::BTreeSet;
use std::time::Instant;

use itertools::Itertools;
use yui_hfk::g_zero::{Corner, Region};
use yui_hfk::{G0Calc, G0CalcConfig, GridDiagram};

/// `G0(K)` as a set of regions, each a set of corners.
pub type Regions = BTreeSet<BTreeSet<(isize, isize)>>;

pub fn regions_of(rs: &[Region]) -> Regions {
    rs.iter()
        .map(|r| r.corners().iter().map(|c| (c.i, c.j)).collect())
        .collect()
}

/// `G0(K*)` from `G0(K)` by [SS] Cor. A.7: the minimal regions among
/// `cl({-p₁, …, -p_N})`, over all choices of a corner `pᵢ` of each region `Rᵢ`.
#[allow(dead_code)] // note: not every test crate uses it
pub fn dual(g0: &Regions) -> Regions {
    let candidates: Vec<Region> = g0.iter()
        .map(|r| r.iter().copied().collect_vec())
        .multi_cartesian_product()
        .map(|ps| Region::new(ps.into_iter().map(|(i, j)| Corner::new(-i, -j))))
        .collect();
    let minimal = candidates.iter()
        .filter(|r| !candidates.iter().any(|s| s != *r && r.includes(s)))
        .cloned()
        .collect_vec();
    regions_of(&minimal)
}

/// Checks one knot, returning the number of regions left undecided. The known
/// `genus`, `g4` and `tau` are passed to `G0Calc::init`.
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
