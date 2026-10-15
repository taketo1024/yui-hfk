//! Shared helpers for the `G0(K)` regression tests.

pub mod types;

use std::time::Instant;

use yui_hfk::g_zero::Region;
use yui_hfk::g_zero::calc::G0CalcError;
use yui_hfk::{G0Calc, G0CalcConfig, GridDiagram, tau};

use types::{K, Regions, regions_of, to_regions};

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
pub fn run_all<T>(label: &str, knots: &[&T], name: impl Fn(&T) -> &str, check: impl Fn(&T) -> usize) {
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

/// Checks `yui_hfk::tau` of the diagram tested against `tau` of each knot. As in
/// `check_g0`, resigning is "undecided" rather than wrong.
#[allow(dead_code)] // note: not every test crate uses it
pub fn run_tau(label: &str, knots: &[&K]) {
    assert!(!knots.is_empty(), "{label}: no knots selected");

    let undecided = knots.iter().filter(|k| match tau(&k.diagram()) {
        Ok(t) => {
            assert_eq!(t, k.tau, "{}: tau", k.name);
            false
        }
        Err(G0CalcError::Resigned) => {
            eprintln!("  {}: tau undecided", k.name);
            true
        }
    }).count();

    eprintln!("{label}: {} knots, {undecided} undecided", knots.len());
}
