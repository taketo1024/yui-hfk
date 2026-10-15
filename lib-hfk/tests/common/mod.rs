//! Shared helpers for the `G0(K)` regression tests.
#![allow(dead_code)] // note: each test crate uses only some of these

pub mod types;

use std::time::Instant;

use itertools::Itertools;
use yui_hfk::g_zero::{Candidates, Region, dual};
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

/// Checks `G0(K)` = the expected type `T` of `k` through the mirror `K*` ([SS] Cor. A.7), when
/// the side of `k` leaves some regions undecided. On `K` (with `config`), `G0(K) = Min(A ∪ V)`
/// for the accepted regions `A` and the realizable ones `V` among the undecided ones `U`, and
/// `T` must be one of these. Then the tests of `mirror_tests` on `K*` (with `mirror_config`)
/// exclude all the others.
pub fn check_via_mirror(k: &K, config: G0CalcConfig, mirror_config: G0CalcConfig) {
    init_logger();

    let t = to_regions(&k.expected());

    // note: a block, so that the calculation on K is dropped before the one on K*.
    let (accept, reject) = {
        let mut calc = G0Calc::init(&k.diagram(), Some(k.g), Some(k.g4), Some(k.tau), config);
        calc.set_expected(t.clone());
        calc.run();

        // note: no region inside a rejected one is realizable.
        let in_rejected = |r: &Region| calc.rejected().iter().any(|s| s.includes(r));
        let missed: Vec<_> = t.iter().filter(|e| in_rejected(e)).collect();
        assert!(missed.is_empty(), "{}: rejected realizable {missed:?}", k.name);

        let u: Vec<_> = calc.resigned().iter().filter(|r| !in_rejected(r)).cloned().collect();
        mirror_tests(calc.accepted(), &u, &t)
    };
    if accept.is_empty() && reject.is_empty() {
        return;
    }

    let mut candidates = Candidates::new();
    accept.iter().chain(&reject).for_each(|r| candidates.add(r.clone()));

    let mut calc = G0Calc::init_with_candidates(&k.diagram().rotate90(), candidates, mirror_config);
    calc.set_expected(accept.clone());
    calc.run();

    let missed: Vec<_> = accept.iter().filter(|q| !calc.accepted().contains(q)).collect();
    assert!(missed.is_empty(), "{}: not accepted on the mirror {missed:?}", k.name);

    let kept: Vec<_> = reject.iter()
        .filter(|p| calc.accepted().contains(p) || calc.resigned().contains(p))
        .collect();
    assert!(kept.is_empty(), "{}: not rejected on the mirror {kept:?}", k.name);
}

/// Regions to accept and to reject on `K*`, chosen greedily, that exclude every alternative
/// `G0(K) = Min(A ∪ V) ≠ T` with `V ⊆ U` closed under transposition. `T` must be one of these.
///
/// With `D = dual(T)` and `D'` the dual of an alternative ([SS] Cor. A.7), `D ≠ D'`. So either
/// some region of `D` contains no region of `D'` (to accept: realizable only under `D`), or
/// some region of `D'` contains no region of `D` (to reject: realizable only under `D'`).
pub fn mirror_tests(a: &[Region], u: &[Region], t: &[Region]) -> (Vec<Region>, Vec<Region>) {
    assert!(u.len() <= 16, "too many undecided regions: {}", u.len());

    let subset = |mask: usize| -> Vec<Region> {
        u.iter().enumerate().filter(|(i, _)| mask >> i & 1 == 1).map(|(_, r)| r.clone()).collect()
    };
    let possible: Vec<Vec<Region>> = (0..1usize << u.len())
        .map(subset)
        .filter(|v| v.iter().all(|r| v.contains(&r.transposed())))
        .map(|v| min_regions(&[a, &v].concat()))
        .filter(|g| !g.is_empty())
        .collect();
    assert!(possible.iter().any(|g| regions_of(g) == regions_of(t)), "T is not a possible G0(K)");

    let d = dual(t);
    let alternatives: Vec<Vec<Region>> = possible.iter()
        .filter(|g| regions_of(g) != regions_of(t))
        .map(|g| dual(g))
        .collect();

    // note: a test (r, accept) excludes D' iff the realizability of r under D' differs from accept.
    let realizable = |r: &Region, g0: &[Region]| g0.iter().any(|p| r.includes(p));
    let pool: Vec<(Region, bool)> = d.iter().map(|q| (q.clone(), true))
        .chain(alternatives.iter().flatten().filter(|p| !realizable(p, &d)).unique().map(|p| (p.clone(), false)))
        .collect();
    let excludes = |(r, accept): &(Region, bool), i: usize| realizable(r, &alternatives[i]) != *accept;

    let mut remaining: Vec<usize> = (0..alternatives.len()).collect();
    let (mut accept, mut reject) = (vec![], vec![]);
    while !remaining.is_empty() {
        let test = pool.iter()
            .max_by_key(|test| (remaining.iter().filter(|&&i| excludes(test, i)).count(), test.1))
            .expect("D is nonempty");
        assert!(remaining.iter().any(|&i| excludes(test, i)), "no test excludes the rest");
        remaining.retain(|&i| !excludes(test, i));
        if test.1 {
            accept.push(test.0.clone());
        } else {
            reject.push(test.0.clone());
        }
    }
    (accept, reject)
}

fn init_logger() {
    let _ = env_logger::builder().is_test(true).try_init();
}

fn min_regions(rs: &[Region]) -> Vec<Region> {
    rs.iter()
        .filter(|r| !rs.iter().any(|s| s != *r && r.includes(s)))
        .cloned()
        .collect()
}

/// Runs `check` on each knot in turn, with one line per knot and a summary.
/// Set `RUST_LOG` (e.g. `yui_hfk=debug`) to log the calculation as well.
pub fn run_all<T>(label: &str, knots: &[&T], name: impl Fn(&T) -> &str, check: impl Fn(&T) -> usize) {
    assert!(!knots.is_empty(), "{label}: no knots selected");
    init_logger();

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
