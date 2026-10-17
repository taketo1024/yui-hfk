//! `G0(K)` for the prime knots with up to 11 crossings that are homologically thick,
//! non-slice and have `tau = 0`: 14 of the 44 knots of Table 1 in [SS], of grid number
//! 8 to 11. All of them have `G0(K) = {R(0,0)}`, the `G0`-type `O` of the unknot.
//!
//! The grid diagrams are inlined, as in `thick_11.rs`. The type `O` is self-dual,
//! so the mirror has the same `G0`. Where `mirror` is set, the test uses the mirror of the
//! diagram recorded here, which has no candidate region to reject.
//!
//! With `G0CalcConfig::default()`, every region of every knot is decided.
//!
//! Cost is governed by the grid number (generator enumeration is O(n!)), so
//! only the smallest grids run by default. For the rest:
//!
//!     cargo test --release --test tau0_11 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::{run_all, run_tau};
use common::types::{K, O};

// note: tau comes from a separate search, so it is checked against the same knots.
fn run_tau_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    run_tau(label, &selected, G0CalcConfig::default());
}

#[test]
fn tau_small() {
    run_tau_over("tau_small", |k| k.grid_number() <= 9);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn tau_all() {
    run_tau_over("tau_all", |_| true);
}

fn check(k: &K) -> usize {
    k.check(G0CalcConfig::default())
}

fn run_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    run_all(label, &selected, |k| k.name, check);
}

#[test]
fn g0_small() {
    run_over("g0_small", |k| k.grid_number() <= 9);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_all() {
    run_over("g0_all", |_| true);
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "9_42", os: &[5,6,0,1,7,3,4,2], xs: &[7,3,4,5,2,0,1,6], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "10_136", os: &[8,4,5,0,6,2,3,1,7], xs: &[6,7,3,4,1,8,0,5,2], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "11n6", os: &[10,7,8,6,0,3,2,1,9,4,5], xs: &[8,9,4,2,5,10,7,6,3,0,1], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n20", os: &[7,8,4,5,0,6,9,3,1,2], xs: &[9,6,7,3,4,1,2,0,5,8], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "11n24", os: &[9,8,7,1,0,2,6,4,3,5], xs: &[7,6,4,5,3,9,1,0,8,2], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n34", os: &[2,0,1,4,7,5,6,3,8,9,10], xs: &[5,3,9,0,2,1,10,7,4,6,8], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n38", os: &[3,2,4,0,8,1,5,7,6], xs: &[7,5,1,3,2,6,0,4,8], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "11n45", os: &[8,6,4,3,5,0,1,10,7,9,2], xs: &[10,9,8,7,1,2,6,3,0,4,5], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n79", os: &[8,1,0,2,5,4,3,7,6,9], xs: &[2,4,3,9,1,0,6,5,8,7], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "11n92", os: &[9,8,7,1,2,0,6,3,4,5], xs: &[7,6,3,4,9,5,1,8,0,2], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n96", os: &[1,8,6,0,9,4,3,2,5,7], xs: &[4,3,2,5,1,0,6,7,8,9], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n138", os: &[9,8,7,6,0,2,1,4,3,5], xs: &[7,6,4,1,3,9,5,0,8,2], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "11n143", os: &[7,5,3,6,9,0,2,1,8,4], xs: &[9,8,7,1,4,5,6,3,2,0], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "11n145", os: &[9,8,7,6,0,1,5,2,3,4], xs: &[7,5,3,2,4,9,0,8,6,1], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
];
