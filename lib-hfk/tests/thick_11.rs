//! `G0(K)` for the prime knots with up to 11 crossings that are homologically thick,
//! non-slice and have `tau != 0`: 30 of the 44 knots of Table 1 in [SS], of grid number
//! 7 to 11. Those of `G0`-type `O` (`tau = 0`) are in `tau0_11.rs`.
//!
//! Both the expected values and the grid diagrams are inlined, so these tests
//! need no external data. Where `mirror` is set, the test uses the mirror of the
//! diagram recorded here, and `tau` refers to the mirror. This is so where the Swift
//! result came from a mirrored diagram.
//!
//! With `G0CalcConfig::default()`, every region of every knot is decided, except
//! for 11n31, whose remaining systems exceed `max_size` and are left undecided in
//! `g0_all`. `g0_11n31` verifies it fully with `max_size = 2_000_000` and `try_fast_ng`.
//!
//! Cost is governed by the grid number (generator enumeration is O(n!)), so
//! only the smallest grids run by default. For the rest:
//!
//!     cargo test --release --test thick_11 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::{run_all, run_tau};
use common::types::{K, T2_3, T2_3_2_1, T2_5, T3_4, T3_5};

fn check(k: &K) -> usize {
    k.check(G0CalcConfig::default())
}

fn run_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    run_all(label, &selected, |k| k.name, check);
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

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_11n31() {
    let k = KNOTS.iter().find(|k| k.name == "11n31").expect("no such knot");
    let config = G0CalcConfig { try_fast_ng: true, max_size: 2_000_000, ..Default::default() };
    run_all("g0_11n31", &[k], |k| k.name, |k| k.check(config));
}

// note: tau comes from a separate search, so it is checked against the same knots.
fn run_tau_over(label: &str, keep: impl Fn(&K) -> bool) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| keep(k)).collect();
    run_tau(label, &selected);
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
    K { name: "8_19", os: &[3,2,1,0,6,5,4], xs: &[6,5,4,3,2,1,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "10_124", os: &[4,3,2,1,0,7,6,5], xs: &[7,6,5,4,3,2,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "10_128", os: &[2,1,0,8,5,4,3,7,6], xs: &[7,4,3,2,1,0,6,5,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "10_132m", os: &[1,0,2,5,4,3,6,7,8], xs: &[4,3,7,1,0,8,2,5,6], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "10_139", os: &[6,4,2,1,0,8,7,5,3], xs: &[8,7,6,5,4,3,1,2,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "10_145", os: &[0,2,1,4,3,5,6,7,8], xs: &[3,7,5,0,8,2,1,4,6], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "10_152", os: &[7,5,3,2,9,1,8,0,6,4], xs: &[9,8,7,6,5,4,2,3,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "10_154", os: &[7,5,3,2,9,1,0,6,8,4], xs: &[9,8,7,6,4,5,2,1,3,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "10_161", os: &[6,4,2,1,8,0,7,5,3], xs: &[8,7,6,5,3,4,1,2,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n9", os: &[8,9,4,2,1,0,10,7,5,3,6], xs: &[10,7,8,6,5,4,3,1,2,9,0], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n12", os: &[9,8,7,1,2,0,4,6,3,5], xs: &[7,6,4,5,9,3,1,2,8,0], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n19", os: &[6,0,1,2,4,5,3,7,8], xs: &[2,3,7,8,0,1,6,4,5], mirror: true, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "11n27", os: &[7,5,2,1,0,6,9,8,4,3], xs: &[9,8,7,4,3,2,5,1,0,6], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n31", os: &[8,9,4,1,2,0,10,5,7,3,6], xs: &[10,7,8,6,5,4,3,1,2,9,0], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "11n57", os: &[7,8,4,2,1,0,9,6,5,3], xs: &[9,6,7,5,4,3,2,1,0,8], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n61", os: &[2,1,8,0,9,5,4,3,6,7], xs: &[5,4,3,6,2,1,0,7,8,9], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "11n70", os: &[9,8,7,6,1,0,2,5,4,3], xs: &[7,6,5,2,4,3,9,1,0,8], mirror: false, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "11n77", os: &[8,6,4,3,10,2,1,7,0,9,5], xs: &[10,9,8,7,6,5,3,2,4,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "11n80", os: &[10,7,8,6,0,3,1,2,9,4,5], xs: &[8,9,4,2,5,10,6,7,3,0,1], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n81", os: &[9,8,1,0,7,6,3,2,5,4], xs: &[6,3,5,4,2,1,9,8,0,7], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n88", os: &[7,5,4,2,0,9,8,1,3,6], xs: &[9,8,7,6,5,3,2,4,0,1], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n102", os: &[5,0,8,2,1,3,4,7,6,9], xs: &[1,4,3,5,6,0,9,2,8,7], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n104", os: &[7,8,2,1,9,0,6,5,4,3], xs: &[9,6,7,5,4,3,1,2,0,8], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n111", os: &[7,5,3,2,0,9,4,8,6,1], xs: &[9,8,7,6,5,1,0,2,3,4], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n126", os: &[7,5,3,9,2,1,0,8,6,4], xs: &[9,8,6,4,7,5,3,2,1,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "11n133", os: &[7,8,2,1,0,9,5,4,6,3], xs: &[9,5,6,4,3,2,1,7,0,8], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "11n135", os: &[7,8,2,9,0,6,5,4,3,1], xs: &[9,6,7,4,3,1,2,0,8,5], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "11n151", os: &[10,9,8,3,6,1,0,7,4,2,5], xs: &[8,6,4,7,10,5,3,2,1,9,0], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n152", os: &[10,9,8,3,6,5,0,7,4,2,1], xs: &[8,6,4,7,10,2,3,1,0,9,5], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "11n183", os: &[7,5,2,3,1,9,0,8,6,4], xs: &[9,8,6,7,4,2,5,3,1,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
];
