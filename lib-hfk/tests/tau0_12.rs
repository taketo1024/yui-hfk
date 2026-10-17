//! `G0(K)` for the prime knots with 12 crossings that are homologically thick, non-slice
//! and have `tau = 0`. By [SS] (Section 6), all of them have `G0(K) = {R(0,0)}`, the
//! `G0`-type of the unknot. Each check passes the genus, `g4 = 1` and `tau = 0`.
//!
//! This `G0`-type is self-dual, so the mirror has the same `G0`. For each knot, one side
//! has no candidate region to reject, so only `R(0,0)` must be accepted. `mirror` picks
//! that side (if both sides are so, the one with the smaller system for `R(0,0)`).
//!
//! The grid diagrams are KnotInfo's, of grid number 10, 11 or 12 (one test each).
//! All tests are slow:
//!
//!     cargo test --release --test tau0_12 -- --ignored --nocapture --test-threads=1

mod common;

use yui_hfk::G0CalcConfig;
use common::{run_all, run_tau};
use common::types::{K, O};

fn config() -> G0CalcConfig {
    G0CalcConfig { try_fast_ng: true, max_size: 2_000_000, ..Default::default() }
}

fn run_grid(label: &str, n: usize) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| k.grid_number() == n).collect();
    run_all(label, &selected, |k| k.name, |k| k.check(config()));
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_tau0_grid10() {
    run_grid("g0_tau0_grid10", 10);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_tau0_grid11() {
    run_grid("g0_tau0_grid11", 11);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_tau0_grid12() {
    run_grid("g0_tau0_grid12", 12);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn tau_all() {
    let knots: Vec<_> = KNOTS.iter().collect();
    run_tau("tau_all", &knots, config());
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "12n25", os: &[7,8,6,0,1,9,3,5,2,4], xs: &[9,5,3,4,7,2,0,1,6,8], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n475", os: &[6,0,4,2,3,5,1,7,8,9], xs: &[2,3,1,8,9,0,6,4,5,7], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n488", os: &[0,6,3,2,4,5,1,7,8,9], xs: &[4,2,1,8,9,0,6,3,5,7], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n835", os: &[9,8,7,6,0,3,1,4,2,5], xs: &[6,4,3,2,5,9,7,0,8,1], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n7", os: &[1,7,0,4,2,6,3,5,8,9,10], xs: &[5,2,3,1,9,10,7,0,4,6,8], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n10", os: &[8,9,5,6,0,1,10,7,3,4,2], xs: &[10,7,8,3,4,5,2,0,9,1,6], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n20", os: &[6,0,1,5,7,3,4,2,8,9,10], xs: &[3,4,9,0,2,10,1,6,5,7,8], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n41", os: &[0,1,6,7,10,4,2,3,5,9,8], xs: &[4,9,0,2,3,1,5,6,8,7,10], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n65", os: &[10,7,8,4,5,0,6,2,3,1,9], xs: &[8,9,6,7,3,4,1,10,0,5,2], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n198", os: &[10,7,8,6,1,0,2,5,4,3,9], xs: &[8,9,5,2,4,3,10,1,0,7,6], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n258", os: &[8,7,4,6,0,10,2,1,9,3,5], xs: &[10,9,8,2,5,3,7,6,4,0,1], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n265", os: &[8,5,0,10,2,1,3,7,4,9,6], xs: &[10,9,7,4,6,5,8,2,0,3,1], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n314", os: &[7,5,6,0,8,2,3,1,4,9,10], xs: &[2,9,3,4,1,10,0,5,7,6,8], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n322", os: &[1,3,0,2,7,4,10,6,5,9,8], xs: &[4,6,5,9,1,0,3,2,8,7,10], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n351", os: &[10,6,8,7,1,3,0,9,5,4,2], xs: &[8,9,5,3,4,10,6,2,1,0,7], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n352", os: &[8,9,5,2,10,0,7,1,4,3,6], xs: &[10,7,8,6,4,3,2,5,0,9,1], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n354", os: &[8,9,4,5,7,0,10,6,2,3,1], xs: &[10,6,8,2,3,4,1,0,7,9,5], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n362", os: &[8,9,5,6,10,1,0,7,3,4,2], xs: &[10,7,8,3,4,5,2,1,9,0,6], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n377", os: &[8,7,4,5,10,0,2,9,1,3,6], xs: &[10,9,8,2,3,6,7,4,5,0,1], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n442", os: &[9,1,0,2,6,4,5,3,8,7,10], xs: &[2,4,3,10,1,7,0,6,5,9,8], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n449", os: &[8,9,5,2,4,0,10,1,7,6,3], xs: &[10,7,8,6,1,3,2,5,4,0,9], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n483", os: &[1,9,2,0,10,5,7,4,3,6,8], xs: &[5,4,7,6,1,0,3,2,8,9,10], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n519", os: &[8,9,4,3,10,6,7,1,0,2,5], xs: &[10,6,8,7,5,1,2,4,3,9,0], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n552", os: &[8,9,5,6,10,0,7,2,4,3,1], xs: &[10,7,8,2,4,3,1,5,0,9,6], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n579", os: &[8,9,5,2,10,0,7,1,3,4,6], xs: &[10,7,8,6,4,3,2,5,9,0,1], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n630", os: &[8,9,4,2,0,10,1,6,7,3,5], xs: &[10,6,8,7,5,3,9,2,4,0,1], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n634", os: &[10,6,8,5,0,2,1,4,3,9,7], xs: &[8,9,4,1,3,10,7,0,6,5,2], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n650", os: &[10,9,8,2,0,1,6,7,4,3,5], xs: &[8,6,4,7,5,10,0,3,2,9,1], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n810", os: &[10,6,8,0,4,7,5,3,2,9,1], xs: &[8,9,3,5,10,2,1,0,6,4,7], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n812", os: &[10,6,8,3,1,2,0,9,7,4,5], xs: &[8,9,4,7,5,10,6,3,1,0,2], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n815", os: &[5,0,7,3,6,2,1,4,8,9,10], xs: &[2,4,1,9,0,10,5,7,3,6,8], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n868", os: &[8,9,7,6,2,4,0,1,10,5,3], xs: &[10,6,4,1,5,8,3,7,2,0,9], mirror: false, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n28", os: &[11,8,9,7,5,4,3,6,10,0,1,2], xs: &[9,10,6,3,11,1,8,2,5,4,7,0], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n31", os: &[9,10,6,4,2,0,11,8,1,3,5,7], xs: &[11,8,9,7,5,4,3,2,6,10,0,1], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n34", os: &[11,8,9,7,0,3,2,6,10,4,1,5], xs: &[9,10,6,4,5,11,8,1,3,0,7,2], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n118", os: &[11,8,9,4,1,2,0,3,5,7,10,6], xs: &[9,10,7,8,6,5,4,11,1,2,3,0], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n124", os: &[11,8,9,4,2,0,3,1,6,10,5,7], xs: &[9,10,6,7,5,4,11,8,2,3,0,1], mirror: true, g: 2, g4: 1, tau: 0, g0_type: O },
    K { name: "12n129", os: &[9,10,6,7,11,1,3,0,8,5,4,2], xs: &[11,8,9,3,5,4,6,2,1,10,0,7], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n132", os: &[0,6,4,10,2,3,1,11,7,5,8,9], xs: &[2,1,0,5,7,8,4,6,3,9,10,11], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n141", os: &[9,10,5,4,3,1,11,7,8,0,2,6], xs: &[11,7,9,8,6,4,2,1,3,5,10,0], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n210", os: &[9,10,6,4,2,5,11,8,0,1,7,3], xs: &[11,8,9,7,6,1,3,4,5,10,2,0], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n221", os: &[0,10,3,6,4,2,5,1,7,9,8,11], xs: &[7,4,8,2,1,11,0,6,3,5,10,9], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n264", os: &[3,0,2,4,8,7,11,6,1,5,10,9], xs: &[6,5,10,1,3,2,4,0,7,9,8,11], mirror: false, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n423", os: &[11,7,9,8,3,4,1,2,0,10,5,6], xs: &[9,10,5,4,7,11,6,8,3,1,0,2], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n437", os: &[0,4,2,1,11,3,5,7,8,6,10,9], xs: &[5,7,10,6,2,0,1,3,4,9,8,11], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n456", os: &[1,0,2,10,3,5,6,7,4,9,8,11], xs: &[5,4,7,1,0,11,2,3,8,6,10,9], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
    K { name: "12n681", os: &[1,0,2,7,4,6,5,3,8,9,10,11], xs: &[5,3,10,1,0,2,11,7,4,6,8,9], mirror: true, g: 4, g4: 1, tau: 0, g0_type: O },
    K { name: "12n690", os: &[11,8,9,4,1,3,0,2,6,10,7,5], xs: &[9,10,6,8,5,7,4,11,1,3,2,0], mirror: true, g: 3, g4: 1, tau: 0, g0_type: O },
];
