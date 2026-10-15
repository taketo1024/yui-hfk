//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, whose linear systems exceed `max_size = 2_000_000`, against the `G0`-types
//! of Tables 1 and 2 in [SS] (v2). The others are in `thick_12_t*.rs`, by `G0`-type.
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! There is one test per knot, with `try_fast_ng` and the `max_size` given there, with
//! which every region is decided (12n96 and 12n203: not decided yet). The `ignore` reason
//! gives the run time, one knot at a time on an Apple M1 Ultra (20 cores). "No hint" marks
//! times measured without the expected type as a hint (12n693 takes 19m without it).
//! Each test is slow:
//!
//!     cargo test --release --test thick_12_heavy g0_12n59 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T2_3, T2_5, T2_7, T3_4};

fn run(name: &str, max_size: usize) {
    let k = KNOTS.iter().find(|k| k.name == name).expect("no such knot");
    let config = G0CalcConfig { try_fast_ng: true, max_size, ..Default::default() };
    run_all(name, &[k], |k| k.name, |k| k.check(config));
}

#[test]
#[ignore = "slow: 49m (no hint)"]
fn g0_12n45() {
    run("12n45", 5_000_000);
}

#[test]
#[ignore = "slow: 1h48m (no hint)"]
fn g0_12n382() {
    run("12n382", 5_000_000);
}

#[test]
#[ignore = "slow: 1h25m (no hint)"]
fn g0_12n549() {
    run("12n549", 5_000_000);
}

#[test]
#[ignore = "slow: 1h41m (no hint)"]
fn g0_12n220() {
    run("12n220", 5_000_000);
}

#[test]
#[ignore = "slow: 22m (no hint)"]
fn g0_12n228() {
    run("12n228", 5_000_000);
}

#[test]
#[ignore = "slow: 9m"]
fn g0_12n693() {
    run("12n693", 5_000_000);
}

#[test]
#[ignore = "slow: 5h45m (no hint)"]
fn g0_12n696() {
    run("12n696", 7_000_000);
}

#[test]
#[ignore = "slow: 11h55m"]
fn g0_12n59() {
    run("12n59", 14_000_000);
}

#[test]
#[ignore = "slow: not decided yet"]
fn g0_12n96() {
    run("12n96", 17_000_000);
}

#[test]
#[ignore = "slow: not decided yet"]
fn g0_12n203() {
    run("12n203", 17_000_000);
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "12n45", os: &[10,9,8,2,0,1,5,3,7,4,6], xs: &[8,7,5,6,3,4,2,10,0,9,1], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n382", os: &[8,7,4,0,10,1,2,9,3,5,6], xs: &[10,9,8,6,3,7,5,4,0,1,2], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n549", os: &[8,7,5,3,10,1,2,0,9,4,6], xs: &[10,9,8,6,4,5,7,3,1,0,2], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n220", os: &[9,10,5,3,2,1,0,11,6,8,4,7], xs: &[11,8,9,7,6,5,4,2,0,3,10,1], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n228", os: &[9,10,6,4,2,1,11,0,8,5,7,3], xs: &[11,8,9,7,6,5,3,4,1,0,2,10], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n693", os: &[9,10,5,3,11,1,8,0,7,6,4,2], xs: &[11,8,9,6,4,5,3,2,1,0,10,7], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n696", os: &[9,10,5,3,1,0,11,6,8,4,2,7], xs: &[11,8,9,7,6,5,4,2,3,1,10,0], mirror: true, g: 4, g4: 3, tau: -3, g0_type: T3_4 },
    K { name: "12n59", os: &[9,10,5,3,2,1,11,0,6,8,4,7], xs: &[11,8,9,7,6,5,4,2,1,3,10,0], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n96", os: &[1,0,7,10,6,2,5,3,4,9,8,11], xs: &[5,3,4,2,1,11,0,7,8,6,10,9], mirror: true, g: 3, g4: 3, tau: -3, g0_type: T2_7 },
    K { name: "12n203", os: &[1,0,7,10,6,2,5,3,4,9,8,11], xs: &[5,4,2,3,1,11,0,7,8,6,10,9], mirror: true, g: 3, g4: 3, tau: -3, g0_type: T2_7 },
];
